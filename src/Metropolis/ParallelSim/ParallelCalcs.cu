/*
	Contains the methods required to calculate energies in parallel.
*/

#include "ParallelCalcs.h"
#include "ParallelCalcs.cuh"
#include "ParallelBox.cuh"
#include <string>
#include "Metropolis/Utilities/FileUtilities.h"
#include "Metropolis/Box.h"
#include "Metropolis/SimulationArgs.h"
#include "Metropolis/SerialSim/SerialCalcs.h"

#define NO -1
#define MAX_WARP 32
#define MOL_BLOCK 256
#define BATCH_BLOCK 512
#define AGG_BLOCK 512

using namespace std;

Box* ParallelCalcs::createBox(string inputPath, InputFileType inputType, long* startStep, long* steps)
{
	ParallelBox* box = new ParallelBox();
	if (!loadBoxData(inputPath, inputType, box, startStep, steps))
	{
		if (inputType != InputFile::Unknown)
		{
			std::cerr << "Error: Could not build from file: " << inputPath << std::endl;
			return NULL;
		}
		else
		{
			std::cerr << "Error: Can not build environment with unknown file: " << inputPath << std::endl;
			return NULL;
		}
	}
	box->copyDataToDevice();
	return (Box*) box;
}

void ParallelCalcs::runParallelSteps(int simSteps, Box *box, Real &systemEnergy, int &accepted, int &rejected)
{
	ParallelBox *pBox = (ParallelBox*) box;
	
	if (pBox != NULL)
	{
		Molecule *molecules = pBox->getMolecules();
		Environment *enviro = pBox->getEnvironment();
		
		Real *oldEnergyConts, *newEnergyConts, wrongE, rightE;
		Real  kT = ParallelCalcs::kBoltz * enviro->temp;
		bool accept;
		
		//calculate old energy
		if (systemEnergy == 0)
		{
			systemEnergy = ParallelCalcs::calcSystemEnergy(box);
		}
		
		int move = 0, stepsTaken;
		while (move < simSteps)
		{
			if (move + ParallelCalcs::MAX_PAR_STEPS > simSteps)
			{
				stepsTaken = pBox->chooseMolecules(simSteps - move);
			}
			else
			{
				stepsTaken = pBox->chooseMolecules(ParallelCalcs::MAX_PAR_STEPS);
			}
			
			oldEnergyConts = ParallelCalcs::calcMolecularEnergyContributions(pBox);
			
			//TODO: pBox->changeMolecules should clear and then change the
			//molecules currently specified by pBox->changedIndices.
			pBox->changeMolecules();
			
			newEnergyConts = ParallelCalcs::calcMolecularEnergyContributions(pBox);
			
			for (int mol1 = 0; mol1 < stepsTaken; mol1++)
			{
				accept = false;
				if (newEnergyConts[mol1] < oldEnergyConts[mol1])
				{
					accept = true;
				}
				else
				{
					if(exp(-(newEnergyConts[mol1] - oldEnergyConts[mol1])) / kT >= randomReal(0.0, 1.0))
					{
						accept = true;
					}
				}
				
				if (accept)
				{
					accepted++;
					systemEnergy += newEnergyConts[mol1] - oldEnergyConts[mol1];
					for (int mol2 = mol1 + 1; mol2 < stepsTaken; mol2++)
					{
						if (pBox->changedMolsWithinCutoff(mol1, mol2))
						{
							pBox->toggleChange(mol1);//undo mol1 change
							pBox->toggleChange(mol2);//undo mol2 change
							wrongE = SerialCalcs::calcInterMolecularEnergy(pBox->molecules, pBox->changedIndices[mol1], pBox->changedIndices[mol2], pBox->environment);
							pBox->toggleChange(mol1);//redo mol1 change
							rightE = SerialCalcs::calcInterMolecularEnergy(pBox->molecules, pBox->changedIndices[mol1], pBox->changedIndices[mol2], pBox->environment);
							oldEnergyConts[mol2] += rightE - wrongE;
							pBox->toggleChange(mol2);//redo mol2 change
						}
					}
				}
				else
				{
					rejected++;
					for (int mol2 = mol1 + 1; mol2 < stepsTaken; mol2++)
					{
						if (pBox->changedMolsWithinCutoff(mol1, mol2))
						{
							wrongE = SerialCalcs::calcInterMolecularEnergy(pBox->molecules, pBox->changedIndices[mol1], pBox->changedIndices[mol2], pBox->environment);
							pBox->toggleChange(mol1);//undo mol1 change
							rightE = SerialCalcs::calcInterMolecularEnergy(pBox->molecules, pBox->changedIndices[mol1], pBox->changedIndices[mol2], pBox->environment);
							newEnergyConts[mol2] += rightE - wrongE;
							pBox->toggleChange(mol1);//redo mol1 change
						}
					}
					pBox->toggleChange(mol1);//permanently undo mol1 change for rejection
				}
			}
			move += stepsTaken;
		}
	}
}

Real ParallelCalcs::calcSystemEnergy(Box *box)
{
	Real totalEnergy = 0;
	
	//for each molecule
	for (int mol = 0; mol < box->moleculeCount; mol++)
	{
		//use startIdx parameter to prevent double-calculating energies (Ex mols 3->5 and mols 5->3)
		totalEnergy += calcMolecularEnergyContribution(box, mol, mol + 1);
	}

    return totalEnergy;
}

Real* ParallelCalcs::calcMolecularEnergyContributions(ParallelBox *pBox)
{
	int i, N = pBox->numChanged;
	cudaStream_t streams[N];
	
	for (i = 0; i < N; i++)
	{
		cudaStreamCreate(&(streams[i]));
	}
	
	//initialize neighbor molecule slots to NO
	cudaMemset(pBox->nbrMolsD, NO, pBox->moleculeCount * sizeof(int));
	
	//check molecule distances in parallel, conditionally filling in pBox->nbrMolsD
	for (i = 0; i < N; i++)
	{
		checkMoleculeDistances<<<pBox->moleculeCount / MOL_BLOCK + 1, MOL_BLOCK, 0, streams[i]>>>(pBox->moleculesD, pBox->atomsD, pBox->changedIndices[i], 0, pBox->environmentD, pBox->nbrMolsD + pBox->moleculeCount * i);
	}
	
	//get sparse neighbor list
	cudaMemcpy(pBox->nbrMolsH, pBox->nbrMolsD, N * pBox->moleculeCount * sizeof(int), cudaMemcpyDeviceToHost);
	
	//innitialize neighbor batch to NO
	memset(pBox->molBatchH, NO, N * pBox->moleculeCount * sizeof(int));
	
	//sparse list compaction to create batch of valid neighbor molecules to send to GPU
	//NOTE: It is possible that this can be done more efficiently on the GPU, without ever involving the CPU.
	//      We tested the thrust library's implementation (included with Cuda), but it was slower than using the CPU.
	//      Perhaps look into writing our own Parallel Stream Compaction kernel? (Looking at you, future group!)
	int batchSizes[N], j;
	for (i = 0; i < N; i++)
	{
		batchSizes[i] = 0;
		for (j = 0; j < pBox->moleculeCount; j++)
		{
			if (pBox->nbrMolsH[j + i * pBox->moleculeCount] != NO)
			{
				pBox->molBatchH[batchSizes[i]++ + i * pBox->moleculeCount] = j;
			}
		}
	}
	
	//initialize first energy to 0
	//All other energies will have been reset from the previous aggregation run.
	//Technically, this should already be 0, but this is a safeguard against implementation changes.
	cudaMemset(pBox->energiesD, 0, sizeof(Real));
	
	//copy neighbor batch to device
	cudaMemcpy(pBox->molBatchD, pBox->molBatchH, N * pBox->moleculeCount * sizeof(int), cudaMemcpyHostToDevice);

	//There will only be as many energy segments filled in as there are molecules in the batch.
	int validEnergies[N], segmentSize = pBox->moleculeCount * pBox->maxMolSize * pBox->maxMolSize;
	Real *contributions = (Real*) malloc(N * sizeof(Real));
	
	for (i = 0; i < N; i++)
	{
		validEnergies[i] = batchSizes[i] * pBox->maxMolSize * pBox->maxMolSize;
				
		//calculate interatomic energies between changed molecule and all molecules in batch
		calcInterAtomicEnergy<<<validEnergies[i] / BATCH_BLOCK + 1, BATCH_BLOCK, 0, streams[i]>>>
		(pBox->moleculesD, pBox->atomsD, pBox->changedIndices[i], pBox->environmentD, pBox->energiesD + i * segmentSize, validEnergies[i], pBox->molBatchD + i * pBox->moleculeCount, pBox->maxMolSize);
	}
	
	for (i = 0; i < N; i++)
	{
		//a batch size of 3 seems to offer the best tradeoff between parallelization and locality
		int batchSize = 3, blockSize = AGG_BLOCK;
		int numBaseThreads = validEnergies[i] / (batchSize);
		for (int i = 1; i < validEnergies[i]; i *= batchSize)
		{
			//there is no need for giant blocks when only a few threads are being run
			if (blockSize > MAX_WARP && numBaseThreads / i + 1 < blockSize)
			{
				blockSize /= 2;
			}
			
			//do the next aggregation pass
			aggregateEnergies<<<numBaseThreads / (i * blockSize) + 1, blockSize, 0, streams[i]>>>
			(pBox->energiesD + i * segmentSize, validEnergies[i], i, batchSize);
		}
	}
	
	for (i = 0; i < N; i++)
	{
		//the total energy will be stored in the first index of pBox->energiesD
		cudaMemcpy(&(contributions[i]), pBox->energiesD + i * segmentSize, sizeof(Real), cudaMemcpyDeviceToHost);
		
		//reset final energy to 0, resulting in a clear energies array for the next calculation
		cudaMemset(pBox->energiesD + i * segmentSize, 0, sizeof(Real));
	}
	
	return contributions;
}

Real ParallelCalcs::calcMolecularEnergyContribution(Box *box, int molIdx, int startIdx)
{
	ParallelBox *pBox = (ParallelBox*) box;
	
	if (pBox == NULL)
	{
		return 0;
	}
	
	return calcBatchEnergy(pBox, createMolBatch(pBox, molIdx, startIdx), molIdx);
}

Real ParallelCalcs::calcBatchEnergy(ParallelBox *box, int numMols, int molIdx)
{
	if (numMols > 0)
	{
		//There will only be as many energy segments filled in as there are molecules in the batch.
		int validEnergies = numMols * box->maxMolSize * box->maxMolSize;
		
		//initialize first energy to 0
		//All other energies will have been reset from the previous aggregation run.
		//Technically, this should already be 0, but this is a safeguard against implementation changes.
		cudaMemset(box->energiesD, 0, sizeof(Real));
		
		//copy neighbor batch to device
		cudaMemcpy(box->molBatchD, box->molBatchH, box->moleculeCount * sizeof(int), cudaMemcpyHostToDevice);
		
		//calculate interatomic energies between changed molecule and all molecules in batch
		calcInterAtomicEnergy<<<validEnergies / BATCH_BLOCK + 1, BATCH_BLOCK>>>
		(box->moleculesD, box->atomsD, molIdx, box->environmentD, box->energiesD, validEnergies, box->molBatchD, box->maxMolSize);
		
		return getEnergyFromDevice(box, validEnergies);
	}
	else
	{
		return 0;
	}
}

int ParallelCalcs::createMolBatch(ParallelBox *box, int currentMol, int startIdx)
{
	//initialize neighbor molecule slots to NO
	cudaMemset(box->nbrMolsD, NO, box->moleculeCount * sizeof(int));
	
	//check molecule distances in parallel, conditionally filling in box->nbrMolsD
	checkMoleculeDistances<<<box->moleculeCount / MOL_BLOCK + 1, MOL_BLOCK>>>(box->moleculesD, box->atomsD, currentMol, startIdx, box->environmentD, box->nbrMolsD);
	
	//get sparse neighbor list
	cudaMemcpy(box->nbrMolsH, box->nbrMolsD, box->moleculeCount * sizeof(int), cudaMemcpyDeviceToHost);
	
	//innitialize neighbor batch to NO
	memset(box->molBatchH, NO, box->moleculeCount * sizeof(int));
	
	//sparse list compaction to create batch of valid neighbor molecules to send to GPU
	//NOTE: It is possible that this can be done more efficiently on the GPU, without ever involving the CPU.
	//      We tested the thrust library's implementation (included with Cuda), but it was slower than using the CPU.
	//      Perhaps look into writing our own Parallel Stream Compaction kernel? (Looking at you, future group!)
	int batchSize = 0;
	for (int i = startIdx; i < box->moleculeCount; i++)
	{
		if (box->nbrMolsH[i] != NO)
		{
			box->molBatchH[batchSize++] = i;
		}
	}
	
	return batchSize;
}

Real ParallelCalcs::getEnergyFromDevice(ParallelBox *box, int validEnergies)
{
	Real totalEnergy = 0;
	
	//a batch size of 3 seems to offer the best tradeoff between parallelization and locality
	int batchSize = 3, blockSize = AGG_BLOCK;
	int numBaseThreads = validEnergies / (batchSize);
	for (int i = 1; i < validEnergies; i *= batchSize)
	{
		//there is no need for giant blocks when only a few threads are being run
		if (blockSize > MAX_WARP && numBaseThreads / i + 1 < blockSize)
		{
			blockSize /= 2;
		}
		
		//do the next aggregation pass
		aggregateEnergies<<<numBaseThreads / (i * blockSize) + 1, blockSize>>>
		(box->energiesD, validEnergies, i, batchSize);
	}
	
	//the total energy will be stored in the first index of box->energiesD
	cudaMemcpy(&totalEnergy, box->energiesD, sizeof(Real), cudaMemcpyDeviceToHost);
	
	//reset final energy to 0, resulting in a clear energies array for the next calculation
	cudaMemset(box->energiesD, 0, sizeof(Real));
	
	return totalEnergy;
}

__global__ void ParallelCalcs::checkMoleculeDistances(MoleculeData *molecules, AtomData *atoms, int currentMol, int startIdx, Environment *enviro, int *inCutoff)
{
	int otherMol = blockIdx.x * blockDim.x + threadIdx.x;
	
	//checks validity of molecule pair
	if (otherMol < molecules->moleculeCount && otherMol >= startIdx && otherMol != currentMol)
	{
		//find primary atom indices for this pair of molecules
		int atom1 = molecules->atomsIdx[currentMol] + enviro->primaryAtomIndex;
		int atom2 = molecules->atomsIdx[otherMol] + enviro->primaryAtomIndex;
			
		//calculate periodic difference in coordinates
		Real deltaX = makePeriodic(atoms->x[atom1] - atoms->x[atom2], enviro->x);
		Real deltaY = makePeriodic(atoms->y[atom1] - atoms->y[atom2], enviro->y);
		Real deltaZ = makePeriodic(atoms->z[atom1] - atoms->z[atom2], enviro->z);
		
		Real r2 = (deltaX * deltaX) +
					(deltaY * deltaY) + 
					(deltaZ * deltaZ);
		
		//if within curoff, write index to inCutoff
		if (r2 < enviro->cutoff * enviro->cutoff)
		{
			inCutoff[otherMol] = otherMol;
		}
	}
}

__global__ void ParallelCalcs::calcInterAtomicEnergy(MoleculeData *molecules, AtomData *atoms, int currentMol, Environment *enviro, Real *energies, int energyCount, int *molBatch, int maxMolSize)
{
	int energyIdx = blockIdx.x * blockDim.x + threadIdx.x, segmentSize = maxMolSize * maxMolSize;
	
	//check validity of thread
	if (energyIdx < energyCount and molBatch[energyIdx / segmentSize] != NO)
	{
		//get other molecule index
		int otherMol = molBatch[energyIdx / segmentSize];
		
		//get atom pair for this thread
		int x = (energyIdx % segmentSize) / maxMolSize, y = (energyIdx % segmentSize) % maxMolSize;
		
		//check validity of atom pair
		if (x < molecules->numOfAtoms[currentMol] && y < molecules->numOfAtoms[otherMol])
		{
			//get atom indices
			int atom1 = molecules->atomsIdx[currentMol] + x;
			int atom2 = molecules->atomsIdx[otherMol] + y;
			
			//check validity of atoms (ensure they are not dummy atoms)
			if (atoms->sigma[atom1] >= 0 && atoms->epsilon[atom1] >= 0 && atoms->sigma[atom2] >= 0 && atoms->epsilon[atom2] >= 0)
			{
				Real totalEnergy = 0;
			  
				//calculate periodic distance between atoms
				Real deltaX = makePeriodic(atoms->x[atom1] - atoms->x[atom2], enviro->x);
				Real deltaY = makePeriodic(atoms->y[atom1] - atoms->y[atom2], enviro->y);
				Real deltaZ = makePeriodic(atoms->z[atom1] - atoms->z[atom2], enviro->z);
				
				Real r2 = (deltaX * deltaX) +
					 (deltaY * deltaY) + 
					 (deltaZ * deltaZ);
				
				//calculate interatomic energies
				totalEnergy += calc_lj(atoms, atom1, atom2, r2);
				totalEnergy += calcCharge(atoms->charge[atom1], atoms->charge[atom2], sqrt(r2));
				
				//store energy
				energies[energyIdx] = totalEnergy;
			}
		}
	}
}

__global__ void ParallelCalcs::aggregateEnergies(Real *energies, int energyCount, int interval, int batchSize)
{
	//calculate starting index for this thread
	int idx = batchSize * interval * (blockIdx.x * blockDim.x + threadIdx.x), i;
	
	for (i = 1; i < batchSize; i++)
	{
		//make strides of size 'interval'
		if (idx + i * interval < energyCount)
		{
			//add energy, and then reset it
			energies[idx] += energies[idx + i * interval];
			energies[idx + i * interval] = 0;
		}
	}
}

__device__ Real ParallelCalcs::calc_lj(AtomData *atoms, int atom1, int atom2, Real r2)
{
    //store LJ constants locally
    Real sigma = calcBlending(atoms->sigma[atom1], atoms->sigma[atom2]);
    Real epsilon = calcBlending(atoms->epsilon[atom1], atoms->epsilon[atom2]);
    
    if (r2 == 0.0)
    {
        return 0.0;
    }
    else
    {
    	//calculate terms
    	const Real sig2OverR2 = (sigma*sigma) / r2;
		const Real sig6OverR6 = (sig2OverR2*sig2OverR2*sig2OverR2);
    	const Real sig12OverR12 = (sig6OverR6*sig6OverR6);
    	const Real energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);
        return energy;
    }
}

__device__ Real ParallelCalcs::calcCharge(Real charge1, Real charge2, Real r)
{  
    if (r == 0.0)
    {
        return 0.0;
    }
    else
    {
    	// conversion factor below for units in kcal/mol
    	const Real e = 332.06;
        return (charge1 * charge2 * e) / r;
    }
}

__device__ Real ParallelCalcs::makePeriodic(Real x, Real boxDim)
{
    
    while(x < -0.5 * boxDim)
    {
        x += boxDim;
    }

    while(x > 0.5 * boxDim)
    {
        x -= boxDim;
    }

    return x;

}

__device__ Real ParallelCalcs::calcBlending(Real d1, Real d2)
{
    return sqrt(d1 * d2);
}

__device__ int ParallelCalcs::getXFromIndex(int idx)
{
    int c = -2 * idx;
    int discriminant = 1 - 4 * c;
    int qv = (-1 + sqrtf(discriminant)) / 2;
    return qv + 1;
}

__device__ int ParallelCalcs::getYFromIndex(int x, int idx)
{
    return idx - (x * x - x) / 2;
}