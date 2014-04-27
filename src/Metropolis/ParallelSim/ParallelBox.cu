/*
	Implements methods related to managing data between the host and device.
	Subclass of Box.
*/

#include "ParallelBox.cuh"
#include "ParallelCalcs.cuh"
#include "Metropolis/SerialSim/SerialCalcs.h"

using namespace std;

ParallelBox::ParallelBox(): Box()
{
	nextMol = 0;
	numChanged = 0;
}

ParallelBox::~ParallelBox()
{
	// TODO: free device memory
}

int ParallelBox::chooseMolecules(const int N)
{
	if (changedIndices == NULL || changedMols == NULL)
	{
		changedIndices = new int[N];
		changedMols = new Molecule[N];
	}
	
	if (N > 0)
	{
		int i = 0, j;
		
		//if the size needs to change
		if (numChanged != N)
		{
			delete[] changedIndices;
			delete[] changedMols;
			changedIndices = new int[N];
			changedMols = new Molecule[N];
			numChanged = N;
		}
		
		//if the next molecule has already been chosen
		if (nextMol != 0)
		{
			changedIndices[i++] = nextMol;
			nextMol = 0;
		}
		
		for (i = i; i < N; i++)
		{
			nextMol = chooseMolecule();
			
			//make sure this molecule is not already in the batch
			for (j = 0; j < i; j++)
			{
				//if this molecule is already in the batch,
				//stop here and start with this molecule next time
				if (changedIndices[j] == nextMol)
				{
					numChanged = i;
					return numChanged;
				}
			}
			changedIndices[i] = nextMol;
		}
		//if batch had no duplicates, reset nextMol to 0
		nextMol = 0;
		return N;
	}
	else
	{
		return 0;
	}
}

void ParallelBox::changeMolecules()
{
	for (nextChangeIdx = 0; nextChangeIdx < numChanged; nextChangeIdx++)
	{
		changeMolecule(changedIndices[nextChangeIdx]);
	}
}

void ParallelBox::saveChangedMol(int molIdx)
{
	Molecule *sourceMol = &molecules[molIdx];

	//free memory of changedMol before allocate memory
	deleteMolMemberArrays(changedMols + nextChangeIdx);

	memcpy(changedMols + nextChangeIdx, sourceMol, sizeof(Molecule));

	createMolMemberArrays(changedMols + nextChangeIdx, sourceMol);

	copyMolecule(changedMols + nextChangeIdx, sourceMol);
}

void ParallelBox::createMolMemberArrays(Molecule *mol, Molecule *sourceMol)
{
	mol->atoms = new Atom[sourceMol->numOfAtoms];
	mol->bonds = new Bond[sourceMol->numOfBonds];
	mol->angles = new Angle[sourceMol->numOfAngles];
	mol->dihedrals = new Dihedral[sourceMol->numOfDihedrals];
	mol->hops = new Hop[sourceMol->numOfHops];
}

void ParallelBox::deleteMolMemberArrays(Molecule *mol)
{
	delete[] mol->atoms;
	delete[] mol->bonds;
	delete[] mol->angles;
	delete[] mol->dihedrals;
	delete[] mol->hops;
}

void ParallelBox::toggleChange(int changeIdx)
{
	//create temporary copy of backup Molecule
	Molecule *backupMol = (Molecule*) malloc(sizeof(Molecule));
	createMolMemberArrays(backupMol, changedMols + changeIdx);
	copyMolecule(backupMol, changedMols + changeIdx);
	
	//save working copy
	nextChangeIdx = changeIdx;
	saveChangedMol(changedIndices[changeIdx]);
	nextChangeIdx = 0;
	
	//overwrite working copy with backup
	copyMolecule(molecules + changedIndices[changeIdx], backupMol);
	
	//de-allocate temporary copy of backup
	deleteMolMemberArrays(backupMol);
}

int ParallelBox::changeMolecule(int molIdx)
{
	Box::changeMolecule(molIdx);
	writeChangeToDevice(molIdx);
	
	return molIdx;
}

int ParallelBox::rollback(int molIdx)
{
	Box::rollback(molIdx);
	writeChangeToDevice(molIdx);
	
	return molIdx;
}

void ParallelBox::copyDataToDevice()
{
	//create AtomData on host, and fill atomic data arrays on device
	atomsH = new AtomData(atoms, atomCount);
	cudaMalloc(&xD, atomCount * sizeof(Real));
	cudaMalloc(&yD, atomCount * sizeof(Real));
	cudaMalloc(&zD, atomCount * sizeof(Real));
	cudaMalloc(&sigmaD, atomCount * sizeof(Real));
	cudaMalloc(&epsilonD, atomCount * sizeof(Real));
	cudaMalloc(&chargeD, atomCount * sizeof(Real));
	cudaMemcpy(xD, atomsH->x, atomCount * sizeof(Real), cudaMemcpyHostToDevice);
	cudaMemcpy(yD, atomsH->y, atomCount * sizeof(Real), cudaMemcpyHostToDevice);
	cudaMemcpy(zD, atomsH->z, atomCount * sizeof(Real), cudaMemcpyHostToDevice);
	cudaMemcpy(sigmaD, atomsH->sigma, atomCount * sizeof(Real), cudaMemcpyHostToDevice);
	cudaMemcpy(epsilonD, atomsH->epsilon, atomCount * sizeof(Real), cudaMemcpyHostToDevice);
	cudaMemcpy(chargeD, atomsH->charge, atomCount * sizeof(Real), cudaMemcpyHostToDevice);
	
	//create device AtomData struct with pointers to filled-in atomic data arrays
	AtomData *tempAD = (AtomData*) malloc(sizeof(AtomData));
	tempAD->x = xD;
	tempAD->y = yD;
	tempAD->z = zD;
	tempAD->sigma = sigmaD;
	tempAD->epsilon = epsilonD;
	tempAD->charge = chargeD;
	tempAD->atomCount = atomsH->atomCount;
	cudaMalloc(&atomsD, sizeof(AtomData));
	cudaMemcpy(atomsD, tempAD, sizeof(AtomData), cudaMemcpyHostToDevice);
	
	//create MoleculeData on host, and fill molecular data arrays on device
	moleculesH = new MoleculeData(molecules, moleculeCount);
	cudaMalloc(&atomsIdxD, moleculeCount * sizeof(int));
	cudaMalloc(&numOfAtomsD, moleculeCount * sizeof(int));
	cudaMemcpy(atomsIdxD, moleculesH->atomsIdx, moleculeCount * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(numOfAtomsD, moleculesH->numOfAtoms, moleculeCount * sizeof(int), cudaMemcpyHostToDevice);
	
	//create device MoleculeData struct with pointers to filled-in molecular data arrays
	MoleculeData *tempMD = (MoleculeData*) malloc(sizeof(MoleculeData));
	tempMD->atomsIdx = atomsIdxD;
	tempMD->numOfAtoms = numOfAtomsD;
	tempMD->moleculeCount = moleculesH->moleculeCount;
	cudaMalloc(&moleculesD, sizeof(MoleculeData));
	cudaMemcpy(moleculesD, tempMD, sizeof(MoleculeData), cudaMemcpyHostToDevice);
	
	/***********************************************************************
		Without Parallel Steps
	************************************************************************
	
	//data structures for neighbor batch in energy calculation
	nbrMolsH = (int*) malloc(moleculeCount * sizeof(int));
	molBatchH = (int*) malloc(moleculeCount * sizeof(int));
	cudaMalloc(&(nbrMolsD), moleculeCount * sizeof(int));
	cudaMalloc(&(molBatchD), moleculeCount * sizeof(int));
	
	//upper bound on number of atoms in any molecule
	maxMolSize = 0;
	for (int i = 0; i < moleculesH->moleculeCount; i++)
	{
		if (moleculesH->numOfAtoms[i] > maxMolSize)
		{
			maxMolSize = moleculesH->numOfAtoms[i];
		}
	}
	
	//energies array on device has one segment for each molecule
	//where each segment has the maximum number of
	//possible interatomic energies for one pair of molecules
	energyCount = moleculesH->moleculeCount * maxMolSize * maxMolSize;
	cudaMalloc(&(energiesD), energyCount * sizeof(Real));
	
	************************************************************************
		Without Parallel Steps
	***********************************************************************
		With Parallel Steps
	***********************************************************************/
	
	//data structures for neighbor batch in energy calculation
	nbrMolsH = (int*) malloc(ParallelCalcs::MAX_PAR_STEPS * moleculeCount * sizeof(int));
	molBatchH = (int*) malloc(ParallelCalcs::MAX_PAR_STEPS * moleculeCount * sizeof(int));
	cudaMalloc(&(nbrMolsD), ParallelCalcs::MAX_PAR_STEPS * moleculeCount * sizeof(int));
	cudaMalloc(&(molBatchD), ParallelCalcs::MAX_PAR_STEPS * moleculeCount * sizeof(int));
	
	//upper bound on number of atoms in any molecule
	maxMolSize = 0;
	for (int i = 0; i < moleculesH->moleculeCount; i++)
	{
		if (moleculesH->numOfAtoms[i] > maxMolSize)
		{
			maxMolSize = moleculesH->numOfAtoms[i];
		}
	}
	
	//energies array on device has one segment for each molecule
	//where each segment has the maximum number of
	//possible interatomic energies for one pair of molecules
	energyCount = ParallelCalcs::MAX_PAR_STEPS * moleculesH->moleculeCount * maxMolSize * maxMolSize;
	cudaMalloc(&(energiesD), energyCount * sizeof(Real));
	
	/***********************************************************************
		With Parallel Steps
	***********************************************************************/
	
	//initialize energies to 0
	cudaMemset(energiesD, 0, energyCount * sizeof(Real));
	
	//copy Environment to device
	cudaMalloc(&(environmentD), sizeof(Environment));
	cudaMemcpy(environmentD, environment, sizeof(Environment), cudaMemcpyHostToDevice);
}

void ParallelBox::writeChangeToDevice(int changeIdx)
{
	//update AtomData atomsH (MoleculeData will not change)
	int startIdx = moleculesH->atomsIdx[changeIdx];
	for (int i = 0; i < molecules[changeIdx].numOfAtoms; i++)
	{
		atomsH->x[startIdx + i] = molecules[changeIdx].atoms[i].x;
		atomsH->y[startIdx + i] = molecules[changeIdx].atoms[i].y;
		atomsH->z[startIdx + i] = molecules[changeIdx].atoms[i].z;
		//sigma, epsilon, and charge will not change, so there is no need to update those arrays
	}

	//copy changed atom data to device
	cudaMemcpy(xD + startIdx, atomsH->x + startIdx, moleculesH->numOfAtoms[changeIdx] * sizeof(Real), cudaMemcpyHostToDevice);
	cudaMemcpy(yD + startIdx, atomsH->y + startIdx, moleculesH->numOfAtoms[changeIdx] * sizeof(Real), cudaMemcpyHostToDevice);
	cudaMemcpy(zD + startIdx, atomsH->z + startIdx, moleculesH->numOfAtoms[changeIdx] * sizeof(Real), cudaMemcpyHostToDevice);
	//sigma, epsilon, and charge will not change, so there is no need to update those arrays
}

bool ParallelBox::changedMolsWithinCutoff(int mol1, int mol2)
{
	Atom atom1 = molecules[changedIndices[mol1]].atoms[environment->primaryAtomIndex];
	Atom atom2 = molecules[changedIndices[mol2]].atoms[environment->primaryAtomIndex];
		
	//calculate periodic difference in coordinates
	Real deltaX = SerialCalcs::makePeriodic(atom1.x - atom2.x, environment->x);
	Real deltaY = SerialCalcs::makePeriodic(atom1.y - atom2.y, environment->y);
	Real deltaZ = SerialCalcs::makePeriodic(atom1.z - atom2.z, environment->z);
	
	Real r2 = (deltaX * deltaX) +
				(deltaY * deltaY) + 
				(deltaZ * deltaZ);
	
	return r2 < environment->cutoff * environment->cutoff;
}