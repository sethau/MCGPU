/*
	Driver for the simulation. Takes in a SimulationArgs object and creates the
	the necessary Box type, state file output path, etc.

	Author: Nathan Coleman
	Created: February 21, 2014
	
	-> February 26, by Albert Wallace
	-> March 28, by Joshua Mosby
	-> April 21, by Nathan Coleman
*/

#include <string>
#include <iostream>
#include <fstream>
#include <time.h>
#include <stdio.h>

#include "Simulation.h"
#include "SimulationArgs.h"
#include "Box.h"
#include "Metropolis/Utilities/MathLibrary.h"
#include "Metropolis/Utilities/Parsing.h"
#include "SerialSim/SerialBox.h"
#include "SerialSim/SerialCalcs.h"
#include "ParallelSim/ParallelCalcs.h"
#include "Utilities/FileUtilities.h"

#define RESULTS_FILE_DEFAULT "run"
#define RESULTS_FILE_EXT ".results"
#define PAR_STEPS 2

//Constructor & Destructor
Simulation::Simulation(SimulationArgs simArgs)
{
	args = simArgs;

	stepStart = 0;

	if (simArgs.simulationMode == SimulationMode::Parallel)
		box = ParallelCalcs::createBox(args.filePath, args.fileType, &stepStart, &simSteps);
	else
		box = SerialCalcs::createBox(args.filePath, args.fileType, &stepStart, &simSteps);

	if (box == NULL)
	{
		std::cerr << "Error: Unable to initialize simulation Box" << std::endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		std::cout << "Using seed: " << box->environment->randomseed << std::endl;
		seed(box->environment->randomseed);
	}

	if (args.stepCount > 0)
		simSteps = args.stepCount;
}

Simulation::~Simulation()
{
	if(box != NULL)
	{
		delete box;
		box = NULL;
	}
}

/*
systemE = calcSystemE();
accepted = 0;
rejected = 0;
For each Nth simulation step:
changeMols = chooseRandomMols(N);
preChangeEs = calculateContributions(changeMols);
changeRandomMols(changeMols);
postChangeEs = calculateContributions(changeMols);
for M1 in changeMols:
	if accept(M1):
		accepted++;
		systemE += postChangeEs[M1] – preChangeEs[M1];
		for M2 after M1 in changeMols:
			if distance(M1, M2) < cutoff: 
				rollback(M1);
				rollback(M2) ;
				wrongE = interE(M1, M2);
				redoChange(M1);
				rightE = interE(M1, M2);
				preChangeEs[M2] += rightE – wrongE;
				redoChange(M2);
	else:
		rejected++;
		for M2 after M1 in changeMols:
			if distance(M1, M2) < cutoff:
				wrongE = interE(M1, M2);
				rollback(M1);
				rightE = interE(M1, M2);
				postChangeEs[M2] += rightE – wrongE;
				redoChange(M1);
		rollback(M1);
clearMoleculeBackups();

*/
void Simulation::runV2()
{
	std::cout << "Simulation Name: " << args.simulationName << std::endl;
	
	Molecule *molecules = box->getMolecules();
	Environment *enviro = box->getEnvironment();
	
	Real systemEnergy = 0;
	Real *newEnergyConts, *oldEnergyConts, wrongE, rightE;
	Real  kT = kBoltz * enviro->temp;
	int accepted = 0, rejected = 0;
	bool accept;

	clock_t startTime, endTime;
    startTime = clock();
	
	//calculate old energy
	if (systemEnergy == 0)
	{
		systemEnergy = ParallelCalcs::calcSystemEnergy(box);
	}
	
	std::cout << std::endl << "Running " << simSteps << " steps" << std::endl << std::endl;
	
	int move = stepStart;
	while (move < stepStart + simSteps)
	{
		if (args.statusInterval > 0 && (move - stepStart) % args.statusInterval == 0)
		{
			std::cout << "Step " << move << ":\n--Current Energy: " << systemEnergy << std::endl;	
		}
		
		if (args.stateInterval > 0 && move > stepStart && (move - stepStart) % args.stateInterval == 0)
		{
			std::cout << std::endl;
			saveState(move);
			std::cout << std::endl;
		}
		
		if (move + PAR_STEPS > stepStart + simSteps)
		{
			box->chooseMolecules(stepStart + simSteps - move);
		}
		else
		{
			box->chooseMolecules(PAR_STEPS);
		}
		
		oldEnergyConts = ParallelCalcs::calcMolecularEnergyContributions(box);
		
		//TODO: box->changeMolecules should clear and then re
		box->changeMolecules();
		
		newEnergyConts = ParallelCalcs::calcMolecularEnergyContributions(box);
		
		for (int mol1 = 0; mol1 < PAR_STEPS; mol1++)
		{
			accept = false;
			if (newEnergyConts[mol1] < oldEnergyConts[mol1])
			{
				accept = true;
			}
			else
			{
				if(exp(-(newEnergyConts[mol1] - oldEnergyConts[mol1]) / kT >= randomReal(0.0, 1.0))
				{
					accept = true;
				}
			}
			
			if (accept)
			{
				accepted++;
				systemEnergy += newEnergyConts[mol1] - oldEnergyConts[mol1];
				for (int mol2 = mol1 + 1; mol2 < PAR_STEPS; mol2++)
				{
					if (box->changedMolsWithinCutoff(mol1, mol2))
					{
						//TODO: ParallelBox::toggleChange() should treat its argument as an index in changedIndices.
						//It should either rollback a change and store the new version,
						//or redo the change and store the old version each time it is called,
						//depending on what is currently in the backup for that molecule
						box->toggleChange(mol1);//undo mol1 change
						box->toggleChange(mol2);//undo mol2 change
						wrongE = SerialCalcs::calcInterMolecularEnergy(box->molecules, box->changedIndices[mol1], box->changedIndices[mol2], box->environment);
						box->toggleChange(mol1);//redo mol1 change
						rightE = SerialCalcs::calcInterMolecularEnergy(box->molecules, box->changedIndices[mol1], box->changedIndices[mol2], box->environment);
						oldEnergyConts[mol2] += rightE - wrongE;
						box->toggleChange(mol2);//redo mol2 change
					}
				}
			}
			else
			{
				rejected++;
				for (int mol2 = mol1 + 1; mol2 < PAR_STEPS; mol2++)
				{
					if (box->changedMolsWithinCutoff(mol1, mol2))
					{
						wrongE = SerialCalcs::calcInterMolecularEnergy(box->molecules, box->changedIndices[mol1], box->changedIndices[mol2], box->environment);
						box->toggleChange(mol1);//undo mol1 change
						rightE = SerialCalcs::calcInterMolecularEnergy(box->molecules, box->changedIndices[mol1], box->changedIndices[mol2], box->environment);
						newEnergyConts[mol2] += rightE - wrongE;
						box->toggleChange(mol1);//redo mol1 change
					}
				}
				box->toggleChange(mol1);//permanently undo mol1 change for rejection
			}
		}
		
		
		
		
		if(accept)
		{
			accepted++;
			systemEnergy += newEnergyCont - oldEnergyCont;
		}
		else
		{
			rejected++;
			//restore previous configuration
			box->rollback(changeIdx);
		}
		move += PAR_STEPS;
	}

	endTime = clock();
    double diffTime = difftime(endTime, startTime) / CLOCKS_PER_SEC;

	std::cout << "Step " << (stepStart + simSteps) << ":\r\n--Current Energy: " << systemEnergy << std::endl;
	systemEnergy = systemEnergy;

	// Save the final state of the simulation
	if (args.stateInterval >= 0)
	{
		saveState((stepStart + simSteps));
	}
	
	std::cout << std::endl << "Finished running " << simSteps << " steps" << std::endl;
	std::cout << "Final Energy: " << systemEnergy << std::endl;
	std::cout << "Run Time: " << diffTime << " seconds" << std::endl;
	std::cout << "Accepted Moves: " << accepted << std::endl;
	std::cout << "Rejected Moves: " << rejected << std::endl;
	std::cout << "Acceptance Ratio: " << 100.0 * accepted / (accepted + rejected) << '\%' << std::endl;
	
	saveResults(systemEnergy, diffTime, accepted, rejected);
}

//Utility
void Simulation::run()
{
	std::cout << "Simulation Name: " << args.simulationName << std::endl;
	//declare variables common to both parallel and serial
	Molecule *molecules = box->getMolecules();
	Environment *enviro = box->getEnvironment();
	
	Real systemEnergy = 0;
	Real newEnergyCont, oldEnergyCont;
	Real  kT = kBoltz * enviro->temp;
	int accepted = 0, rejected = 0;

	clock_t startTime, endTime;
    startTime = clock();
	
	//calculate old energy
	if (systemEnergy == 0)
	{
		if (args.simulationMode == SimulationMode::Parallel)
		{
			systemEnergy = ParallelCalcs::calcSystemEnergy(box);
		}
		else
		{
			systemEnergy = SerialCalcs::calcSystemEnergy(molecules, enviro);
		}
	}
	
	std::cout << std::endl << "Running " << simSteps << " steps" << std::endl << std::endl;
	
	for(int move = stepStart; move < (stepStart + simSteps); move++)
	{
		if (args.statusInterval > 0 && (move - stepStart) % args.statusInterval == 0)
		{
			std::cout << "Step " << move << ":\n--Current Energy: " << systemEnergy << std::endl;	
		}
		
		if (args.stateInterval > 0 && move > stepStart && (move - stepStart) % args.stateInterval == 0)
		{
			std::cout << std::endl;
			saveState(move);
			std::cout << std::endl;
		}
		
		int changeIdx = box->chooseMolecule();
		
		if (args.simulationMode == SimulationMode::Parallel)
		{
			oldEnergyCont = ParallelCalcs::calcMolecularEnergyContribution(box, changeIdx);
		}
		else
		{
			oldEnergyCont = SerialCalcs::calcMolecularEnergyContribution(molecules, enviro, changeIdx);
		}
			
		box->changeMolecule(changeIdx);
		
		if (args.simulationMode == SimulationMode::Parallel)
		{
			newEnergyCont = ParallelCalcs::calcMolecularEnergyContribution(box, changeIdx);
		}
		else
		{
			newEnergyCont = SerialCalcs::calcMolecularEnergyContribution(molecules, enviro, changeIdx);
		}
		
		bool accept = false;
		
		if(newEnergyCont < oldEnergyCont)
		{
			accept = true;
		}
		else
		{
			Real x = exp(-(newEnergyCont - oldEnergyCont) / kT);
			
			if(x >= randomReal(0.0, 1.0))
			{
				accept = true;
			}
			else
			{
				accept = false;
			}
		}
		
		if(accept)
		{
			accepted++;
			systemEnergy += newEnergyCont - oldEnergyCont;
		}
		else
		{
			rejected++;
			//restore previous configuration
			box->rollback(changeIdx);
		}
	}

	endTime = clock();
    double diffTime = difftime(endTime, startTime) / CLOCKS_PER_SEC;

	std::cout << "Step " << (stepStart + simSteps) << ":\r\n--Current Energy: " << systemEnergy << std::endl;
	systemEnergy = systemEnergy;

	// Save the final state of the simulation
	if (args.stateInterval >= 0)
	{
		saveState((stepStart + simSteps));
	}
	
	std::cout << std::endl << "Finished running " << simSteps << " steps" << std::endl;
	std::cout << "Final Energy: " << systemEnergy << std::endl;
	std::cout << "Run Time: " << diffTime << " seconds" << std::endl;
	std::cout << "Accepted Moves: " << accepted << std::endl;
	std::cout << "Rejected Moves: " << rejected << std::endl;
	std::cout << "Acceptance Ratio: " << 100.0 * accepted / (accepted + rejected) << '\%' << std::endl;
	
	saveResults(systemEnergy, diffTime, accepted, rejected);
}

void Simulation::saveState(int simStep)
{
	//determine where we want the state file to go
	std::string baseStateFile;	
	if (!args.simulationName.empty())
	{
		baseStateFile = args.simulationName;
	}
	else
	{
		baseStateFile = "untitled";
	}
	
	StateScanner statescan = StateScanner("");
	std::string stateOutputPath = baseFileName;
	std::string stepCount;

	if (!toString<int>(simStep, stepCount))
		return;

	stateOutputPath.append("_");
	stateOutputPath.append(stepCount); //add the step number to the name of the output file
	stateOutputPath.append(".state");

	std::cout << "Saving state file " << stateOutputPath << std::endl;

	statescan.outputState(box->getEnvironment(), box->getMolecules(), box->getMoleculeCount(), simStep, stateOutputPath);
}

void Simulation::saveResults()
{
	std::string resultsName;
	if (args.simulationName.empty())
		resultsName = RESULTS_FILE_DEFAULT;
	else
		resultsName = args.simulationName;
	resultsName.append(RESULTS_FILE_EXT);

	// Save the simulation results.
	std::ofstream resultsFile;
	resultsFile.open(resultsName.c_str());

	resultsFile << "######### MCGPU Results File #############" << std::endl;
	resultsFile << "[Information]" << std::endl;
	resultsFile << "Timestamp = " << currentDateTime() << std::endl;
	if (!args.simulationName.empty())
		resultsFile << "Simulation-Name = " << args.simulationName << std::endl;

	if (args.simulationMode == SimulationMode::Parallel)
		resultsFile << "Simulation-Mode = GPU" << std::endl;
	else
		resultsFile << "Simulation-Mode = CPU" << std::endl;
	resultsFile << "Starting-Step = " << stepStart << std::endl;
	resultsFile << "Steps = " << simSteps << std::endl;
	resultsFile << "Molecule-Count = " << box->environment->numOfMolecules << std::endl << std::endl;
	resultsFile << "[Results]" << std::endl;
	resultsFile << "Final-Energy = " << currentEnergy << std::endl;
	resultsFile << "Run-Time = " << diffTime << " seconds" << std::endl;
	resultsFile << "Accepted-Moves = " << accepted << std::endl;
	resultsFile << "Rejected-Moves = " << rejected << std::endl;
	resultsFile << "Acceptance-Rate = " << 100.0f * accepted / (float) (accepted + rejected) << '\%' << std::endl;

	resultsFile.close();
}

const std::string Simulation::currentDateTime()
{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);

    return buf;
}
