/*
	Driver for the simulation. Takes in a SimulationArgs object and creates the
	the necessary Box type, state file output path, etc.
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

void Simulation::runParallelSteps()
{
	std::cout << "Simulation Name: " << args.simulationName << std::endl;
	
	Real systemEnergy = 0;
	int accepted = 0, rejected = 0;

	clock_t startTime, endTime;
    startTime = clock();
	
	//calculate old energy
	if (systemEnergy == 0)
	{
		systemEnergy = ParallelCalcs::calcSystemEnergy(box);
	}
	
	std::cout << std::endl << "Running " << simSteps << " steps" << std::endl << std::endl;
	
	ParallelCalcs::runParallelSteps(simSteps, box, systemEnergy, accepted, rejected);
	
	/*int move = stepStart, stepsUntilNextBreak;
	while (move < stepStart + simSteps)
	{
		if (args.statusInterval > 0 && (move - stepStart) % args.statusInterval == 0)
		{
			std::cout << "Step " << move << ":\n--Current Energy: " << systemEnergy << std::endl;	
		}
		
		if (args.stateInterval > 0 && move > stepStart && (move - stepStart) % args.stateInterval == 0)
		{
			saveState(move);
		}
		
		stepsUntilNextBreak = PAR_STEPS;//If we want to do it this way, this will not be PAR_STEPS
		
		ParallelCalcs::runParallelSteps(stepsUntilNextBreak, &systemEnergy, &accepted, &rejected);
		move += PAR_STEPS;
	}*/

	endTime = clock();
    double diffTime = difftime(endTime, startTime) / CLOCKS_PER_SEC;

	std::cout << "Step " << (stepStart + simSteps) << ":\r\n--Current Energy: " << systemEnergy << std::endl;

	// Save the final state of the simulation
	if (args.stateInterval >= 0)
	{
		saveState(stepStart + simSteps);
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
	if (args.simulationMode == SimulationMode::Parallel)
	{
		runParallelSteps();
		return;
	}
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
		saveState(stepStart + simSteps);
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
	std::cout << std::endl;
	
	//determine where we want the state file to go
	std::string stateOutputPath;	
	if (!args.simulationName.empty())
	{
		stateOutputPath = args.simulationName;
	}
	else
	{
		stateOutputPath = "untitled";
	}
	
	StateScanner statescan = StateScanner("");
	std::string stepCount;

	if (!toString<int>(simStep, stepCount))
		return;

	stateOutputPath.append("_");
	stateOutputPath.append(stepCount); //add the step number to the name of the output file
	stateOutputPath.append(".state");

	std::cout << "Saving state file " << stateOutputPath << std::endl;

	statescan.outputState(box->getEnvironment(), box->getMolecules(), box->getMoleculeCount(), simStep, stateOutputPath);
	
	std::cout << std::endl;
}

void Simulation::saveResults(Real systemEnergy, double diffTime, int accepted, int rejected)
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
	resultsFile << "Final-Energy = " << systemEnergy << std::endl;
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
