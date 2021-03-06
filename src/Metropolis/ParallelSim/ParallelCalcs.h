/*
	Declares ParallellCalcs methods accessed from C++ compiled code (no Cuda).

	Created: February 21, 2014
*/

#ifndef PARALLELCALCS_H
#define PARALLELCALCS_H

#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <math.h>
#include <string>
#include "Metropolis/Box.h"
#include "Metropolis/DataTypes.h"
#include "Metropolis/SimulationArgs.h"

namespace ParallelCalcs
{
	const double kBoltz = 0.00198717;
	
	/// Factory method for creating a Box from a configuration file.
	/// @param configpath The path to the configuration file.
	/// @param steps The number of steps desired in the simulation,
	/// @return Returns a pointer to the filled-in Box.
	/// @note This functionality should ideally reside in ParallelBox,
	///   but it was placed here due to time constraints.
	///   TODO for future group.
	Box* createBox(std::string inputPath, InputFileType inputType, long* startStep, long* steps);
	
	/// Runs the simulation using step parallelism. Although this
	///   method closely mirrors Simulation::run(), it needs access
	///   to methods specific to ParallelBox, which cannot be included
	///   in the cpp-restricted Simulation class.
	void runParallelSteps(int simSteps, Box *box, Real &systemEnergy, int &accepted, int &rejected);
	
	/// Calculates the system energy using consecutive calls to
	///   calcMolecularEnergyContribution.
	/// @param box A ParallelBox cast as a Box, passed from Simulation.
	/// @return Returns total system energy.
	Real calcSystemEnergy(Box *box);
	
	/// Calculates the inter-molecular energy contribution of a given molecule,
	///   without intramolecular energy, using a batch method.
	/// @param box A ParallelBox cast as a Box, passed from Simulation.
	/// @param molIdx the index of the current changed molecule.
	/// @param startIdx The optional starting index for other molecules.
	///   Used for system energy calculation.
	/// @return Returns total molecular energy contribution, without
	///   intramolecular energy.
	Real calcMolecularEnergyContribution(Box *box, int molIdx, int startIdx = 0);
}

#endif