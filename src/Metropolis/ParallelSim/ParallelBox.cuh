/*
	Implements methods related to managing data between the host and device.
	Subclass of Box.

	Created: February 21, 2014
	
	-> February 26, by Albert Wallace
	-> March 28, by Joshua Mosby
	-> April 21, by Nathan Coleman
*/

#ifndef PARALLELBOX_H
#define PARALLELBOX_H

#include "Metropolis/Box.h"
#include "Metropolis/DataTypes.h"
#include "Metropolis/Utilities/Coalesced_Structs.h"

class ParallelBox : public Box
{
	private:
		Real *xD, *yD, *zD, *sigmaD, *epsilonD, *chargeD;
		int *atomsIdxD, *numOfAtomsD;
		
		/// Copies a specified molecule and all of its atoms
		///   over to the device. Called after changing a
		///   molecule in the simulation.
		/// @param changeIdx The index of the changed molecule.
		void writeChangeToDevice(int changeIdx);

	public:
		int *changedIndices, numChanged;
		Molecule *changedMols;
		AtomData *atomsH, *atomsD;
		Environment *environmentD;
		MoleculeData *moleculesH, *moleculesD;
		int *nbrMolsH, *nbrMolsD, *molBatchH, *molBatchD;
		Real *energiesD;
		int energyCount, maxMolSize;
		
		ParallelBox();
		~ParallelBox();
		
		/// Chooses N random molecules to be changed for a given
		///   batch of parallel simulation steps.
		/// @param N The number of molecules to be chosen.
		void chooseMolecules(int N);
		
		/// Changes all changed molecules in a random way.
		void changeMolecules();
		
		/// Given two indices in the changedIndices array,
		///   check the distance against the cutoff.
		/// @param mol1 The first index.
		/// @param mol2 The second index.
		/// @return Returns true if the molecules are within
		///   the cutoff distance of each other.
		bool changedMolsWithinCutoff(int mol1, int mol2);
		
		/// Given one index in the changedIndices array,
		///   toggle the status of the molecule between
		///   changed and rolled-back, using the backups.
		/// @param molIdx The index to be toggled.
		void toggleChange(int molIdx);
		
		/// Changes a specified molecule in a random way.
		///   This method overrides the virtual method of the
		///   same name in the parent class, Box, to add a
		///   call to writeChangeToDevice() after the change.
		/// @param molIdx The index of the molecule to be
		///   changed.
		/// @return Returns the index of the changed molecule.
		virtual int changeMolecule(int molIdx);
		
		/// Rolls back the previous changes to the specified
		///   molecule. This method overrides the virtual method
		///   of the same name in the parent class, Box, to add a
		///   call to writeChangeToDevice() after the rollback.
		/// @param molIdx The index of the molecule that was
		///   changed.
		/// @return Returns the index of the changed molecule.
		/// @note This method should ideally not take in any
		///   parameters, as only one molecule change is pending
		///   at any time. Box should store the index of the
		///   changed molecule.
		virtual int rollback(int molIdx);
		
		/// Allocates and copies over all Box and ParallelBox data
		///   to the device as needed for the simulation.
		/// @note This functionality logically belongs in the
		///   constructor for ParallelBox, but, at the time of
		///   instantiation, the Box does not yet have all of the
		///   necessary data to copy over, and so it is a separate
		///   method for now.
		void copyDataToDevice();
};

#endif