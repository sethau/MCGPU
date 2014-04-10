/*
	New version of SimBox
	Minimized to include only Atoms and Molecules

	Author: Nathan Coleman
	Last Changed: February 21, 2014
	-> February 26, by Albert Wallace
*/

#ifndef SERIALBOX_H
#define SERIALBOX_H

#include "Metropolis/Box.h"

class SerialBox : public Box
{
	public:
		SerialBox();
		~SerialBox();

		int molecTypenum;
		Table *tables;
};

#endif