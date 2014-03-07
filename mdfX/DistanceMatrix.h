#pragma once
#include "Structs.h"
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include "gzstream.h"

class DistanceMatrix
{
public:
	DistanceMatrix(void);
	~DistanceMatrix(void);
	static void Calculate(Structs::FrameGeometric *framedata, Structs::GenericParameters *me);
	static void WriteResults(Structs::FrameGeometric *framedata, Structs::GenericParameters *me, gz::ogzstream &outfile);
};