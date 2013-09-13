#pragma once
#include "Structs.h"
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include "gzstream.h"

class PeriodicityCheck
{
public:
    PeriodicityCheck(void);
    ~PeriodicityCheck(void);
    static void NearestImagesFinder(Structs::FrameGeometric *framedata, Structs::GenericParameters *me, Structs::FrameoutReferenceGrids *grids);
    static void WriteResults(Structs::FrameGeometric *framedata, Structs::GenericParameters *me, gz::ogzstream &outfile);
    //static void NearestImagesFinder(Structs::FrameGeometric *framedata, Structs::GenericParameters *me, std::vector<std::vector<int>> *grid_list, int startIndex, int endIndex);
};

