#pragma once
#include <Eigen/Dense>
#include "Structs.h"

class Gather
{
public:
    Gather(void);
    ~Gather(void);
    static void FirstAtomBasedBoxShifter(Structs::FrameGeometric *framedata, int atom_number, Structs::InputParametersFrameout *me, Eigen::Matrix<int, 3, Eigen::Dynamic> *grid);
    static void SoluteMolecule(Structs::FrameGeometric *framedata, Structs::InputParametersFrameout *me, std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grid);
    static void SoluteCenterOfGeometry(Structs::FrameGeometric *framedata, Structs::InputParametersFrameout *me, std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grid);
    static void IonsCenterOfGeometry(Structs::FrameGeometric *framedata, Structs::InputParametersFrameout *me, std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grid);
    static void Solvent(Structs::FrameGeometric *framedata, Structs::InputParametersFrameout *me, std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grid);
    static void Solvent(Structs::FrameGeometric *framedata, Structs::InputParametersFrameout *me, std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grid, double solvent_sphere_radius);
};

