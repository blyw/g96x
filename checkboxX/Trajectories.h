#pragma once
#include "Structs.h"
#include <Eigen/Dense>
#include "FrameGeometry.h"
#include <vector>
#include <string>

class Trajectories
{
public:
    std::vector<Structs::FrameGeometric> activeFrames;
    std::vector<Structs::FrameNewtonian> framesNewtonian;
    std::string state;
    int activeFrame_counter;
    //holds the active frames which is equal to the number of frames handle per thread
    int activeFrame_count;    

    Trajectories(void);

    Trajectories(Structs::GenericParameters *me, 
        std::vector<std::string> *prefix, Eigen::MatrixXd *coordinates,
        bool *firstPass, int *frame_counter, double *frame_time, bool *processThisFrame);
    ~Trajectories(void);

    void ReadGeometric(gz::igzstream &file, gz::ogzstream &outfile);

private:
    Structs::GenericParameters *me;
    std::vector<std::string> *prefix;
    Eigen::MatrixXd *coordinates;
    bool *firstPass;
    int *frame_counter;
    double *frame_time;
    bool *processThisFrame;
    Structs::FrameGeometric currentFrame;
};