#pragma once
#define EIGEN_NO_DEBUG  
#define EIGEN_MPL2_ONLY
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include "Structs.h"
#include <string>
#include <vector>
#include <Eigen/Dense>

class Structs
{
public:
    Structs(void);
    ~Structs(void);

    struct FrameoutReferenceGrids {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> firstAtomBasedBoxShifter;
        std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> SolutesGatherer;
        std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> SoluteCOGGatherer;
        std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> IonsGatherer;
        std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> SolventGatherer;
        Eigen::Matrix<int, 3, 26> Checkbox;
    };
    
    //Holds settings for FRAMEOUT derived from parsing the user-provided input file 
    struct GenericParameters {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        //hardware parameters
        //number of threads to use
        int num_thread;
        //number of frames each thread should handle
        int num_frame_per_thread;
        //number of threads available supported by processor
        int num_thread_real;

        //output parameters
        //when to split output file into smaller files. done by counting frames
        //IMPLEMENT!!
        int output_fragment_size;
        //how many frames to skip before writing out the next frame.
        int output_fragment_skipframes;
        //how much time to skip before writing out the next frame
        double output_fragment_skiptime;
        //the format of the output file, either pdb or cnf
        std::string outformat;
        //the name of the output file, a absolute path is also allowed
        std::string outfilename;

        //topology parameters - solutes
        //the number of solute molecules in the simulation
        int solute_count;
        //the number of solutes (solutes_molecules) gather using nearest-image clustering method 
        //index: first atom, last atom, number of images to search, init_atom, skip or not
        Eigen::Matrix<int, 5, Eigen::Dynamic> solute_molecules;
        //the number of solutes (solute_cog_molecules) gather using nearest-image clustering method based on the center of geometry of a given set of solute molecules (solute_molecules)
        //index: first atom, last atom, number of atoms per solute_cog molecule, number of images to search, init_atom, skip or not
        Eigen::Matrix<int, 6, Eigen::Dynamic> solute_cog_molecules;
        //topology parameters - ions
        //ion molecules (ion_molecules) to be gather using nearest-image clustering method based on the center of geometry of all solute molecules (solute_molecules +  solute_cog_molecules) 
        //index: first atom, last atom, number of atoms per ion molecule, number of images to search, init_atom, skip or not
        Eigen::Matrix<int, 6, Eigen::Dynamic> ion_molecules;
        //to gather ions
        bool ions_skip;
        //topology parameters - solvent
        //index: first atom, last atom, number of atoms per solvent molecule, number of images to search
        Eigen::Matrix<int, 4, 1>  solvent_molecules;
        //to gather solvent
        bool solvent_skip;

        //input paramters
        //the input trajector format. either cnf or trc.
        std::string informat;
        //a reference frame in cnf format, from which the program can fall back to 
        //to get the atom description that correspond with the coordinates
        std::string trc_reference;
        //a list of trajectory files
        std::vector<std::string> input_files;

        //shared parameters
        //the number of atoms records in the frame. uses the same definition as the PDB format.
        int atomrecords;
        //the cut-off distance used for gathering. this improves performance by preventing 
        //searching through all previously defined atoms using the clustering gathering method 
        double distance_cut_off;
        //the verbosity of the out logs if any available or if multiple level available.
        //e.g. determines what checkbox will write out
        int verbosity;

        //frameoutX parameters
        //reference coordinates for gathering the first solute molecules atom based
        //on the coordinates of the same atom in the previous frame if any exist.
        Eigen::Vector3d ref_coords;
        //write out the center of geometry of the solute molecules
        bool cog_write;
        //shift the center of geometry to (0,0,0)
        bool correction_translation;
        //shift the center of geometry to (0,0,0)
        bool correction_rotation;
        //disable the gathering routing completely
        bool gather;

        //solvent sphere
        bool solvent_sphere;
        double solvent_sphere_cut_off;

    };

    //Holds space-time-related properties of a given simulation frame 
    struct FrameGeometric {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        //variables for the timestep block
        //the step number of the current tracjectory
        long timestep;
        //the time elapsed since the first timestep
        double time;

        //variables with atom specific data
        //this contains the first 24 characters of a line with coordinates in the POSITION block
        std::vector<std::string> prefix;
        //add some variable to hold data contain in within these 24 characters

        //variable with cartesian coordinates of atoms
        //coordinates are stored in a column-major matrix provided by the Eigen library
        Eigen::Matrix<double, 3, Eigen::Dynamic> coordinates;

        //the center of geometry calculated from all solute molecules (as defined by user input)
        Eigen::Vector3d solute_cog;
        //raw sum of all coordinates included so far
        Eigen::Vector3d solute_cog_sum;
        //raw number of set of coordinates included so far
        int solute_cog_count;

        //simulation box parameters
        //the shape of the box
        int boxtype;
        //the edge length in each dimension
        Eigen::Vector3d box_length;
        //the angles of the box corners
        Eigen::Vector3d box_angle;
        Eigen::Vector3d box_3;
        Eigen::Vector3d box_4;

        //how much the current box should be shifted to be in the similar space as the previous box
        Eigen::Vector3i init_shift;

        //use to keep tracck of all periodic interactions between atoms in a frame
        std::vector<std::vector<double>> periodic_interactions;
    };
    
    //Holds Newtonian properties of a given simulation frame 
    struct FrameNewtonian {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Eigen::Matrix<double, 3, Eigen::Dynamic> velocities;
        Eigen::Matrix<double, 3, Eigen::Dynamic> forces;
    };
};

