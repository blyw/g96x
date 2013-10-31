#define DEBUG
#define EIGEN_NO_DEBUG  
#define EIGEN_MPL2_ONLY
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include "FrameGeometry.h"
#include "Gather.h"
#include "InputParameters.h"
#include "SearchGrid.h"
#include "Structs.h"

#include "gzstream.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <chrono>
#include <thread>
#include <future>
#include <cmath>
#include <Eigen/Dense>
#include <vector>

//main program
int main(int argc, char* argv[])
{
    Eigen::initParallel();

#ifdef DEBUG
    std::cout << "--->get CLI arguments" << std::endl;
#endif // DEBUG

    //CLI user input arguments
    std::string job_id = argv[1];
    std::string param_file = argv[2];

    //performance log - set start time
    auto start = std::chrono::system_clock::now();

#ifdef DEBUG
    std::cout << "--->parse frameout input parameters" << std::endl;
#endif // DEBUG

    //parse input parameters
    InputParameters ip (param_file);
    Structs::InputParametersFrameout me;
    me.ref_coords.setZero();
    ip.ParseInputFile(&me, job_id);
    ip.PrintInputParametersFrameOut(me);

#ifdef DEBUG
    std::cout << "--->determine output file extension" << std::endl;
#endif // DEBUG

    //define output file
    gz::ogzstream outfile;
    if (me.outformat == "pdb")
    {
        outfile.open((me.outfilename + "." + me.outformat + ".gz").c_str(), std::ios_base::out);
    }
    else if (me.outformat == "cnf")
    {
        outfile.open((me.outfilename + "." + me.outformat + ".gz").c_str(), std::ios_base::out);
    }

#ifdef DEBUG
    std::cout << "--->print frameout input parameters to file" << std::endl;
#endif // DEBUG

    //print input parameters to output file
    ip.PrintInputParametersFrameOut(me, outfile);

#ifdef DEBUG
    std::cout << "--->determine how many atom records are valid" << std::endl;
#endif // DEBUG

    //decide based on input, how many lines to write out in final output file
    int writeAtomRecordsCount = 0;
    if (me.solvent_skip)
    {
        int solutes_cog_count = me.solute_cog_molecules.cols();
        if (solutes_cog_count > 0)
        {                            
            writeAtomRecordsCount= me.solute_cog_molecules(1,solutes_cog_count-1);
        }
        else
        {
            writeAtomRecordsCount = me.solute_molecules(1,me.solute_count-1);
        }
    }
    else {
        writeAtomRecordsCount = me.atomrecords;
    }

    //build the search grids in advance and only reference to them in the rest of the code
    Structs::FrameoutReferenceGrids grids;

#ifdef DEBUG
    std::cout << "-->search space per molecule:" << std::endl;  
#endif // DEBUG

    //search grids for solute molecules
    for (int x = 0; x < me.solute_molecules.cols(); x++)
    {        
        SearchGrid::BuildRectangular(&grids.SolutesGatherer, me.solute_molecules(2,x));

#ifdef DEBUG
        std::cout << "solute molecule " << x << std::endl;
        std::cout << grids.SolutesGatherer[x] << "\n" << std::endl; 
#endif // DEBUG

    }

    //search grid for gathering the n-th atom or reference atom of each frame with respect to its previous frame
    //this search space is many folds larger than the largest user defined search space for solute molelecules
    // will use a default value of 2 additional images in each direction from the center box of the solute which specifies the
    //largest number of images
    if (me.solute_molecules.cols() <= 0)
    {
        SearchGrid::BuildRectangular(&grids.firstAtomBasedBoxShifter, 2);
    }
    else 
    {
        Eigen::MatrixXd::Index maxIndex;
        me.solute_molecules.row(2).maxCoeff(&maxIndex);

        SearchGrid::BuildRectangularExtended(&grids.firstAtomBasedBoxShifter, me.solute_molecules(2,maxIndex), 4);
    }

#ifdef DEBUG
    std::cout << "first atom of the frame" << std::endl;  
    std::cout << grids.firstAtomBasedBoxShifter[0] << std::endl;
#endif // DEBUG

    //search grids for solute molecules that are to be gather with respect to the COG of other solute molecules
    for (int x = 0; x < me.solute_cog_molecules.cols(); x++)
    {
        SearchGrid::BuildRectangular(&grids.SoluteCOGGatherer, me.solute_cog_molecules(3,x));

#ifdef DEBUG
        std::cout << "solute cog molecule " << x << std::endl;
        std::cout << grids.SoluteCOGGatherer[x] << "\n" << std::endl;  
#endif // DEBUG

    }
    //search grids for ion molecules that are to be gather with respect to the COG of all solute molecules
    for (int x = 0; x < me.ion_molecules.cols(); x++)
    {        
        SearchGrid::BuildRectangular(&grids.IonsGatherer, me.ion_molecules(3,x));

#ifdef DEBUG
        std::cout << "ion molecule " << x << std::endl;
        std::cout << grids.IonsGatherer[x] << "\n" << std::endl;  
#endif // DEBUG

    }
    //search grids for solvent molecules that are to be gather with respect to all other molecules
    for (int x = 0; x < me.solvent_molecules.cols(); x++)
    {
        SearchGrid::BuildRectangular(&grids.SolventGatherer, me.solvent_molecules(3,x));

#ifdef DEBUG
        std::cout << "solvent molecule " << x << std::endl;
        std::cout << grids.SolventGatherer[x] << "\n" << std::endl;  
#endif // DEBUG

    }

#ifdef DEBUG
    std::cout << "--->initialize a bunch of variables" << std::endl;
#endif // DEBUG

    //holds the active frames which is equal to the number of frames handle per thread
    int activeFrame_count = me.num_thread_real * me.num_frame_per_thread;

    std::vector<Structs::FrameGeometric> activeFrames (activeFrame_count);
    std::vector<Structs::FrameGeometric> activeFramesCopy (activeFrame_count);

    //defined if this is the first frame to be read
    //e.g. can be used to parse blocks that should only be processed once
    bool firstPass = true;
    Structs::FrameGeometric rotationalFitFrame;

    //frame counters which can be used for skipping frames
    //keep track of the number of frames
    int frame_counter = 0;
    //keep track of the time elapsed
    double frame_time = 0;
    //use to keep track which frames to include for processing
    bool processThisFrame = false;

    //variables for processing the POSITION block
    //the first 24 columnns of the line
    std::vector<std::string> prefix (me.atomrecords);
    //3xn matrix for holding the actual coordinates
    Eigen::MatrixXd coordinates(0,0);
    //set the number of rows and columns for the matrix
    coordinates.resize(3,me.atomrecords);

    //remember if box has been shift to gather first atom of a frame
    //this is only used for gathering the first atom of each frame
    Eigen::Vector3d init_shift(0,0,0);

    //a bunch of variables use for loops
    unsigned int i_input_files = 0;

#ifdef DEBUG
    std::cout << "--->if raw trajector file, parse a reference file for the first 24 columns" << std::endl;
#endif // DEBUG

    //parse reference file for prefix if it is a trajectory file in TRC format
    if (me.informat=="trc")
    {
        FrameGeometry::TrcReferenceFrame(&prefix, me.trc_reference);
    }

#ifdef DEBUG
    std::cout << "--->initialize the output file thread" << std::endl;
#endif // DEBUG

    //output thread
    std::thread outfileThread = std::thread([](){return 0;});

#ifdef DEBUG
    std::cout << "      processsing input file: " << me.input_files[i_input_files] << std::endl;
#endif // DEBUG

    //data files one-by-one as specified in the input parameters
    //then apply additional processing to the data
    for (i_input_files = 0; i_input_files < me.input_files.size(); i_input_files++)
    {

#ifdef DEBUG
        std::cout << "      initialize addition variable: " << std::endl;
#endif // DEBUG

        //define file to read
        gz::igzstream file(me.input_files[i_input_files].c_str());

        //boolean specifying current active block i.e. which block is being read
        bool isTitleBlock = false;
        bool isTimestepBlock = false;
        bool isPositionBlock = false;
        bool isGenboxBlock = false;

        //variables for holding the string content of the respective blocks
        std::string titleBlock, timestepBlock, positionBlock, genboxBlock;

        //additional counters
        //keep track of which line from the GENBOX block to process
        int positionBlock_counter = 0;
        //keep track how frames have been added
        int activeFrame_counter = 0;
        //keep track of which line from the GENBOX block to process
        int genBox_counter = 0;

#ifdef DEBUG
        std::cout << "      checking if data file is accessible" << std::endl;
#endif // DEBUG

        //check if it is possible to read file
        if (!file)
        {
            std::cerr << "          cannot open data file" << "\n";
        }
        else 
        {
            //define holder for current frame
            Structs::FrameGeometric currentFrame;

#ifdef DEBUG
            std::cout << "      reading data file" << std::endl;
#endif // DEBUG

            //read line-by-line as string while not end-of-file
            std::string line;
            while (!file.eof()) {
                std::getline(file, line);

                //ignore comments and empty lines
                if (line[0] != '#' && line.length() > 0)
                {
                    //detect block type found and set state for specific block
                    if (line.substr(0,6) == "TITLE")
                    {
                        isTitleBlock = true;
                    }
                    else if (line.substr(0,8) == "TIMESTEP")
                    {
                        isTimestepBlock = true;
                    }
                    else if (line.substr(0,8) == "POSITION")
                    {
                        positionBlock_counter = 0;
                        isPositionBlock = true;
                    }
                    else if (line.substr(0,6) == "GENBOX")
                    {
                        genBox_counter = 0;
                        isGenboxBlock = true;
                    }
                    //end of a block found, reset state to neutral
                    //  and do some post-processing if needed
                    else if (line.substr(0,3) == "END")
                    {
                        //what to do if end of a TITLE block
                        if (isTitleBlock)
                        {
                            //TITLE block
                            //only wirte something to the output file if the user desires CNF format
                            if (me.outformat == "cnf" && firstPass)
                            {
                                outfile << "TITLE" << "\n";
                                outfile << titleBlock;
                                outfile << "END" << "\n";
                            }
                        }
                        //what to do if end of a TIMESTEP block
                        if (isTimestepBlock)
                        {
#ifdef DEBUG
                            if (processThisFrame)
                            {
                                std::cout << "--->frame number " << frame_counter << std::endl;
                            }
#endif // DEBUG
                            //here we evaluate exclusion of the current frame
                            double time_interval;

                            //always process the first frame
                            if (frame_counter == 0)
                            {
                                processThisFrame = true;
                                //process the TIMESTEP block based on the input format
                                if (me.informat == "trc")
                                {
                                    frame_time = std::stod(timestepBlock.substr(15,15));
                                }
                                else if (me.informat == "cnf")
                                {
                                    frame_time = std::stod(timestepBlock.substr(19,19));                                    
                                }
                            }
                            //if skip-by-frames and nth frame and the current frame should be processed
                            else if (me.output_fragment_skipframes > 0 && (frame_counter % me.output_fragment_skipframes) == 0)
                            {
                                processThisFrame = true;
                            } 
                            //if skip-by-time-interval and interval and the current frame should be processed
                            else if (me.output_fragment_skiptime > 0)
                            {   
                                //calculate the time interval between the current and the last frame
                                //again, this is stupidly output program dependent
                                if (me.informat == "trc")
                                {
                                    time_interval = std::stod(timestepBlock.substr(15,15)) - frame_time;
                                }
                                else if (me.informat == "cnf")
                                {
                                    time_interval = std::stod(timestepBlock.substr(19,19)) - frame_time;
                                }

                                if (fmod(time_interval, me.output_fragment_skiptime) < 1e-16)
                                {
                                    processThisFrame = true;
                                }
                                else
                                {
                                    processThisFrame = false;
                                }
                            }
                            //if all output is desired
                            else if (me.output_fragment_skiptime <= 0 && me.output_fragment_skipframes <= 0)
                            {
                                processThisFrame = true;                                
                            }
                            else
                            {
                                processThisFrame = false;
                            }

                            //we can parse the TIMESTEP block, it should only contain one inline
                            //again, this is stupidly output program dependent
                            if (me.informat == "trc")
                            {
                                currentFrame.time = std::stod(timestepBlock.substr(15,15));
                                currentFrame.timestep = std::stol(timestepBlock.substr(0,15));
                            }
                            else if (me.informat == "cnf")
                            {
                                currentFrame.time = std::stod(timestepBlock.substr(18,20));
                                currentFrame.timestep = std::stol(timestepBlock.substr(0,18));
                            }

#ifdef DEBUG
                            if (processThisFrame)
                            {
                                std::cout << "      got the time and step ( " << currentFrame.time << " / " << frame_counter << " )" << std::endl;
                            }
#endif // DEBUG
                        }
                        //what to do if end of a POSITION block
                        if (isPositionBlock)
                        {
                            //reassign prefix and coordinates
                            currentFrame.prefix = prefix;
                            currentFrame.coordinates = coordinates;

#ifdef DEBUG
                            if (processThisFrame)
                            {
                                std::cout << "      got the coordinates " << std::endl;
                            }
#endif // DEBUG
                        } 
                        //what to do if end of a GENBOX block i.e. after a full frame has been read inWriteOutFrame
                        if (isGenboxBlock)
                        {
                            //the first frame has been collected, it's no longer the first pass through the loops
                            firstPass = false;

                            //first atom which is used to define reference between frames in the code below
                            //is defined as the first atom of the first solute molecule specified in the input file
                            //currently the user can specify any order of gathering
                            //MAKE IT POSSIBLE TO ALLOW THE USER TO SPECIFY AN ALTERNATIVE ATOM AS REFERENCE
                            int reference_solute_atom = 0;
                            if (me.solute_count >  0)
                            {
                                //THIS IS WRONG
                                reference_solute_atom = me.solute_molecules(3,0) - 1;
                                //reference_solute_atom = me.solute_molecules(0,0) - 1 + me.shift_reference_atom;
                            }

                            //skip frame[0] from the reference process, because it is the first frame 
                            //and has no reference coordinates for the first atom
                            //also get the time of the first frame as reference for time-based skipping
                            if (frame_counter==0)
                            {
                                currentFrame.init_shift.setZero();
                                me.ref_coords = coordinates.col(reference_solute_atom);
                            }
                            else 
                            {
                                //do not combine the following lines
                                //gather should be kept separate from reference coordinates assignment
                                Gather::FirstAtomBasedBoxShifter(&currentFrame, reference_solute_atom, &me, &grids.firstAtomBasedBoxShifter[0]);
                                //update reference coordinates using shifted coordinates.. do not use raw coordinates
                                me.ref_coords = currentFrame.coordinates.col(reference_solute_atom);
                            }
                            //set other frame properties
                            currentFrame.solute_cog_sum.setZero();
                            currentFrame.solute_cog_count = 0;
                            currentFrame.solute_cog.setZero();

                            //until n frames have been read in, store any frame that will be processed
                            if (processThisFrame)
                            {
                                activeFrames[activeFrame_counter] = currentFrame;
                                activeFrame_counter += 1;
                            }

                            //this is need for rotational fit, a reference frame is required
                            if (frame_counter==0 && me.correction_rotation && me.gather)
                            {
                                //we need to gather the first frame first and use it as a reference if we can to do rotational fit
                                Gather::SoluteMolecule(&currentFrame, &me, &grids.SolutesGatherer);
                                Gather::SoluteCenterOfGeometry(&currentFrame, &me, &grids.SoluteCOGGatherer);
                                Gather::IonsCenterOfGeometry(&currentFrame, &me, &grids.IonsGatherer);
                                Gather::Solvent(&currentFrame, &me, &grids.SolventGatherer, 0.5);
                                me.correction_translation = true;
                                currentFrame.coordinates.colwise() -= currentFrame.solute_cog;
                                rotationalFitFrame = currentFrame;
                            }

                            //if after processing the GENBOX block, n frames are stored than do something
                            if (activeFrame_counter == activeFrame_count)
                            {
#ifdef DEBUG
                                ////print out all time and corresponding timestep of n frames stored
                                //for (int x = 0; x < activeFrames.size(); x++)
                                //{
                                //    std::cout << activeFrames[x].time << " " << activeFrames[x].timestep << std::endl;  
                                //}                  
                                std::cout << "--->process the batch contains ( " << activeFrame_counter << " of " << activeFrame_count << " )" << std::endl;
#endif // DEBUG
                                //the calculation part starts here
                                //C++11 threading (domain decomposition)
                                std::thread *myThreads = new std::thread[me.num_thread_real];

#ifdef DEBUG
                                std::cout << "      calculation threads initialized" << std::endl;
#endif // DEBUG
                                if (me.gather)
                                {
                                    try
                                    {
                                        for (int i_threads = 0; i_threads < me.num_thread_real; i_threads++)
                                        {
                                            myThreads[i_threads] = std::thread([&activeFrames, &rotationalFitFrame, &me, i_threads, &grids]()
                                            {
                                                for (int i_frames = (i_threads * me.num_frame_per_thread); i_frames < ((i_threads + 1) * me.num_frame_per_thread); i_frames++)
                                                {
                                                    Gather::SoluteMolecule(&activeFrames[i_frames], &me, &grids.SolutesGatherer);
                                                    Gather::SoluteCenterOfGeometry(&activeFrames[i_frames], &me, &grids.SoluteCOGGatherer);
                                                    Gather::IonsCenterOfGeometry(&activeFrames[i_frames], &me, &grids.IonsGatherer);
                                                    if (me.solvent_sphere)
                                                    {                                                        
                                                        Gather::Solvent(&activeFrames[i_frames], &me, &grids.SolventGatherer, me.solvent_sphere_cut_off);
                                                    }
                                                    if (!me.solvent_sphere)
                                                    {                                                        
                                                        Gather::Solvent(&activeFrames[i_frames], &me, &grids.SolventGatherer);
                                                    }
                                                    if (me.correction_translation)
                                                    {
                                                        activeFrames[i_frames].coordinates.colwise() -= activeFrames[i_frames].solute_cog;
                                                        activeFrames[i_frames].solute_cog.setZero();
                                                        activeFrames[i_frames].solute_cog_sum.setZero();
                                                    }
                                                    if (me.correction_rotation)
                                                    {
                                                        FrameGeometry::CorrectionRotational(&activeFrames[i_frames], &rotationalFitFrame, &me);
                                                    }
                                                }
                                            });	
                                        }

                                        //wait for all threads to finish
                                        for (int i_threads = 0; i_threads < me.num_thread_real; i_threads++)
                                        {
                                            myThreads[i_threads].join();
#ifdef DEBUG
                                            std::cout << "--->joining thread: " << i_threads << std::endl;
#endif // DEBUG
                                        }

#ifdef DEBUG
                                        std::cout << "--->all thread finished successfully" << std::endl;
#endif // DEBUG
                                    }
                                    //in case there is an error
                                    catch (const std::exception &e) {
                                        std::wcout << "\nEXCEPTION (calculation threads): " << e.what() << std::endl;
                                    }

                                    //remove threads
                                    delete[] myThreads;
                                }

                                //check if the output writing thread has already finished
                                try
                                {
                                    if (outfileThread.joinable())
                                    {
                                        outfileThread.join();
                                    }
                                }
                                //in case there is an error
                                catch (const std::exception &e) {
                                    std::wcout << "\nEXCEPTION (join output thread): " << e.what() << std::endl;
                                }

                                //write out all processed frames sequentially
                                try                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
                                {
                                    //copy processed frames into a holder
                                    activeFramesCopy = activeFrames;
                                    outfileThread = std::thread([&activeFramesCopy, activeFrame_count, &me, &outfile, writeAtomRecordsCount]()
                                    {
                                        for (int x = 0; x < activeFrame_count; x++)
                                        {
                                            //DO SOMETHING
                                            FrameGeometry::WriteOutFrame(&activeFramesCopy[x], outfile, &me);
                                        }                                     
                                    });
                                }
                                catch (const std::exception &e) {
                                    std::cout << "\nEXCEPTION (restart output thread): " << e.what() << std::endl;
                                }
                                //reset activeframe_counter to 0
                                //i.e. get the next n frames for processing
                                activeFrame_counter = 0;
                            }

                            //clear all existing block content
                            titleBlock = "";
                            timestepBlock = "";
                            //increment frame_counter by 1, for keeping track which frame is being read in.
                            frame_counter += 1;
                        }

                        //reset booleans 
                        isTitleBlock = false;
                        isTimestepBlock = false;
                        isPositionBlock = false;
                        isGenboxBlock = false;
                    }
                    else 
                    {
                        //what to do if currently in a TITLE block
                        if (isTitleBlock)
                        {
                            titleBlock += line + "\n";
                        }
                        //what to do if currently in a TIMESTEP block
                        if (isTimestepBlock)
                        {
                            timestepBlock += line;
                        }  //what to do if currently in a POSITION block
                        if (isPositionBlock)
                        {
                            if (processThisFrame)
                            {
                                FrameGeometry::TrajectoryPositionBlockLineParser(me, line, positionBlock_counter, &prefix, &coordinates);
                            }
                            //a frame is not process the data of the previous frame is retained,
                            //but the data of the first atom is changed to match that of the current frame
                            else if (!processThisFrame && positionBlock_counter == 0)
                            {
                                FrameGeometry::TrajectoryPositionBlockLineParser(me, line, positionBlock_counter, &prefix, &coordinates);
                            }
                            else if (!processThisFrame) {
                            }
                            positionBlock_counter += 1;
                        }
                        //what to do if currently in a GENBOX block
                        if (isGenboxBlock)
                        {
                            FrameGeometry::TrajectoryGenboxBlockLineParser(&currentFrame, genBox_counter, line);
                            genBox_counter += 1;
                        }
                    }
                }
            }

            //what to do with remaining frames if end-of-file is reach
            if (activeFrame_counter > 0)
            {
                int perThread = activeFrame_counter / me.num_thread_real; 
#ifdef DEBUG
                ////print out all time and corresponding timestep of n frames stored
                //for (int x = 0; x < activeFrame_counter; x++)
                //{
                //    std::cout << activeFrames[x].time << " " << activeFrames[x].timestep << std::endl;  
                //}                  
                std::cout << "--->process the remaining frames ( " << activeFrame_counter << " of " << activeFrame_count << " )" << std::endl;                  
                std::cout << "--->  assigning " << perThread << " frames to each thread" << std::endl;
#endif // DEBUG

                //the calculation part starts here
                //C++11 threading (domain decomposition)
                std::thread *myThreads = new std::thread[me.num_thread_real];

#ifdef DEBUG
                std::cout << "      calculation threads initialized" << std::endl;
#endif // DEBUG


                if (me.gather)
                {
                    try
                    {
                        for (int i_threads = 0; i_threads < me.num_thread_real; i_threads++)
                        {
                            myThreads[i_threads] = std::thread([&activeFrames, &rotationalFitFrame, &me, i_threads, &perThread, &grids]()
                            {
                                for (int i_frames = (i_threads * perThread); i_frames < ((i_threads + 1) * perThread); i_frames++)
                                {
                                    Gather::SoluteMolecule(&activeFrames[i_frames], &me, &grids.SolutesGatherer);
                                    Gather::SoluteCenterOfGeometry(&activeFrames[i_frames], &me, &grids.SoluteCOGGatherer);
                                    Gather::IonsCenterOfGeometry(&activeFrames[i_frames], &me, &grids.IonsGatherer);
                                    if (me.solvent_sphere)
                                    {
                                        Gather::Solvent(&activeFrames[i_frames], &me, &grids.SolventGatherer, me.solvent_sphere_cut_off);
                                    }
                                    if (!me.solvent_sphere)
                                    {
                                        Gather::Solvent(&activeFrames[i_frames], &me, &grids.SolventGatherer);
                                    }
                                    if (me.correction_translation)
                                    {
                                        activeFrames[i_frames].coordinates.colwise() -= activeFrames[i_frames].solute_cog;
                                        activeFrames[i_frames].solute_cog.setZero();
                                        activeFrames[i_frames].solute_cog_sum.setZero();
                                    }
                                    if (me.correction_rotation)
                                    {
                                        FrameGeometry::CorrectionRotational(&activeFrames[i_frames], &rotationalFitFrame, &me);
                                    }
                                }
                            });	
                        }

                        //wait for all threads to finish
                        for (int i_threads = 0; i_threads < me.num_thread_real; i_threads++)
                        {
                            myThreads[i_threads].join();
#ifdef DEBUG
                            std::cout << "--->joining thread: " << i_threads << std::endl;
#endif // DEBUG
                        }

#ifdef DEBUG
                        std::cout << "--->all thread finished successfully" << std::endl;
#endif // DEBUG
                    }
                    //in case there is an error
                    catch (const std::exception &e) {
                        std::wcout << "\nEXCEPTION (calculation threads): " << e.what() << std::endl;
                    }

                    //remove threads
                    delete[] myThreads;

#ifdef DEBUG
                    std::cout << "--->spawn last threads to finish the job" << std::endl;
#endif // DEBUG
                    //spawn one thread per remaining frame
                    std::vector<std::thread> myThreadsRemainders;
                    for (int i_remainders = (perThread * me.num_thread_real); i_remainders < activeFrame_counter; i_remainders++)
                    {
                        myThreadsRemainders.push_back((std::thread([&activeFrames, &rotationalFitFrame, &me, i_remainders, &grids]()
                        {
                            Gather::SoluteMolecule(&activeFrames[i_remainders], &me, &grids.SolutesGatherer);
                            //20130824: some bug here
                            Gather::SoluteCenterOfGeometry(&activeFrames[i_remainders], &me, &grids.SoluteCOGGatherer);
                            Gather::IonsCenterOfGeometry(&activeFrames[i_remainders], &me, &grids.IonsGatherer);
                            if (me.solvent_sphere)
                            {
                                Gather::Solvent(&activeFrames[i_remainders], &me, &grids.SolventGatherer, me.solvent_sphere_cut_off);
                            }
                            if (!me.solvent_sphere)
                            {
                                Gather::Solvent(&activeFrames[i_remainders], &me, &grids.SolventGatherer);
                            }
                            if (me.correction_translation)
                            {
                                activeFrames[i_remainders].coordinates.colwise() -= activeFrames[i_remainders].solute_cog;
                                activeFrames[i_remainders].solute_cog.setZero();
                                activeFrames[i_remainders].solute_cog_sum.setZero();
                            }
                            if (me.correction_rotation)
                            {
                                FrameGeometry::CorrectionRotational(&activeFrames[i_remainders], &rotationalFitFrame, &me);
                            }
                        })));	
                    }

#ifdef DEBUG
                    std::cout << "--->wait for last "<< myThreadsRemainders.size() <<" threads to finish the job" << std::endl;
#endif // DEBUG
                    //wait for all threads to finish
                    for(int i_remainder = 0; i_remainder < myThreadsRemainders.size(); i_remainder++)
                    {
                        myThreadsRemainders[i_remainder].join();
                    }
                    myThreadsRemainders.clear();
                }
            }

            //make sure the output thread is idle
            if (outfileThread.joinable())
            {
                outfileThread.join();
            }

            //write out all processed frames sequentially
            try                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
            {
                //copy processed frames into a holder
                outfileThread = std::thread([&activeFrames, activeFrame_counter, &me, &outfile, writeAtomRecordsCount]()
                {
                    for (int x = 0; x < activeFrame_counter; x++)
                    {
                        //DO SOMETHING
                        FrameGeometry::WriteOutFrame(&activeFrames[x], outfile, &me);
                    }
                });
            }
            catch (const std::exception &e) {
                std::cout << "\nEXCEPTION (restart output thread): " << e.what() << std::endl;
            }

            //make sure the output thread is idle
            if (outfileThread.joinable())
            {
                outfileThread.join();
            }
        }
        file.close();
#ifdef DEBUG
        std::cout << "--->done with current input file, moving to next input file" << std::endl;  
#endif // DEBUG
    }

#ifdef DEBUG
    std::cout << "--->close output file" << std::endl;
#endif // DEBUG

    //make sure the output thread is idle
    if (outfileThread.joinable())
    {
        outfileThread.join();
    }
    //close output file
    outfile.close();

#ifdef DEBUG
    std::cout << "--->execution main program completed in:" << std::endl;
#endif // DEBUG

    //performance log - calculate execution time
    auto end = std::chrono::system_clock::now();
    auto diff = end - start;
    std::cout << "        " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
    return 0;
}