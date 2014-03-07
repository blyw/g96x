#define DEBUG
#define EIGEN_NO_DEBUG  
#define EIGEN_MPL2_ONLY
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include "Trajectories.h"

Trajectories::Trajectories(void) {
}

Trajectories::Trajectories(Structs::GenericParameters *me,
	std::vector<std::string> *prefix, Eigen::MatrixXd *coordinates,
	bool *firstPass, int *frame_counter, double *frame_time, bool *processThisFrame)
{
	Trajectories::me = me;
	Trajectories::prefix = prefix;
	Trajectories::coordinates = coordinates;
	Trajectories::firstPass = firstPass;
	Trajectories::frame_counter = frame_counter;
	Trajectories::frame_time = frame_time;
	Trajectories::processThisFrame = processThisFrame;
	Trajectories::state = "new";
	//holds the active frames which is equal to the number of frames handle per thread
	Trajectories::activeFrame_count = me->num_thread_real * me->num_frame_per_thread;
	Trajectories::activeFrames = std::vector<Structs::FrameGeometric>(Trajectories::activeFrame_count);
}


Trajectories::~Trajectories(void)
{
}

void Trajectories::ReadGeometric(gz::igzstream &file, gz::ogzstream &outfile)
{
#ifdef DEBUG
	std::cout << "      initialize addition variable: " << std::endl;
#endif // DEBUG

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
	activeFrame_counter = 0;
	//keep track of which line from the GENBOX block to process
	int genBox_counter = 0;

#ifdef DEBUG
	std::cout << "      checking if data file is accessible" << std::endl;
#endif // DEBUG

	//read line-by-line as string while not end-of-file
	std::string line;
	while (!file.eof()) {
		std::getline(file, line);

		//ignore comments and empty lines
		if (line[0] != '#' && line.length() > 0)
		{
			//detect block type found and set state for specific block
			if (line.substr(0, 6) == "TITLE")
			{
				isTitleBlock = true;
			}
			else if (line.substr(0, 8) == "TIMESTEP")
			{
				isTimestepBlock = true;
			}
			else if (line.substr(0, 8) == "POSITION")
			{
				positionBlock_counter = 0;
				isPositionBlock = true;
			}
			else if (line.substr(0, 6) == "GENBOX")
			{
				genBox_counter = 0;
				isGenboxBlock = true;
			}
			//end of a block found, reset state to neutral
			//  and do some post-processing if needed
			else if (line.substr(0, 3) == "END")
			{
				//what to do if end of a TITLE block
				if (isTitleBlock)
				{
				}
				//what to do if end of a TIMESTEP block
				if (isTimestepBlock)
				{
#ifdef DEBUG
					if (*processThisFrame)
					{
						std::cout << "--->frame number " << *frame_counter << std::endl;
					}
#endif // DEBUG
					//here we evaluate exclusion of the current frame
					double time_interval;

					//always process the first frame
					if (*frame_counter == 0)
					{
						*processThisFrame = true;
						//process the TIMESTEP block based on the input format
						*frame_time = std::stod(timestepBlock.substr(19, 19));
					}
					//if skip-by-frames and nth frame and the current frame should be processed
					else if (me->output_fragment_skipframes > 0 && (*frame_counter % me->output_fragment_skipframes) == 0)
					{
						*processThisFrame = true;
					}
					//if skip-by-time-interval and interval and the current frame should be processed
					else if (me->output_fragment_skiptime > 0)
					{
						//calculate the time interval between the current and the last frame
						//again, this is stupidly output program dependent
						time_interval = std::stod(timestepBlock.substr(19, 19)) - *frame_time;

						if (fmod(time_interval, me->output_fragment_skiptime) < 1e-16)
						{
							*processThisFrame = true;
						}
						else
						{
							*processThisFrame = false;
						}
					}
					//if all output is desired
					else if (me->output_fragment_skiptime <= 0 && me->output_fragment_skipframes <= 0)
					{
						*processThisFrame = true;
					}
					else
					{
						*processThisFrame = false;
					}

					//we can parse the TIMESTEP block, it should only contain one inline
					//again, this is stupidly output program dependent
					currentFrame.time = std::stod(timestepBlock.substr(18, 20));
					currentFrame.timestep = std::stol(timestepBlock.substr(0, 18));
#ifdef DEBUG
					if (processThisFrame)
					{
						std::cout << "      got the time and step ( " << currentFrame.time << " / " << *frame_counter << " )" << std::endl;
					}
#endif // DEBUG
				}
				//what to do if end of a POSITION block
				if (isPositionBlock)
				{
					//reassign prefix and coordinates
					currentFrame.prefix = *prefix;
					currentFrame.coordinates = *coordinates;

#ifdef DEBUG
					if (*processThisFrame)
					{
						std::cout << "      got the coordinates " << std::endl;
					}
#endif // DEBUG
				}
				//what to do if end of a GENBOX block i.e. after a full frame has been read inWriteOutFrame
				if (isGenboxBlock)
				{
					//the first frame has been collected, it's no longer the first pass through the loops
					//*firstPass = false;

					//until n frames have been read in, store any frame that will be processed
					if (processThisFrame)
					{
						activeFrames[activeFrame_counter] = currentFrame;
						activeFrame_counter += 1;
					}

					//clear all existing block content
					titleBlock = "";
					timestepBlock = "";
					//increment frame_counter by 1, for keeping track which frame is being read in.
					*frame_counter += 1;

					//if after processing the GENBOX block, n frames are stored than do something
					if (activeFrame_counter == activeFrame_count)
					{
#ifdef DEBUG
						////print out all time and corresponding timestep of n frames stored        
						std::cout << "--->process the batch contains ( " << activeFrame_counter << " of " << activeFrame_count << " )" << std::endl;
#endif // DEBUG

						break;
					}
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
					if (*processThisFrame)
					{
						FrameGeometry::TrajectoryPositionBlockLineParser(*me, line, positionBlock_counter, prefix, coordinates);
					}
					//a frame is not process the data of the previous frame is retained,
					//but the data of the first atom is changed to match that of the current frame
					else if (!*processThisFrame && positionBlock_counter == 0)
					{
						FrameGeometry::TrajectoryPositionBlockLineParser(*me, line, positionBlock_counter, prefix, coordinates);
					}
					else if (!*processThisFrame) {
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
	state = "read";
}
