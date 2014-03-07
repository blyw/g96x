#define DEBUG
#define EIGEN_NO_DEBUG  
#define EIGEN_MPL2_ONLY
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include "FrameGeometry.h"
#include "InputParameters.h"
#include "SearchGrid.h"
#include "Structs.h"
#include "Trajectories.h"
#include "gzstream.h"
#include "DistanceMatrix.h"

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

	//some other parallelization parameters
	bool done_reading = false;
	bool done_calculating = false;
	bool done_writting = false;

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
	InputParameters ip(param_file);
	Structs::GenericParameters me;
	me.ref_coords.setZero();
	ip.ParseInputFile(&me, job_id);
	ip.PrintGenericParameters(me);

#ifdef DEBUG
	std::cout << "--->determine output file extension" << std::endl;
#endif // DEBUG

	//define output file
	gz::ogzstream outfile;
	outfile.open((me.outfilename + ".dismat.dat.gz").c_str(), std::ios_base::out);

#ifdef DEBUG
	std::cout << "--->print frameout input parameters to file" << std::endl;
#endif // DEBUG

	//print input parameters to output file
	ip.PrintGenericParameters(me, outfile);

	//print program data output header
	outfile << "#" << std::right << std::setw(9) << "min" << std::setw(9) << "max" << std::setw(10) << " time" << std::setw(10) << " timestep" << std::endl;

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
			writeAtomRecordsCount = me.solute_cog_molecules(1, solutes_cog_count - 1);
		}
		else
		{
			writeAtomRecordsCount = me.solute_molecules(1, me.solute_count - 1);
		}
	}
	else {
		writeAtomRecordsCount = me.atomrecords;
	}

	//build the search grids in advance and only reference to them in the rest of the code
	Structs::FrameoutReferenceGrids grids;
	SearchGrid::BuildRectangularForCheckbox(&grids.Checkbox);

#ifdef DEBUG
	std::cout << "checkbox search grid" << std::endl;
	std::cout << grids.Checkbox << std::endl;
#endif // DEBUG

#ifdef DEBUG
	std::cout << "--->initialize a bunch of variables" << std::endl;
#endif // DEBUG

	//remember if box has been shift to gather first atom of a frame
	//this is only used for gathering the first atom of each frame
	Eigen::Vector3d init_shift(0, 0, 0);

	int trajs_buffer_size = 6;
	std::vector<Trajectories> trajs(trajs_buffer_size);

#ifdef DEBUG
	std::cout << "--->initialize the output file thread" << std::endl;
#endif // DEBUG

	//output thread
	std::thread outfileThread = std::thread([&me, &outfile, &trajs, &trajs_buffer_size, &done_calculating, &done_reading, &done_writting](){

		std::chrono::milliseconds dura(2000);
		for (int i_trajs_buffer = 0; i_trajs_buffer < trajs_buffer_size; i_trajs_buffer++)
		{
#ifdef DEBUG
			for (int q = 0; q < trajs_buffer_size; q++)
			{
				std::cout << "### " << i_trajs_buffer << " " << q << " " << (trajs[q].state == "calculated") << " " << trajs[q].state << " " << trajs[q].activeFrame_counter << " " << std::endl;
			}
#endif // DEBUG

			while (!(trajs[i_trajs_buffer].state == "calculated"))
			{
				if (done_reading && done_calculating)
				{
					done_writting = true;
					return;
				}
				std::this_thread::sleep_for(dura);
			}

			for (int i_frames = 0; i_frames < trajs[i_trajs_buffer].activeFrame_counter; i_frames++)
			{
				//DO SOMETHING HERE
				DistanceMatrix::WriteResults(&trajs[i_trajs_buffer].activeFrames[i_frames], &me, outfile);
			}

			trajs[i_trajs_buffer].state = "written";

			if (i_trajs_buffer == (trajs_buffer_size - 1))
			{
				i_trajs_buffer = -1;
			}
		}
	});

	//calculation thread
	std::thread calculationThread = std::thread([&me, &outfile, &trajs, &trajs_buffer_size, &done_reading, &grids, &done_calculating](){

		std::thread *myThreads = new std::thread[me.num_thread_real];

		std::chrono::milliseconds dura(2000);
		for (int i_trajs_buffer = 0; i_trajs_buffer < trajs_buffer_size; i_trajs_buffer++)
		{
#ifdef DEBUG
			for (int q = 0; q < trajs_buffer_size; q++)
			{
				std::cout << "###### " << i_trajs_buffer << " " << q << " " << (trajs[q].state == "read") << " " << trajs[q].state << " " << trajs[q].activeFrame_counter << " " << std::endl;
			}
#endif // DEBUG

			while (!(trajs[i_trajs_buffer].state == "read"))
			{

				if (done_reading)
				{
					done_calculating = true;
					return;
				}
				std::this_thread::sleep_for(dura);
			}

			try
			{
				int perThread = me.num_frame_per_thread;
				if (trajs[i_trajs_buffer].activeFrame_counter < trajs[i_trajs_buffer].activeFrame_count)
				{
					perThread = trajs[i_trajs_buffer].activeFrame_counter / me.num_thread_real;
				}

				for (int i_threads = 0; i_threads < me.num_thread_real; i_threads++)
				{
					myThreads[i_threads] = std::thread([i_threads, &me, &grids, &i_trajs_buffer, &trajs, &perThread, &dura]()
					{
						for (int i_frames = (i_threads * perThread); i_frames < ((i_threads + 1) * perThread); i_frames++)
						{
							//DO SOMETHING HERE
							DistanceMatrix::Calculate(&trajs[i_trajs_buffer].activeFrames[i_frames], &me);
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

				for (int i_remainders = (perThread * me.num_thread_real); i_remainders < trajs[i_trajs_buffer].activeFrame_counter; i_remainders++)
				{
					//DO SOMETHING HERE
					DistanceMatrix::Calculate(&trajs[i_trajs_buffer].activeFrames[i_remainders], &me);
				}

#ifdef DEBUG
				std::cout << "--->all thread finished successfully" << std::endl;
#endif // DEBUG
			}
			//in case there is an error
			catch (const std::exception &e) {
				std::wcout << "\nEXCEPTION (calculation threads): " << e.what() << std::endl;
			}

			trajs[i_trajs_buffer].state = "calculated";

			if (i_trajs_buffer == (trajs_buffer_size - 1))
			{
				i_trajs_buffer = -1;
			}
		}
	});

	//input thread
	std::thread infileThread = std::thread([&me, &outfile, &trajs, &trajs_buffer_size, &done_reading](){

#ifdef DEBUG
		std::cout << "      created read thread" << std::endl;
#endif // DEBUG

		std::chrono::milliseconds dura(500);
		//defined if this is the first frame to be read
		//e.g. can be used to parse blocks that should only be processed once
		bool firstPass = true;
		//frame counters which can be used for skipping frames
		//keep track of the number of frames
		int frame_counter = 0;
		//keep track of the time elapsed
		double frame_time = 0;
		//use to keep track which frames to include for processing
		bool processThisFrame = false;
		//variables for processing the POSITION block
		//the first 24 columnns of the line
		std::vector<std::string> prefix(me.atomrecords);
		//3xn matrix for holding the actual coordinates
		Eigen::MatrixXd coordinates(0, 0);
		//set the number of rows and columns for the matrix
		coordinates.resize(3, me.atomrecords);

#ifdef DEBUG
		std::cout << "      read variables initialized" << std::endl;
#endif // DEBUG

		for (int i = 0; i < trajs_buffer_size; i++)
		{
			trajs[i] = Trajectories(&me, &prefix, &coordinates, &firstPass, &frame_counter, &frame_time, &processThisFrame);
		}

#ifdef DEBUG
		std::cout << "      reader data buffer objects created" << std::endl;
#endif // DEBUG

		//data files one-by-one as specified in the input parameters
		//then apply additional processing to the data

		//keeps track of which file is being read
		int file_counter = 0;
		//keeps track of which buffer is being filled
		int i_trajs_buffer = 0;
		//loop through all input data files
		while (file_counter < me.input_files.size())
		{
			//read a file containing data
			gz::igzstream file(me.input_files[file_counter].c_str());

			//loops through the buffers that hold the files
			//never reset the buffer counter as it is needed to tell the next data file which buffer to start filling
			for (i_trajs_buffer; i_trajs_buffer < trajs_buffer_size; i_trajs_buffer++)
			{
#ifdef DEBUG
				std::cout << "#####" << file_counter << " " << i_trajs_buffer << std::endl;
				std::cout << "      reading data file" << std::endl;
#endif // DEBUG

				//wait until the current buffer is ready for reading new data into it
				while (!(trajs[i_trajs_buffer].state == "new" || trajs[i_trajs_buffer].state == "written"))
				{
					std::this_thread::sleep_for(dura);
				}

				//read new data into current buffer 
				trajs[i_trajs_buffer].ReadGeometric(file, outfile);

				//keep looping over the buffers until break
				if (i_trajs_buffer == trajs_buffer_size - 1)
				{
					i_trajs_buffer = -1;
				}

				//if the end of a file is reached, update the to next file
				if (file.eof())
				{
					//specify the next file
					i_trajs_buffer += 1;
					file_counter += 1;
					//exit the for-loop without resetting the trajectory buffer count
					break;
				}
			}
			//close input file
			file.close();
		}

		done_reading = true;
	});

#ifdef DEBUG
	std::cout << "--->close output file" << std::endl;
#endif // DEBUG

	std::chrono::milliseconds dura(500);
	while (!done_writting)
	{
		std::this_thread::sleep_for(dura);
	}

	//make sure the output thread is idle
	if (infileThread.joinable())
	{
		infileThread.join();
	}


	if (calculationThread.joinable())
	{
		calculationThread.join();
	}
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
	std::cout << "        " << std::chrono::duration <double, std::milli>(diff).count() << " ms" << std::endl;
	return 0;
}