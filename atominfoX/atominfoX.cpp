#define DEBUG
#define EIGEN_NO_DEBUG  
#define EIGEN_MPL2_ONLY
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include "InputParameters.h"
#include "Structs.h"

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

	auto end = std::chrono::system_clock::now();
	auto diff = end - start;
	std::cout << "        " << std::chrono::duration <double, std::milli>(diff).count() << " ms" << std::endl;
	return 0;
}