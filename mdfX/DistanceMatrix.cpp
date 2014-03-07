#define DEBUG
#define EIGEN_NO_DEBUG  
#define EIGEN_MPL2_ONLY
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include "DistanceMatrix.h"
#include "Structs.h"
#include "gzstream.h"
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>

DistanceMatrix::DistanceMatrix()
{
}


DistanceMatrix::~DistanceMatrix()
{
}

void DistanceMatrix::Calculate(Structs::FrameGeometric *framedata, Structs::GenericParameters *me)
{
	//calculate distances for each set of atoms
	for (int i_columns = 0; i_columns < me->distanceMatrixCutOff.size(); i_columns++)
	{
		std::vector<double> temp;

		//std::cout << "\t\t#######" << i_columns << std::endl;
		//always assign the smallest sets of particles to set A
		for (int i_atoms = 0; i_atoms < me->distanceMatrixSets(1, i_columns) - me->distanceMatrixSets(0, i_columns) + 1; i_atoms++)
		{
			Eigen::Matrix<double, 1, Eigen::Dynamic> tempMat1;
			//reassigning values to temp with each iteration
			//output is an array of half an symmetric matrix
			//matrix.block(i,j,p,q) --> block of size (p,q), starting at (i,j)
			//i = start row
			//j = start column
			//p = 3 --> x,y,z
			//j = last atom - first atom
			tempMat1 = (
				//setB
				framedata->coordinates.block(0, me->distanceMatrixSets(2, i_columns) - 1, 3, me->distanceMatrixSets(3, i_columns) - me->distanceMatrixSets(2, i_columns) + 1).colwise() -
				//setA reference
				framedata->coordinates.col(me->distanceMatrixSets(0, i_columns) - 1 + i_atoms)
				).cwiseAbs2().colwise().sum().cwiseSqrt();
			//if (me->distanceMatrixWriteRaw[i_columns])
			//{
			//	std::cout << "diff " << me->distanceMatrixSets(1, i_columns) + i_atoms << " " << me->distanceMatrixSets(3, i_columns) - me->distanceMatrixSets(2, i_columns) + 1 - i_atoms << " : " << tempMat1 << std::endl;
			//}
			for (int i_temp = 0; i_temp < tempMat1.cols(); i_temp++)
			{
				temp.push_back(tempMat1.data()[i_temp]);
			}
		}


		//LOOK FOR A MORE EFFICIENT WAY TO CALCULATE THE DIFFERENCE BETWEEN THE TWO MATRICES
		//Eigen::Matrix<double, 1, Eigen::Dynamic> tempMat1;
		//tempMat1 = (
		//	//setB
		//	framedata->coordinates.block(0, me->distanceMatrixSets(2, i_columns) - 1, 3, me->distanceMatrixSets(3, i_columns) - me->distanceMatrixSets(2, i_columns) + 1) -
		//	//setA reference
		//	framedata->coordinates.block(0, me->distanceMatrixSets(0, i_columns) - 1, 3, me->distanceMatrixSets(1, i_columns) - me->distanceMatrixSets(0, i_columns) + 1)
		//	).cwiseAbs2().colwise().sum().cwiseSqrt();
		//if (me->distanceMatrixWriteRaw[i_columns])
		//{
		//	std::cout << "diff : \n" << tempMat1 << std::endl;
		//}
		//for (int i_temp = 0; i_temp < tempMat1.cols(); i_temp++)
		//{
		//	temp.push_back(tempMat1.data()[i_temp]);
		//}

		//store results for this set sof atoms
		framedata->interatomicDistances.push_back(Eigen::Map<Eigen::RowVectorXd>(temp.data(), temp.size()));
	}
}

void DistanceMatrix::WriteResults(Structs::FrameGeometric *framedata, Structs::GenericParameters *me, gz::ogzstream &outfile)
{
	//use a for finding minimum distances in matrix/vector
	Eigen::MatrixXd::Index maxIndex;
	Eigen::MatrixXd::Index minIndex;
	std::vector<std::string> temp;

	//i don't remember why we don't use distanceMatrixSets, maybe the size is preset somewhere?
	for (int i_columns = 0; i_columns < me->distanceMatrixCutOff.size(); i_columns++)
	{
		for (int i_ref = me->distanceMatrixSets(0, i_columns) - 1; i_ref < me->distanceMatrixSets(1, i_columns); i_ref++)
		{
			for (int i_target = me->distanceMatrixSets(2, i_columns) - 1; i_target < me->distanceMatrixSets(3, i_columns); i_target++)
			{
				temp.push_back(framedata->prefix[i_ref] + " " + framedata->prefix[i_target]);
			}
		}
	}

	for (int i_columns = 0; i_columns < me->distanceMatrixCutOff.size(); i_columns++)
	{
		//gets the shortest distance
		framedata->interatomicDistances[i_columns].maxCoeff(&maxIndex);
		double longest_distance = framedata->interatomicDistances[i_columns](0, maxIndex);
		framedata->interatomicDistances[i_columns].minCoeff(&minIndex);
		double shortest_distance = framedata->interatomicDistances[i_columns](0, minIndex);

		if (longest_distance > me->distanceMatrixCutOff[i_columns])
		{
			outfile << std::setprecision(4) << std::right;
			outfile << std::setw(10) << temp[minIndex] << " ";
			outfile << std::setw(10) << shortest_distance << " ";
			outfile << std::setw(10) << temp[maxIndex] << " ";
			outfile << std::setw(10) << longest_distance << " ";
			outfile << std::setw(10) << framedata->time << " ";
			outfile << std::setw(10) << std::setprecision(0) << framedata->timestep;

			if (me->distanceMatrixWriteRaw[i_columns])
			{
				for (int i_distance = 0; i_distance < framedata->interatomicDistances[i_columns].cols(); i_distance++)
				{
					outfile << " " << std::setprecision(4) << std::setw(10) << std::fixed << framedata->interatomicDistances[i_columns].col(i_distance) << std::flush;
				}
				outfile << std::endl;
			}
		}
	}
}