#pragma once
#include "Structs.h"
#include <Eigen/Dense>
#include "Structs.h"
#include <exception>

//build search grid
//this creates a table with 3 rows and n-column for progressive search through neighbouring boxes/images
//should consider adding some randomness in this part of the algorithm; this will speed the search process, but only if not using Eigen functions
class SearchGrid
{
public:
    SearchGrid(void);
    ~SearchGrid(void);
    static void BuildRectangular(std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grids, int number_of_neighbour_per_dimension);
    static void BuildRectangularExtended(std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grids, int number_of_neighbour_per_dimension, int extended_by_n);
    static void BuildTruncatedOctahedral(std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grids, int number_of_neighbour_per_dimension);
    static void BuildTruncatedOctahedral(std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grids, int number_of_neighbour_per_dimension, int extended_by_n);
    static void BuildRectangularForCheckbox(Eigen::Matrix<int, 3, 26> *grid);
};

