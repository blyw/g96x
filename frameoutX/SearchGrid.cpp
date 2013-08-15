#define EIGEN_NO_DEBUG  
#define EIGEN_MPL2_ONLY
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include "SearchGrid.h"

//build search grid
//this creates a table with 3 rows and n-column for progressive search through neighbouring boxes/images
//should consider adding some randomness in this part of the algorithm; this will speed the search process, but only if not using Eigen functions
//make the search space even larger than earlier defined
SearchGrid::SearchGrid(void)
{
}

SearchGrid::~SearchGrid(void)
{
}

void SearchGrid::BuildRectangular(std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grids, int number_of_neighbour_per_dimension)
{
    Eigen::Matrix<int, 3, Eigen::Dynamic> ref_x;

    int grid_counter = 0, x = 0, y = 0, z = 0; 
    ref_x.resize(3, (number_of_neighbour_per_dimension*2+1)*(number_of_neighbour_per_dimension*2+1)*(number_of_neighbour_per_dimension*2+1));
    for (x = 0-number_of_neighbour_per_dimension; x <= number_of_neighbour_per_dimension; x++)
    {
        for (y = 0-number_of_neighbour_per_dimension; y <= number_of_neighbour_per_dimension; y++)
        {
            for (z = 0-number_of_neighbour_per_dimension; z <= number_of_neighbour_per_dimension; z++)
            {
                ref_x(0,grid_counter) = x;
                ref_x(1,grid_counter) = y;
                ref_x(2,grid_counter) = z;
                grid_counter += 1;
            }
        }
    }
    grids->push_back(ref_x);
}

void SearchGrid::BuildRectangularExtended(std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grids, int number_of_neighbour_per_dimension, int extended_by_n)
{
    Eigen::Matrix<int, 3, Eigen::Dynamic> ref_x;
    int grid_counter = 0, x = 0, y = 0, z = 0; 
    ref_x.resize(3, ((number_of_neighbour_per_dimension+extended_by_n)*2+1)*((number_of_neighbour_per_dimension+extended_by_n)*2+1)*((number_of_neighbour_per_dimension+extended_by_n)*2+1));
    for (x = 0-(number_of_neighbour_per_dimension+extended_by_n); x <= (number_of_neighbour_per_dimension+extended_by_n); x++)
    {
        for (y = 0-(number_of_neighbour_per_dimension+extended_by_n); y <= (number_of_neighbour_per_dimension+extended_by_n); y++)
        {
            for (z = 0-(number_of_neighbour_per_dimension+extended_by_n); z <= (number_of_neighbour_per_dimension+extended_by_n); z++)
            {
                ref_x(0,grid_counter) = x;
                ref_x(1,grid_counter) = y;
                ref_x(2,grid_counter) = z;
                grid_counter += 1;
            }
        }
    }
    grids->push_back(ref_x);
}

//not implemented
void SearchGrid::BuildTruncatedOctahedral(std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grids, int number_of_neighbour_per_dimension) {
}

//not implemented
void SearchGrid::BuildTruncatedOctahedral(std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grids, int number_of_neighbour_per_dimension, int extended_by_n) {
}