//#define DEBUG
#define EIGEN_NO_DEBUG  
#define EIGEN_MPL2_ONLY
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include "Gather.h"
#include "Structs.h"
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <math.h>

Gather::Gather(void)
{
}


Gather::~Gather(void)
{
}

//simple gathering of a specified atom in frame with respect to the given coordinates of an atom
//input: framedata, atom for correction, input parameters from user, a precalculated search grid
void Gather::FirstAtomBasedBoxShifter(Structs::FrameGeometric *framedata, int atom_number, Structs::InputParametersFrameout *me, Eigen::Matrix<int, 3, Eigen::Dynamic> *grid) {
    //assign a very large initial value
    double distance_shortest = 1E20;

    //reference coordinates
    Eigen::Vector3d *ref_coords = &me->ref_coords;
    //how much to shift the coordinates 
    Eigen::Vector3i min_shift(0,0,0);

    //shift the coordinates of the specified atom in the original framedata based on the previous frame-shift (if required)
    //this is done only for frames after the first frame as the initial reference values for the first frame should be (0,0,0)
    //x.all() only applies to arrays and matrices... not for vectors
    //std::cout << "        shift from previous frame \n" << framedata->init_shift.array() << "\n" << (framedata->init_shift.array() == 0).all() << std::endl;
    if (!(framedata->init_shift.array() == 0).all())
    {
        //    std::cout << framedata->coordinates << std::endl;
        framedata->coordinates.colwise() += framedata->box_length.cwiseProduct(framedata->init_shift.cast<double>());
        //    std::cout << framedata->coordinates << std::endl;
    }

    //gather the n-th atom of the frame
    Eigen::MatrixXd::Index minIndex;
    (((grid->cast<double>().array().colwise() * framedata->box_length.array()).matrix().colwise() + framedata->coordinates.col(atom_number)).colwise() - me->ref_coords).cwiseAbs2().colwise().sum().cwiseSqrt().minCoeff(&minIndex);

    ////debug and explaination code below, uncomment to understand above statement
    //std::cout << framedata->coordinates.col(atom_number) << std::endl << std::endl;
    //std::cout << me->ref_coords << std::endl << std::endl;
    //std::cout << grid.cast<double>().array() << std::endl << std::endl;
    //std::cout << framedata->box_length.array() << std::endl << std::endl;
    //std::cout << (grid.cast<double>().array().colwise() * framedata->box_length.array()).matrix() << std::endl << std::endl;
    //std::cout << (grid.cast<double>().array().colwise() * framedata->box_length.array()).matrix().colwise() + framedata->coordinates.col(atom_number) << std::endl << std::endl;
    //std::cout << temp1 << std::endl << std::endl;
    //std::cout << temp1.cwiseAbs2() << std::endl << std::endl;
    //std::cout << temp1.cwiseAbs2().colwise().sum() << std::endl << std::endl;
    //std::cout << temp1.cwiseAbs2().colwise().sum().cwiseSqrt() << std::endl << std::endl;
    //std::cout << temp1.cwiseAbs2().colwise().sum().cwiseSqrt().minCoeff(&minIndex) << std::endl << std::endl;
    //std::cout << &minIndex << std::endl << std::endl;
    //std::cout << grid.col(minIndex) << std::endl << std::endl;
    //return atom shift corresponding to shortest distance

    min_shift = grid->col(minIndex).cast<int>();

    //std::cout << "        shift index \n" << minIndex << std::endl;
    //std::cout << (((grid->cast<double>().array().colwise() * framedata->box_length.array()).matrix().colwise() + framedata->coordinates.col(atom_number)).colwise() - me->ref_coords).cwiseAbs2().colwise().sum().cwiseSqrt()(&minIndex) << std::endl;
    //std::cout << (((grid->cast<double>().array().colwise() * framedata->box_length.array()).matrix().colwise() + framedata->coordinates.col(atom_number)).colwise() - me->ref_coords).cwiseAbs2().colwise().sum().cwiseSqrt() << std::endl;
    //std::cout << "        shift from current frame \n" << min_shift << "\n" << (min_shift.array() == 0).all() << std::endl;

    //shift again if a frame-shift was required for the positioning the n-th atom properly
    if (!(min_shift.array() == 0).all())
    {
        //    std::cout << framedata->coordinates << std::endl;
        framedata->coordinates.colwise() += framedata->box_length.cwiseProduct(min_shift.cast<double>());
        //    std::cout << framedata->coordinates << std::endl;
    }

    //remember how much shifting is applied to this frame, this will be applied to the next frame, and therefore reduced the chance
    //  that gathering will go wrong
    framedata->init_shift += min_shift;

    //std::cout << "        shift for next frame \n" << framedata->init_shift.array() << "\n" << (framedata->init_shift.array() == 0).all() << std::endl;
}

void Gather::SoluteMolecule(Structs::FrameGeometric *framedata, Structs::InputParametersFrameout *me, std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grid)
{
    //these two varaible are required to ensure that gathering will always occured
    double shortestDistance = 1e20;
    double currentDistance = 0;

    Eigen::Vector3i min_shift(0,0,0);

    //use a for finding minimum distances in matrix/vector
    Eigen::MatrixXd::Index minIndex;

    //cut-off is derived from (longest bond)
    //the cut-off is used to limit the search space
    double cut_off = me->distance_cut_off*me->distance_cut_off;
    //initialize some variables used for defining molecules in gathering process
    int molecule_start = 0, molecule_end = 0, molecule_start_previous = 0, initial_atom = 0;

    //initiallize variables used for loops
    int i_molecule_num = 0, i_atom = 0, j_atom = 0;

    //cross-reference gathering
    for (i_molecule_num = 0; i_molecule_num < me->solute_count; i_molecule_num++)
    {
        //define first atom, last atom and number of periodic copies of a solute molecule
        molecule_start = me->solute_molecules(0,i_molecule_num);
        molecule_end = me->solute_molecules(1,i_molecule_num);
        initial_atom = me->solute_molecules(3,i_molecule_num);
        framedata->solute_cog_count += (molecule_end - molecule_start) + 1;

        //the first atom in a solute molecule which is not the first solute molecule is 
        //gather with respect to all atoms of the previous solute molecule
        if (i_molecule_num > 0)
        { 
            shortestDistance = 1e20;
            currentDistance = 0;
            //reset the value for the shift applied
            min_shift.setZero();
            //molecule_start_previous is re-defined each time a solute molecule is gathered
            //molecule_start - 1 --> array index of the first atom of the current solute molecule
            //molecule_start - 2 --> array index of the last atom of the previous solute molecule
            //do a reverse look-up 
            for (j_atom = (molecule_start - 2); j_atom >= (molecule_start_previous - 1); j_atom--)
            {	
                currentDistance = (((((
                    (*grid)[i_molecule_num]
                ).cast<double>().array().colwise() * framedata->box_length.array()
                    ).matrix().colwise() + framedata->coordinates.col(molecule_start - 1)
                    ).colwise() - framedata->coordinates.col(j_atom)
                    ).cwiseAbs2().colwise().sum().minCoeff(&minIndex)
                    );
                //find the shift required to acquire the shortest distance
                if (currentDistance < shortestDistance) {
                    shortestDistance = currentDistance;
                    min_shift = ((*grid)[i_molecule_num]).col(minIndex).cast<int>();
                    if (shortestDistance < cut_off)
                    {
                        break;
                    }
                }
            }
            //shift the coordinates of the atom in the original frame
            framedata->coordinates.col(molecule_start - 1) += framedata->box_length.cwiseProduct(min_shift.cast<double>());
        }

        //center of geometry 
        //add the coordinates of the first atom to sum(x,y,z) for calculation of COG
        framedata->solute_cog_sum += framedata->coordinates.col(molecule_start - 1);

        //now that the first atom in the current solute is gathered correctly
        //the rest of the atoms will be gathered with respect to the first atom of the current solute.
        //gather solute molecule
        for (i_atom = molecule_start; i_atom < molecule_end; i_atom++)
        {	
            shortestDistance = 1e20;
            currentDistance = 0;			
            min_shift.setZero();
            //reverse search from coordinate closest to already-in-list atoms of current solute molecule
            for (j_atom = (i_atom - 1); j_atom >= (molecule_start - 1); j_atom--)
            {
                currentDistance = (((((
                    (*grid)[i_molecule_num]
                ).cast<double>().array().colwise() * framedata->box_length.array()
                    ).matrix().colwise() + framedata->coordinates.col(i_atom)
                    ).colwise() - framedata->coordinates.col(j_atom)
                    ).cwiseAbs2().colwise().sum().minCoeff(&minIndex)
                    );
                //find the shift required to acquire the shortest distance
                if (currentDistance < shortestDistance) {
                    shortestDistance = currentDistance;
                    min_shift = ((*grid)[i_molecule_num]).col(minIndex).cast<int>();
                    if (shortestDistance < cut_off)
                    {
                        break;
                    }
                }
            }

            //std::cout << i_atom << " " << j_atom << std::endl;
            //std::cout << framedata->coordinates.col(i_atom) << std::endl;
            //std::cout << min_shift << std::endl;
            //shift the coordinates of the atom i according to the frame-shift determined in preceeding code
            //std::cout << "  " << i_atom << std::endl;
            //std::cout << "  " << cut_off << std::endl;    
            //std::cout << "    " << framedata->coordinates.col(i_atom).x() << " "  << framedata->coordinates.col(i_atom).y() << " "  << framedata->coordinates.col(i_atom).z() << " " << std::endl;
            framedata->coordinates.col(i_atom) += framedata->box_length.cwiseProduct(min_shift.cast<double>());
            //std::cout << "    " << framedata->coordinates.col(i_atom).x() << " "  << framedata->coordinates.col(i_atom).y() << " "  << framedata->coordinates.col(i_atom).z() << " " << std::endl;
            //std::cout << framedata->coordinates.col(i_atom) << std::endl << std::endl;

            //add the new coordinate (x,y,z) values of the current atom to the COG calculation
            framedata->solute_cog_sum += framedata->coordinates.col(i_atom);
        }


        //correct for drifting molecule
        if (i_molecule_num > 0)
        { 
            shortestDistance = 1e20;
            currentDistance = 0;
            //reset the value for the shift applied
            min_shift.setZero();
            //molecule_start_previous is re-defined each time a solute molecule is gathered
            //molecule_start - 1 --> array index of the first atom of the current solute molecule
            //molecule_start - 2 --> array index of the last atom of the previous solute molecule
            //do a reverse look-up 
            for (j_atom = (molecule_start - 2); j_atom >= (molecule_start_previous - 1); j_atom--)
            {	
                currentDistance = (((((
                    (*grid)[i_molecule_num]
                ).cast<double>().array().colwise() * framedata->box_length.array()
                    ).matrix().colwise() + framedata->coordinates.col(initial_atom - 1)
                    ).colwise() - framedata->coordinates.col(j_atom)
                    ).cwiseAbs2().colwise().sum().minCoeff(&minIndex)
                    );
                //find the shift required to acquire the shortest distance
                if (currentDistance < shortestDistance) {
                    shortestDistance = currentDistance;
                    min_shift = ((*grid)[i_molecule_num]).col(minIndex).cast<int>();
                    //assuming the distance between individual molecules is at most
                    //3x the intermolecular cut-off... e.g. 3*0.2nm=0.6nm
                    if (shortestDistance < cut_off*3)
                    {
                        break;
                    }
                }
            }
            //shift the coordinates of the atom in the original frame
            if (!(min_shift.array() == 0).all())
            {
                framedata->coordinates.block(3, molecule_end - molecule_start + 1, 0, molecule_start -1) += framedata->box_length.cwiseProduct(min_shift.cast<double>());
                framedata->solute_cog_sum -= (framedata->box_length.cwiseProduct(min_shift.cast<double>()).array() * (molecule_end - molecule_start + 1)).matrix();
            }
        }

        //done gathering for a solute molecule
        //remember the first atom of the current solute before starting to gather the next solute molecule.
        //this is necessary for gathering the first atom of the next solute molecule 
        molecule_start_previous = molecule_start;
    }

    //if solute molecules are defined than calculate the COG for the frame and update the variable
    if (me->solute_count > 0)
    {
        framedata->solute_cog(0) = framedata->solute_cog_sum(0)/(framedata->solute_cog_count);
        framedata->solute_cog(1) = framedata->solute_cog_sum(1)/(framedata->solute_cog_count);
        framedata->solute_cog(2) = framedata->solute_cog_sum(2)/(framedata->solute_cog_count);
    }
    else
    {
        //use the coordinate of the first atom as center of geometry    
        framedata->solute_cog = framedata->coordinates.col(0);
    }
}

void Gather::SoluteCenterOfGeometry(Structs::FrameGeometric *framedata, Structs::InputParametersFrameout *me, std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grid)
{
    //std::cout << framedata->solute_cog_sum.array() << std::endl;
    //std::cout << framedata->solute_cog_count << std::endl;
    //std::cout << framedata->solute_cog << std::endl << std::endl;

    //these two varaible are required to ensure that gathering will always occured
    double shortestDistance = 1e20;
    double currentDistance = 0;

    Eigen::Vector3i min_shift(0,0,0);

    //use a for finding minimum distances in matrix/vector
    Eigen::MatrixXd::Index minIndex;

    //cut-off is derived from (longest bond)
    //the cut-off is used to limit the search space
    double cut_off = me->distance_cut_off*me->distance_cut_off;
    //initialize some variables used for defining molecules in gathering process
    int molecule_start = 0, molecule_end = 0, molecule_start_previous = 0, molecule_size = 0;

    //initiallize variables used for loops
    int i_molecule_num = 0, i_atom = 0, j_atom = 0;

    //cross-reference gathering
    for (i_molecule_num = 0; i_molecule_num < me->solute_cog_molecules.cols(); i_molecule_num++)
    {
        //define first atom, last atom and number of periodic copies of a solute molecule
        molecule_start = me->solute_cog_molecules(0,i_molecule_num);
        molecule_end = me->solute_cog_molecules(1,i_molecule_num);
        molecule_size = me->solute_cog_molecules(2,i_molecule_num);
        framedata->solute_cog_count += (molecule_end - molecule_start) + 1;
        //std::cout << molecule_start << " " << molecule_end << " " << molecule_size << std::endl;

        //the first atom in a molecule with respect to the center of geometry of solutes

        for (int i_set = (molecule_start-1); i_set < molecule_end; i_set+=molecule_size)
        { 
            //reset the value for the shift applied
            min_shift.setZero();
            //molecule_start_previous is re-defined each time a solute molecule is gathered

            (((((
                (*grid)[i_molecule_num]
            ).cast<double>().array().colwise() * framedata->box_length.array()
                ).matrix().colwise() + framedata->coordinates.col(i_set)
                ).colwise() - framedata->solute_cog
                ).cwiseAbs2().colwise().sum().minCoeff(&minIndex)
                );
            //find the shift required to acquire the shortest distance
            min_shift = ((*grid)[i_molecule_num]).col(minIndex).cast<int>();

            //shift the coordinates of the atom in the original frame
            framedata->coordinates.col(i_set) += framedata->box_length.cwiseProduct(min_shift.cast<double>());
            framedata->solute_cog_sum += framedata->coordinates.col(i_atom);

            //now that the first atom in the current solute is gathered correctly
            //the rest of the atoms will be gathered with respect to the first atom of the current solute.
            //gather solute molecule
            for (int i_atom = (i_set+1); i_atom < i_set+molecule_size; i_atom++)
            {	
                shortestDistance = 1e20;
                currentDistance = 0;			
                min_shift.setZero();
                //reverse search from coordinate closest to already-in-list atoms of current solute molecule
                for (j_atom = (i_atom - 1); j_atom >= (i_set); j_atom--)
                {
                    currentDistance = (((((
                        (*grid)[i_molecule_num]
                    ).cast<double>().array().colwise() * framedata->box_length.array()
                        ).matrix().colwise() + framedata->coordinates.col(i_atom)
                        ).colwise() - framedata->coordinates.col(j_atom)
                        ).cwiseAbs2().colwise().sum().minCoeff(&minIndex)
                        );
                    //find the shift required to acquire the shortest distance
                    if (currentDistance < shortestDistance) {
                        shortestDistance = currentDistance;
                        min_shift = ((*grid)[i_molecule_num]).col(minIndex).cast<int>();
                        if (shortestDistance < cut_off)
                        {
                            break;
                        }
                    }
                }

                //shift the coordinates of the atom i according to the frame-shift determined in preceeding code
                framedata->coordinates.col(i_atom) += framedata->box_length.cwiseProduct(min_shift.cast<double>());
                framedata->solute_cog_sum += framedata->coordinates.col(i_atom);
            }
        }
    }

    //if solute molecules are defined than calculate the COG for the frame and update the variable
    if (me->solute_cog_molecules.cols() > 0)
    {
        framedata->solute_cog(0) = framedata->solute_cog_sum(0)/(framedata->solute_cog_count);
        framedata->solute_cog(1) = framedata->solute_cog_sum(1)/(framedata->solute_cog_count);
        framedata->solute_cog(2) = framedata->solute_cog_sum(2)/(framedata->solute_cog_count);
    }
}

void Gather::IonsCenterOfGeometry(Structs::FrameGeometric *framedata, Structs::InputParametersFrameout *me, std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grid)
{
    //std::cout << framedata->solute_cog_sum.array() << std::endl;
    //std::cout << framedata->solute_cog_count << std::endl;
    //std::cout << framedata->solute_cog << std::endl << std::endl;

    //these two varaible are required to ensure that gathering will always occured
    double shortestDistance = 1e20;
    double currentDistance = 0;

    Eigen::Vector3i min_shift(0,0,0);

    //use a for finding minimum distances in matrix/vector
    Eigen::MatrixXd::Index minIndex;

    //cut-off is derived from (longest bond)
    //the cut-off is used to limit the search space
    double cut_off = me->distance_cut_off*me->distance_cut_off;
    //initialize some variables used for defining molecules in gathering process
    int molecule_start = 0, molecule_end = 0, molecule_start_previous = 0, molecule_size = 0;

    //initiallize variables used for loops
    int i_molecule_num = 0, i_atom = 0, j_atom = 0;

    //cross-reference gathering
    for (i_molecule_num = 0; i_molecule_num < me->ion_molecules.cols(); i_molecule_num++)
    {
        //define first atom, last atom and number of periodic copies of a solute molecule
        molecule_start = me->ion_molecules(0,i_molecule_num);
        molecule_end = me->ion_molecules(1,i_molecule_num);
        molecule_size = me->ion_molecules(2,i_molecule_num);
        framedata->solute_cog_count += (molecule_end - molecule_start) + 1;
        //std::cout << molecule_start << " " << molecule_end << " " << molecule_size << std::endl;

        //the first atom in a molecule with respect to the center of geometry of solutes

        for (int i_set = (molecule_start-1); i_set < molecule_end; i_set+=molecule_size)
        { 
            //reset the value for the shift applied
            min_shift.setZero();
            //molecule_start_previous is re-defined each time a solute molecule is gathered

            (((((
                (*grid)[i_molecule_num]
            ).cast<double>().array().colwise() * framedata->box_length.array()
                ).matrix().colwise() + framedata->coordinates.col(i_set)
                ).colwise() - framedata->solute_cog
                ).cwiseAbs2().colwise().sum().minCoeff(&minIndex)
                );
            //find the shift required to acquire the shortest distance
            min_shift = ((*grid)[i_molecule_num]).col(minIndex).cast<int>();

            //shift the coordinates of the atom in the original frame
            framedata->coordinates.col(i_set) += framedata->box_length.cwiseProduct(min_shift.cast<double>());
            framedata->solute_cog_sum += framedata->coordinates.col(i_set);

            //now that the first atom in the current solute is gathered correctly
            //the rest of the atoms will be gathered with respect to the first atom of the current solute.
            //gather solute molecule
            for (int i_atom = (i_set+1); i_atom < i_set+molecule_size; i_atom++)
            {	
                shortestDistance = 1e20;
                currentDistance = 0;			
                min_shift.setZero();
                //reverse search from coordinate closest to already-in-list atoms of current solute molecule
                for (j_atom = (i_atom - 1); j_atom >= (i_set); j_atom--)
                {
                    currentDistance = (((((
                        (*grid)[i_molecule_num]
                    ).cast<double>().array().colwise() * framedata->box_length.array()
                        ).matrix().colwise() + framedata->coordinates.col(i_atom)
                        ).colwise() - framedata->coordinates.col(j_atom)
                        ).cwiseAbs2().colwise().sum().minCoeff(&minIndex)
                        );
                    //find the shift required to acquire the shortest distance
                    if (currentDistance < shortestDistance) {
                        shortestDistance = currentDistance;
                        min_shift = ((*grid)[i_molecule_num]).col(minIndex).cast<int>();
                        if (shortestDistance < cut_off)
                        {
                            break;
                        }
                    }
                }

                //shift the coordinates of the atom i according to the frame-shift determined in preceeding code
                framedata->coordinates.col(i_atom) += framedata->box_length.cwiseProduct(min_shift.cast<double>());
                framedata->solute_cog_sum += framedata->coordinates.col(i_atom);
            }
        }
    }

    //if solute molecules are defined than calculate the COG for the frame and update the variable
    if (me->ion_molecules.cols() > 0)
    {
        framedata->solute_cog(0) = framedata->solute_cog_sum(0)/(framedata->solute_cog_count);
        framedata->solute_cog(1) = framedata->solute_cog_sum(1)/(framedata->solute_cog_count);
        framedata->solute_cog(2) = framedata->solute_cog_sum(2)/(framedata->solute_cog_count);
    }
}

void Gather::Solvent(Structs::FrameGeometric *framedata, Structs::InputParametersFrameout *me, std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grid)
{
    //std::cout << framedata->solute_cog_sum.array() << std::endl;
    //std::cout << framedata->solute_cog_count << std::endl;
    //std::cout << framedata->solute_cog << std::endl << std::endl;

    //these two varaible are required to ensure that gathering will always occured
    double shortestDistance = 1e20;
    double currentDistance = 0;

    Eigen::Vector3i min_shift(0,0,0);

    //use a for finding minimum distances in matrix/vector
    Eigen::MatrixXd::Index minIndex;

    //cut-off is derived from (longest bond)
    //the cut-off is used to limit the search space
    double cut_off = me->distance_cut_off*me->distance_cut_off;
    //initialize some variables used for defining molecules in gathering process
    int molecule_start = 0, molecule_end = 0, molecule_start_previous = 0, molecule_size = 0;

    //initiallize variables used for loops
    int i_molecule_num = 0, i_atom = 0, j_atom = 0;

    //for center of geometry of all solute molecules calculation
    //at the end of the frame this holds sum(x,y,z) 
    Eigen::Vector3d coordinates_sum(0,0,0);

    //cross-reference gathering
    for (i_molecule_num = 0; i_molecule_num < me->solvent_molecules.cols(); i_molecule_num++)
    {
        //define first atom, last atom and number of periodic copies of a solute molecule
        molecule_start = me->solvent_molecules(0,i_molecule_num);
        molecule_end = me->solvent_molecules(1,i_molecule_num);
        molecule_size = me->solvent_molecules(2,i_molecule_num);
        //std::cout << molecule_start << " " << molecule_end << " " << molecule_size << std::endl;

        //the first atom in a molecule with respect to the center of geometry of solutes

        for (int i_set = (molecule_start-1); i_set < molecule_end; i_set+=molecule_size)
        { 
            //reset the value for the shift applied
            min_shift.setZero();
            //molecule_start_previous is re-defined each time a solute molecule is gathered

            (((((
                (*grid)[i_molecule_num]
            ).cast<double>().array().colwise() * framedata->box_length.array()
                ).matrix().colwise() + framedata->coordinates.col(i_set)
                ).colwise() - framedata->solute_cog
                ).cwiseAbs2().colwise().sum().minCoeff(&minIndex)
                );
            //find the shift required to acquire the shortest distance
            min_shift = ((*grid)[i_molecule_num]).col(minIndex).cast<int>();

            //shift the coordinates of the atom in the original frame
            framedata->coordinates.col(i_set) += framedata->box_length.cwiseProduct(min_shift.cast<double>());

            //now that the first atom in the current solute is gathered correctly
            //the rest of the atoms will be gathered with respect to the first atom of the current solute.
            //gather solute molecule
            for (int i_atom = (i_set+1); i_atom < i_set+molecule_size; i_atom++)
            {	
                shortestDistance = 1e20;
                currentDistance = 0;			
                min_shift.setZero();
                //reverse search from coordinate closest to already-in-list atoms of current solute molecule
                for (j_atom = (i_atom - 1); j_atom >= (i_set); j_atom--)
                {
                    currentDistance = (((((
                        (*grid)[i_molecule_num]
                    ).cast<double>().array().colwise() * framedata->box_length.array()
                        ).matrix().colwise() + framedata->coordinates.col(i_atom)
                        ).colwise() - framedata->coordinates.col(j_atom)
                        ).cwiseAbs2().colwise().sum().minCoeff(&minIndex)
                        );
                    //find the shift required to acquire the shortest distance
                    if (currentDistance < shortestDistance) {
                        shortestDistance = currentDistance;
                        min_shift = ((*grid)[i_molecule_num]).col(minIndex).cast<int>();
                        if (shortestDistance < cut_off)
                        {
                            break;
                        }
                    }
                }

                //shift the coordinates of the atom i according to the frame-shift determined in preceeding code
                framedata->coordinates.col(i_atom) += framedata->box_length.cwiseProduct(min_shift.cast<double>());
            }
        }
    }
}

void Gather::Solvent(Structs::FrameGeometric *framedata, Structs::InputParametersFrameout *me, std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grid, double solvent_sphere_radius)
{
    //std::cout << framedata->solute_cog_sum.array() << std::endl;
    //std::cout << framedata->solute_cog_count << std::endl;
    //std::cout << framedata->solute_cog << std::endl << std::endl;

    //where to stop applying correction
    int last_atom = 0;

    //apply correction to the last of the solutes
    if ((me->solute_molecules.cols() > 0 && me->solute_cog_molecules.cols() > 0))
    {
        last_atom = me->solute_cog_molecules(1,me->solute_cog_molecules.cols()-1);
    }
    //apply correction to the essential solutes
    else if ((me->solute_molecules.cols() > 0))
    {
        last_atom = me->solute_molecules(1,me->solute_molecules.cols()-1);
    }
    //if ions are not present, apply correction based on atoms up to the first solvent
    else
    {
        last_atom = me->solvent_molecules(0,0) - 1;
    }

    //get the x, y or z value in solutes and ions
    double cut_off_radius_solvent = (framedata->coordinates.leftCols(last_atom) - framedata->solute_cog).cwiseAbs2().colwise().sum().maxCoeff();
    cut_off_radius_solvent = sqrt(cut_off_radius_solvent);
    cut_off_radius_solvent += solvent_sphere_radius;
    cut_off_radius_solvent = cut_off_radius_solvent*cut_off_radius_solvent;
    bool solvent_include = true;

    //these two varaible are required to ensure that gathering will always occured
    double shortestDistance = 1e20;
    double currentDistance = 0;

    Eigen::Vector3i min_shift(0,0,0);

    //use a for finding minimum distances in matrix/vector
    Eigen::MatrixXd::Index minIndex;

    //cut-off is derived from (longest bond)
    //the cut-off is used to limit the search space
    double cut_off = me->distance_cut_off*me->distance_cut_off;
    //initialize some variables used for defining molecules in gathering process
    int molecule_start = 0, molecule_end = 0, molecule_start_previous = 0, molecule_size = 0;

    //initiallize variables used for loops
    int i_molecule_num = 0, i_atom = 0, j_atom = 0;

    //for center of geometry of all solute molecules calculation
    //at the end of the frame this holds sum(x,y,z) 
    Eigen::Vector3d coordinates_sum(0,0,0);

    //cross-reference gathering
    for (i_molecule_num = 0; i_molecule_num < me->solvent_molecules.cols(); i_molecule_num++)
    {
        //define first atom, last atom and number of periodic copies of a solute molecule
        molecule_start = me->solvent_molecules(0,i_molecule_num);
        molecule_end = me->solvent_molecules(1,i_molecule_num);
        molecule_size = me->solvent_molecules(2,i_molecule_num);
        framedata->solute_cog_count += (molecule_end - molecule_start) + 1;
        //std::cout << molecule_start << " " << molecule_end << " " << molecule_size << std::endl;

        //the first atom in a molecule with respect to the center of geometry of solutes

        for (int i_set = (molecule_start-1); i_set < molecule_end; i_set+=molecule_size)
        { 
            solvent_include = true;
            //reset the value for the shift applied
            min_shift.setZero();
            //molecule_start_previous is re-defined each time a solute molecule is gathered

            currentDistance = (((((
                (*grid)[i_molecule_num]
            ).cast<double>().array().colwise() * framedata->box_length.array()
                ).matrix().colwise() + framedata->coordinates.col(i_set)
                ).colwise() - framedata->solute_cog
                ).cwiseAbs2().colwise().sum().minCoeff(&minIndex)
                );
            //find the shift required to acquire the shortest distance
            min_shift = ((*grid)[i_molecule_num]).col(minIndex).cast<int>();

            //shift the coordinates of the atom in the original frame
            framedata->coordinates.col(i_set) += framedata->box_length.cwiseProduct(min_shift.cast<double>());

            if ((currentDistance) > cut_off_radius_solvent)
            {
                //framedata->coordinates.col(i_set).setZero();
                solvent_include = false;
            }

            //std::cout << "\t\t" << currentDistance << " " << cut_off_radius_solvent << " " << ((currentDistance) > cut_off_radius_solvent) << solvent_include << std::endl;

            //now that the first atom in the current solute is gathered correctly
            //the rest of the atoms will be gathered with respect to the first atom of the current solute.
            //gather solute molecule
            for (int i_atom = (i_set+1); i_atom < i_set+molecule_size; i_atom++)
            {	
                shortestDistance = 1e20;
                currentDistance = 0;			
                min_shift.setZero();
                //reverse search from coordinate closest to already-in-list atoms of current solute molecule
                for (j_atom = (i_atom - 1); j_atom >= (i_set); j_atom--)
                {
                    currentDistance = (((((
                        (*grid)[i_molecule_num]
                    ).cast<double>().array().colwise() * framedata->box_length.array()
                        ).matrix().colwise() + framedata->coordinates.col(i_atom)
                        ).colwise() - framedata->coordinates.col(j_atom)
                        ).cwiseAbs2().colwise().sum().minCoeff(&minIndex)
                        );
                    //find the shift required to acquire the shortest distance
                    if (currentDistance < shortestDistance) {
                        shortestDistance = currentDistance;
                        min_shift = ((*grid)[i_molecule_num]).col(minIndex).cast<int>();
                        if (shortestDistance < cut_off)
                        {
                            break;
                        }
                    }
                }

                //shift the coordinates of the atom i according to the frame-shift determined in preceeding code
                framedata->coordinates.col(i_atom) += framedata->box_length.cwiseProduct(min_shift.cast<double>());

                if (!solvent_include)
                {
                    //framedata->coordinates.col(i_atom).setZero();
                }
            }

            if (!solvent_include)
            {
                framedata->coordinates.block(0,i_set,3, 3).colwise() += framedata->box_length;
            }
            //std::cout << framedata->coordinates.block(0,i_set,3, 3) << std::endl;
        }
    }
}