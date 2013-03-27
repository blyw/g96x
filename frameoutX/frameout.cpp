#define DEBUG
#define POSIX

//jacobi
#define j_abs(x) (temp_abs=x)*((temp_abs>=0)*2-1)
#define j_a(i,j) a[(j-1)*3+i-1]
#define j_v(i,j) v[(j-1)*3+i-1]
#define j_b(i) b[i-1]
#define j_z(i) z[i-1]
#define j_d(i) d[i-1]

#include <iostream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <chrono>
#include <thread>
#include <vector>
#include <future>
#include <cmath>
#include "gzstream.h"
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <algorithm> 
#include <numeric>
#include <complex>

//a single frame containing coordinates data
struct frame {
    long timestep;
    double time;
    std::vector<std::string> prefix;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    double solute_cog_x;
    double solute_cog_y;
    double solute_cog_z;
    int boxtype;
    double box_length_x;
    double box_length_y;
    double box_length_z;
    double box_angle_x;
    double box_angle_y;
    double box_angle_z;
    double box_3_x;
    double box_3_y;
    double box_3_z;
    double box_4_x;
    double box_4_y;
    double box_4_z;
    int init_shift_x;
    int init_shift_y;
    int init_shift_z;
};

//struct for holding parameters needed by the program
struct params {
    //hardware parameters
    int num_thread;
    int num_frame_per_thread;
    int num_thread_real;

    //output parameters
    int output_fragment_size;
    int output_fragment_skipframes;
    double output_fragment_skiptime;
    std::string outformat;
    std::string outfilename;

    //topology parameters - solutes
    int solute_count;
    std::vector<std::vector<int>> solute_molecules;
    std::vector<std::vector<int>> solutes_cog_molecules;

    //topology parameters - solvent
    int solvent_atoms_first;
    int solvent_atoms_last;
    int solvent_dimension_expansion;
    int solvent_size;
    bool solvent_skip;

    //input paramters
    std::string informat;
    std::string trc_reference;
    std::vector<std::string> input_files;

    //frameoutX parameters
    double ref_coords[3];
    bool cog_write;
    bool cog_correction;
    bool rotation_correction;

    //shared parameters
    int atomrecords;    
    double distance_cut_off;
    int verbosity;

    //dmovX parameters essential parameters
    std::vector<std::vector<int>> dmov_atomgroup;
    std::vector<std::vector<double>> dmov_calc_previous_last_eigenvalvec;
    //dmovX parameters data parameters
    std::vector<std::vector<int>> dmov_angles;
    std::vector<std::vector<int>> dmov_dihedral_angles;
    std::vector<std::vector<int>> dmov_atomgroup_inclusion_type;
    std::vector<double> dmov_write_previous_dihedral_angles;
    std::vector<double> dmov_write_previous_angles;

    //Jacobi method parameters
    double jacobi_max_iteration;
    double jacobi_vector_correction_cutoff;
};

//code checked 20130122: CORRECT!
//simple gathering of a specified atom in frame with respect to a given coordinate
void FirstAtomBasedBoxShifter(frame *framedata, int atom_number, params *me) {
    double coords[3] = { 0, 0, 0 };
    double distance_shortest = 1E20;
    int min_shift[3] = { 0, 0, 0 };
    double *ref_coords = me->ref_coords;

    //dimension expansion factor
    int periodic_copies = 2;
    if (me->solute_count > 0)
    {
        periodic_copies = (me->solute_molecules[0][2]+1)*2;
    }

    //shift the coordinates of the specified atom in the original framedata based on the previous frame shift
    if (!(framedata->init_shift_x == 0 && framedata->init_shift_y == 0 && framedata->init_shift_z == 0))
    {
        for (int i = 0; i < framedata->x.size(); i++)
        {        
            framedata->x[i] = framedata->x[i] + (framedata->init_shift_x * framedata->box_length_x);
        }
        for (int i = 0; i < framedata->y.size(); i++)
        {        
            framedata->y[i] = framedata->y[i] + (framedata->init_shift_y * framedata->box_length_y);
        }
        for (int i = 0; i < framedata->z.size(); i++)
        {
            framedata->z[i] = framedata->z[i] + (framedata->init_shift_z * framedata->box_length_z);
        }
    }

    //gather the first atom fo the frame
    for (int x = 0-periodic_copies; x <= periodic_copies; x++)
    {
        for (int y = 0-periodic_copies; y <= periodic_copies; y++)
        {
            for (int z = 0-periodic_copies; z <= periodic_copies; z++)
            {
                double distance =
                    (((framedata->x[atom_number] + x * framedata->box_length_x) - ref_coords[0]) * ((framedata->x[atom_number] + x * framedata->box_length_x) - ref_coords[0])) +
                    (((framedata->y[atom_number] + y * framedata->box_length_y) - ref_coords[1]) * ((framedata->y[atom_number] + y * framedata->box_length_y) - ref_coords[1])) +
                    (((framedata->z[atom_number] + z * framedata->box_length_z) - ref_coords[2]) * ((framedata->z[atom_number] + z * framedata->box_length_z) - ref_coords[2]));
                if (distance < distance_shortest)
                {
                    distance_shortest = distance;
                    min_shift[0] = x;
                    min_shift[1] = y;
                    min_shift[2] = z;
                }
            }
        }
    }

    //shift again if first atom was shifted
    if (!( min_shift[0] == 0 && min_shift[1] == 0 && min_shift[2] == 0))
    {
        for (int i = 0; i < framedata->x.size(); i++)
        {        
            framedata->x[i] = framedata->x[i] + (min_shift[0] * framedata->box_length_x);
        }
        for (int i = 0; i < framedata->y.size(); i++)
        {        
            framedata->y[i] = framedata->y[i] + (min_shift[1] * framedata->box_length_y);
        }
        for (int i = 0; i < framedata->z.size(); i++)
        {
            framedata->z[i] = framedata->z[i] + (min_shift[2] * framedata->box_length_z);
        }
    }

    //maybe useful for debugging?
    framedata->init_shift_x = min_shift[0];
    framedata->init_shift_y = min_shift[1];
    framedata->init_shift_z = min_shift[2];
}

//code checked 20130122: CORRECT!
//solute gathering i.e. cluster algorithm (nearest neighbour)
void SoluteGatherer(frame *framedata, int frameId, params *me){ // std::vector<std::vector<int>> solute_molecules, int solute_count, int frameId, int periodic_copies) {	
    //init variables
    double distance_shortest;
    int min_shift[3] = { 0, 0, 0 };
    //cut-off is derived from (longest bond)^2
    double cut_off = me->distance_cut_off;
    int molecule_start, molecule_end, molecule_start_previous, periodic_copies;

    //for center of geometry of all solute molecules calculation
    double x_sum = 0;
    double y_sum = 0;
    double z_sum = 0;

    //cross-reference gathering
    for (int s = 0; s < me->solute_count; s++)
    {
        //define first atom, last atom and number of periodic copies of a solute molecule
        molecule_start = me->solute_molecules[s][0];
        molecule_end = me->solute_molecules[s][1];
        periodic_copies = me->solute_molecules[s][2];

        //the first atom in a solute molecule which is not the first solute molecule is 
        //gather with respect to all atoms of the previous solute molecule
        if (s > 0)
        { 
            distance_shortest = 1E20;
            min_shift[0] = 0;
            min_shift[1] = 0;
            min_shift[2] = 0;
            //molecule_start_previous is re-defined each time a solute molecule is gathered
            //molecule_start - 1 --> array index of the first atom of the current solute molecule
            //molecule_start - 2 --> array index of the last atom of the previous solute molecule
            for (int i = (molecule_start - 2); i >= (molecule_start_previous - 1); i--)
            {	
                //gather
                for (int x = 0-periodic_copies; x <= periodic_copies; x++)
                {
                    for (int y = 0-periodic_copies; y <= periodic_copies; y++)
                    {
                        for (int z = 0-periodic_copies; z <= periodic_copies; z++)
                        {
                            //defined cutt-off use to prevent unnecessary looping through the full list
                            if (distance_shortest <= cut_off)
                            {
                                break;
                            }

                            double distance =
                                (((framedata->x[molecule_start - 1] + x * framedata->box_length_x) - framedata->x[i]) * ((framedata->x[molecule_start - 1] + x * framedata->box_length_x) - framedata->x[i])) +
                                (((framedata->y[molecule_start - 1] + y * framedata->box_length_y) - framedata->y[i]) * ((framedata->y[molecule_start - 1] + y * framedata->box_length_y) - framedata->y[i])) +
                                (((framedata->z[molecule_start - 1] + z * framedata->box_length_z) - framedata->z[i]) * ((framedata->z[molecule_start - 1] + z * framedata->box_length_z) - framedata->z[i]));
                            if (distance < distance_shortest)
                            {
                                distance_shortest = distance;
                                min_shift[0] = x;
                                min_shift[1] = y;
                                min_shift[2] = z;
                            }
                        }
                    }
                }
            }
            //shift the coordinates of the atom in the original frame
            framedata->x[molecule_start - 1] = framedata->x[molecule_start - 1] + (min_shift[0] * framedata->box_length_x);
            framedata->y[molecule_start - 1] = framedata->y[molecule_start - 1] + (min_shift[1] * framedata->box_length_y);
            framedata->z[molecule_start - 1] = framedata->z[molecule_start - 1] + (min_shift[2] * framedata->box_length_z);
        }

        //dimension sum for cog 
        x_sum += framedata->x[molecule_start - 1];
        y_sum += framedata->y[molecule_start - 1];
        z_sum += framedata->z[molecule_start - 1];

        //double xval, yval, zval;
        //gather solute molecule
        for (int i = molecule_start; i < molecule_end; i++)
        {				
            //atom dependent variables that should be reset for each atom
            distance_shortest = 1E20;
            min_shift[0] = 0;
            min_shift[1] = 0;
            min_shift[2] = 0;
            //reverse search from coordinate closest to already-in-list atoms of current solute molecule
            for (int j = (i - 1); j >= (molecule_start - 1); j--)
            {
                //gather
                for (int x = 0-periodic_copies; x <= periodic_copies; x++)
                {
                    for (int y = 0-periodic_copies; y <= periodic_copies; y++)
                    {
                        for (int z = 0-periodic_copies; z <= periodic_copies; z++)
                        {
                            //defined cutt-off use to prevent unnecessary looping through the full list
                            if (distance_shortest <= cut_off)
                            {
                                break;
                            }

                            double distance =
                                (((framedata->x[i] + x * framedata->box_length_x) - framedata->x[j]) * ((framedata->x[i] + x * framedata->box_length_x) - framedata->x[j])) +
                                (((framedata->y[i] + y * framedata->box_length_y) - framedata->y[j]) * ((framedata->y[i] + y * framedata->box_length_y) - framedata->y[j])) +
                                (((framedata->z[i] + z * framedata->box_length_z) - framedata->z[j]) * ((framedata->z[i] + z * framedata->box_length_z) - framedata->z[j]));
                            if (distance < distance_shortest)
                            {
                                distance_shortest = distance;
                                min_shift[0] = x;
                                min_shift[1] = y;
                                min_shift[2] = z;
                            }
                        }
                    }
                }
            }

            //shift the coordinates of the atom in the original frame	
            framedata->x[i] = framedata->x[i] + (min_shift[0] * framedata->box_length_x);
            framedata->y[i] = framedata->y[i] + (min_shift[1] * framedata->box_length_y);
            framedata->z[i] = framedata->z[i] + (min_shift[2] * framedata->box_length_z);

            //dimension sum for cog 
            x_sum += framedata->x[i];
            y_sum += framedata->y[i];
            z_sum += framedata->z[i];
        }
        //done gathering for a solute molecule
        //set reference for gathering next solute molecule
        molecule_start_previous = molecule_start;
    }

    if (me->solute_count > 0)
    {
        //use this center of geometry for the entire frame
        framedata->solute_cog_x = x_sum / molecule_end;
        framedata->solute_cog_y = y_sum / molecule_end;
        framedata->solute_cog_z = z_sum / molecule_end;
    }
    else
    {
        //use this center of geometry for the entire frame
        framedata->solute_cog_x = framedata->x[0];
        framedata->solute_cog_y = framedata->y[0];
        framedata->solute_cog_z = framedata->z[0];
    }
    //std::cout << frameId << std::endl;
}

//code checked 20130122: CORRECT!
//gather solvent molecules with respect to COG of solutes
void COGGatherer(frame *framedata, int first_atom, int last_atom, int molecule_size, int periodic_copies) {
    double distance_shortest = 1E20;
    int min_shift[3] = { 0, 0, 0 };
    //cut-off is derived from (longest bond)^2
    double cut_off = 0.3 * 0.3;

    for (int i = first_atom - 1; i < last_atom; i+=molecule_size)
    {
        distance_shortest = 1E20;
        min_shift[0] = 0;
        min_shift[1] = 0;
        min_shift[2] = 0;

        //gather the first atom of the solvent molecule
        for (int x = 0-periodic_copies; x <= periodic_copies; x++)
        {
            for (int y = 0-periodic_copies; y <= periodic_copies; y++)
            {
                for (int z = 0-periodic_copies; z <= periodic_copies; z++)
                {                    
                    double distance =
                        (((framedata->x[i] + x * framedata->box_length_x) - framedata->solute_cog_x) * ((framedata->x[i] + x * framedata->box_length_x) - framedata->solute_cog_x)) +
                        (((framedata->y[i] + y * framedata->box_length_y) - framedata->solute_cog_y) * ((framedata->y[i] + y * framedata->box_length_y) - framedata->solute_cog_y)) +
                        (((framedata->z[i] + z * framedata->box_length_z) - framedata->solute_cog_z) * ((framedata->z[i] + z * framedata->box_length_z) - framedata->solute_cog_z));
                    if (distance < distance_shortest)
                    {
                        distance_shortest = distance;
                        min_shift[0] = x;
                        min_shift[1] = y;
                        min_shift[2] = z;
                    }
                }
            }
        }
        //shift the coordinates of the first atom in the original framedata
        framedata->x[i] = framedata->x[i] + (min_shift[0] * framedata->box_length_x);
        framedata->y[i] = framedata->y[i] + (min_shift[1] * framedata->box_length_y);
        framedata->z[i] = framedata->z[i] + (min_shift[2] * framedata->box_length_z);

        ////gather other atoms of the solvent molecule
        for (int j = (i + 1); j < (i + molecule_size); j++)
        {
            //reset for the each additional atom in solvent molecule
            distance_shortest = 1E20;
            min_shift[0] = 0;
            min_shift[1] = 0;
            min_shift[2] = 0;

            for (int k = i; k < j; k++)
            {   
                for (int x = 0-periodic_copies; x <= periodic_copies; x++)
                {
                    for (int y = 0-periodic_copies; y <= periodic_copies; y++)
                    {
                        for (int z = 0-periodic_copies; z < periodic_copies; z++)
                        {
                            //defined cutt-off use to prevent unnecessary looping through the full list
                            if (distance_shortest <= cut_off)
                            {
                                break;
                            }

                            double distance =
                                (((framedata->x[j] + x * framedata->box_length_x) - framedata->x[k]) * ((framedata->x[j] + x * framedata->box_length_x) - framedata->x[k])) +
                                (((framedata->y[j] + y * framedata->box_length_y) - framedata->y[k]) * ((framedata->y[j] + y * framedata->box_length_y) - framedata->y[k])) +
                                (((framedata->z[j] + z * framedata->box_length_z) - framedata->z[k]) * ((framedata->z[j] + z * framedata->box_length_z) - framedata->z[k]));
                            if (distance < distance_shortest)
                            {
                                distance_shortest = distance;
                                min_shift[0] = x;
                                min_shift[1] = y;
                                min_shift[2] = z;
                            }
                        }
                    }
                }
            }
            //	//shift atom correctly
            framedata->x[j] = framedata->x[j] + (min_shift[0] * framedata->box_length_x);
            framedata->y[j] = framedata->y[j] + (min_shift[1] * framedata->box_length_y);
            framedata->z[j] = framedata->z[j] + (min_shift[2] * framedata->box_length_z);
        }
    }
}

//sort function
bool sortFunction (std::vector<double> i, std::vector<double> j) { 
    return (i[0] > j[0]); 
}

//Jacobi method
int Jacobi (double *a , double *d , double *v , int maxsw)
    /* Jacobi's method for diagonalysing matrix a[0..ndim**2-1], returning
    the eigenvalue set d[0..ndim-1] and eigenvector matrix v[0..ndim**2-1].
    Returns the number of required rotations or -1 if unsuccessful. */
{
    int i = 0, j = 0, k = 0, l = 0, nrot = 0, nsweep = 0;
    double msize = 3*3-1;
    double b[3] = { 0, 0, 0 }, z[3] = { 0, 0, 0 };
    double x = 0, temp_abs = 0, sum = 0, tresh = 0;
    double g = 0, h = 0, s = 0, tau = 0, t = 0, c = 0, theta = 0;
    double *vp[3],*xp;
    double vtmp[3*3] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    /* initialise v to identity */

    for ( i = 1 ; i <= 3 ; i++ )
        for ( j = 1 ; j <= 3 ; j++ )
            j_v(i,j) = (i==j);

    /* initialise b and d with diagonal elements and z with zero */

    for ( i = 1 ; i <= 3 ; i++ ) {
        j_d(i) = j_b(i) = j_a(i,i);
        j_z(i) = 0.0;
    }

    for ( nsweep=1 ; nsweep <= maxsw ; nsweep++ ) {

        /* sum is set to the half sum of off diagonal elements absolute value 
        and check for convergence ( OK if underflow rounded to zero ) */

        sum = 0.0;
        for ( i = 1 ; i <= 3-1 ; i++ )
            for ( j = i+1 ; j <= 3 ; j++ )
                sum += j_abs(j_a(i,j));
        if ( sum == 0.0 ) {

            /* restore a matrix to its initial value */
            for ( k=1 ; k<=3-1 ; k++ )
                for ( l=k; l<=3 ; l++ )
                    j_a(k,l) = j_a(l,k);

            /* eigenvalue reordering */
            for ( k=1 ; k<=3 ; k++ )
                vp[k-1] = v+3*(k-1);

            for ( k=1 ; k<=3-1 ; k++ )
                for ( l=k+1; l<=3 ; l++ )
                    if ( j_d(l) < j_d(k) ) {
                        x = j_d(k);
                        j_d(k) = j_d(l);
                        j_d(l) = x;
                        xp = vp[k-1];
                        vp[k-1] = vp[l-1];
                        vp[l-1] = xp;
                    }
                    for ( k=1 ; k<=3 ; k++ ) {
                        xp = vp[k-1];
                        for ( l=0 ; l<3 ; l++ )
                            vtmp[3*(k-1)+l] = *(xp+l);
                    }
                    for ( k=1 ; k<=3*3 ; k++ )
                        v[k-1] = vtmp[k-1];

                    return nrot;
        }
        /* for the first three sweeps tresh <> 0 */

        if ( nsweep <= 3 ) 
            tresh = 0.2*sum/(3*3);
        else
            tresh = 0.0;

        /* perform the systematic sweep */

        for ( i = 1 ; i <= 3-1 ; i++ )
            for ( j = i+1 ; j <= 3 ; j++ ) {

                /* perform the rotations, sweeping all (i,j) off-diagonal pairs of the upper
                triangle */

                g = 100.0*j_abs(j_a(i,j));

                /* after four sweeps, rotation is skipped if off-diagonal element is small */

                if ( nsweep > 4 && abs(j_d(i))+g == abs(j_d(i)) &&
                    j_abs(j_d(j))+g == j_abs(j_d(j)) )
                    j_a(i,j) = 0.0;
                else {

                    /* rotation is then performed only if the off-diagonal element is > tresh */

                    if ( j_abs(j_a(i,j)) > tresh ) {

                        /* rotation */

                        /* t calculation: */
                        h = j_d(j)-j_d(i);
                        if ( j_abs(h)+g == j_abs(h) )
                            /* t = 1/(2*theta), to avoid numerical errors */
                                t = j_a(i,j)/(2*h);
                        else {
                            /* t = sgn(theta)/(abs(theta)+sqrt(1.0+theta*theta) */
                            theta = 0.5*h/j_a(i,j);
                            t = 1.0/(j_abs(theta)+sqrt(1.0+theta*theta))
                                * (2*(theta>0)-1);
                        }

                        c = 1.0/sqrt(1.0+t*t);
                        s = t*c;
                        tau = s/(1.0+c);
                        h = t*j_a(i,j);

                        j_z(i) -= h;
                        j_z(j) += h;
                        j_d(i) -= h;
                        j_d(j) += h;

                        j_a(i,j) = 0.0;

                        for ( k=1 ; k <= i-1 ; k++ ) {
                            g = j_a(k,i);
                            h = j_a(k,j);
                            j_a(k,i) = g-s*(h+g*tau);
                            j_a(k,j) = h+s*(g-h*tau);
                        }
                        for ( k=i+1 ; k <= j-1 ; k++ ) {
                            g = j_a(i,k);
                            h = j_a(k,j);
                            j_a(i,k) = g-s*(h+g*tau);
                            j_a(k,j) = h+s*(g-h*tau);
                        }
                        for ( k=j+1 ; k <= 3 ; k++ ) {
                            g = j_a(i,k);
                            h = j_a(j,k);
                            j_a(i,k) = g-s*(h+g*tau);
                            j_a(j,k) = h+s*(g-h*tau);
                        }

                        for ( k=1 ; k <= 3 ; k++ ) {
                            g = j_v(k,i);
                            h = j_v(k,j);
                            j_v(k,i) = g-s*(h+g*tau);
                            j_v(k,j) = h+s*(g-h*tau);
                        }

                        nrot++;

                    } /* if ... > tresh */
                } /* if element not too small */
            } /* rotation (i,j) done */
            for ( i=1 ; i<=3 ; i++ ) {
                j_d(i) = ( j_b(i) += j_z(i) );
                j_z(i) = 0.0;
            }
    } /* sweep done */

    /*printf ("JACOBI : error, not succeded in %d iter\n",maxsw);*/

    return -1;
}

//calculated the center of geometry, co-variance and use the Jacobi method to determine the eigenvalues and
//eigenvectors that defines defined atom groups.
void CalculateEigenValuesVectors(frame *framedata, params *me, std::vector<std::vector<double>> *eigenValVec, 
    std::vector<std::vector<double>> *cog){ 
        const double pi = 3.1415926535;

        //find the eigenvalue and eigenvector for each specified atom group
        //and put both in a nx4 matrix
        std::vector<std::vector<double>> frame_eigenValVec (me->dmov_atomgroup.size(), std::vector<double> (12));
        std::vector<std::vector<double>> frame_cog (me->dmov_atomgroup.size(), std::vector<double> (3));

        for (int h = 0; h < me->dmov_atomgroup.size(); h++)
        {           
            //matrix structure
            //            0  1  2
            // matrix  =  1  3  4
            //            2  4  5
            double matrix[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            double eigenvalues[] = { 0, 0, 0 };
            double eigenvectors[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            //some more variables
            int c = 0;
            double *x = framedata->x.data();
            double sum_x = 0;
            double avg_x = 0;
            double *y = framedata->y.data();
            double sum_y = 0;
            double avg_y = 0;
            double *z = framedata->z.data();
            double sum_z = 0;
            double avg_z = 0;

            //number of atoms in this atom group
            c = me->dmov_atomgroup[h].size();
            int *atoms = me->dmov_atomgroup[h].data();

            //calculated center of geometry for each dimension separately
            for (int i = 0; i < c; i++)
            {
                sum_x += x[atoms[i]];
            }
            avg_x = sum_x / c;

            for (int i = 0; i < c; i++)
            {
                sum_y += y[atoms[i]];
            }
            avg_y = sum_y / c;

            for (int i = 0; i < c; i++)
            {            
                sum_z += z[atoms[i]];
            }
            avg_z = sum_z / c;

            //center of geometry of atom groups in this frame
            std::vector<double> cog_temp;
            cog_temp.push_back(avg_x);
            cog_temp.push_back(avg_y);
            cog_temp.push_back(avg_z);
            frame_cog[h] = cog_temp;

            for (int i = 0; i < c; i++)
            {
                matrix[0] += (avg_x - x[atoms[i]]) * (avg_x - x[atoms[i]]);
            }

            for (int i = 0; i < c; i++)
            {
                matrix[1] += (avg_x - x[atoms[i]]) * (avg_y - y[atoms[i]]);
            }

            for (int i = 0; i < c; i++)
            {
                matrix[2] += (avg_x - x[atoms[i]]) * (avg_z - z[atoms[i]]);
            }

            matrix[3] = matrix[1];

            for (int i = 0; i < c; i++)
            {
                matrix[4] += (avg_y - y[atoms[i]]) * (avg_y - y[atoms[i]]);
            }

            for (int i = 0; i < c; i++)
            {
                matrix[5] += (avg_y - y[atoms[i]]) * (avg_z - z[atoms[i]]);
            }

            matrix[6] = matrix[2];

            matrix[7] = matrix[5];

            for (int i = 0; i < c; i++)
            {
                matrix[8] += (avg_z - z[atoms[i]]) * (avg_z - z[atoms[i]]);
            }

            for (int i = 0; i < 9; i++)
            {
                matrix[i] = matrix[i]/c;
            }

#ifdef DEBUG
            std::cout << "matrix\n" <<
                "---------------------------------\n" << 
                std::fixed << 
                std::setw(10) << std::setprecision(4) << matrix[0] << " " <<
                std::setw(10) << std::setprecision(4) << matrix[1] << " " <<
                std::setw(10) << std::setprecision(4) << matrix[2] << "\n" <<
                std::setw(10) << std::setprecision(4) << matrix[3] << " " <<
                std::setw(10) << std::setprecision(4) << matrix[4] << " " <<
                std::setw(10) << std::setprecision(4) << matrix[5] << "\n" <<
                std::setw(10) << std::setprecision(4) << matrix[6] << " " <<
                std::setw(10) << std::setprecision(4) << matrix[7] << " " <<
                std::setw(10) << std::setprecision(4) << matrix[8] << std::endl;

            std::cout << "\ncoordinates\n" <<
                "---------------------------------\n" << 
                std::fixed << 
                std::setw(10) << std::setprecision(10) << avg_x << " " <<
                std::setw(10) << std::setprecision(10) << sum_x << "\n" <<
                std::setw(10) << std::setprecision(10) << avg_y << " " <<
                std::setw(10) << std::setprecision(10) << sum_y << "\n" <<
                std::setw(10) << std::setprecision(10) << avg_z << " " <<
                std::setw(10) << std::setprecision(10) << sum_z << std::endl;  
#endif // DEBUG

            //jacobi
            Jacobi(matrix, eigenvalues, eigenvectors, me->jacobi_max_iteration);

            //output matrix
            std::vector<std::vector<double>> test (3, std::vector<double> (4));
            test[0][0] = eigenvalues[0];
            test[0][1] = eigenvectors[0];
            test[0][2] = eigenvectors[1];
            test[0][3] = eigenvectors[2];
            test[1][0] = eigenvalues[1];
            test[1][1] = eigenvectors[3];
            test[1][2] = eigenvectors[4];
            test[1][3] = eigenvectors[5];
            test[2][0] = eigenvalues[2];
            test[2][1] = eigenvectors[6];
            test[2][2] = eigenvectors[7];
            test[2][3] = eigenvectors[8];

            //sort the eigenvectors based on eigenvalues, largest eigenvalue first
            std::stable_sort(test.begin(), test.end(), sortFunction);

            frame_eigenValVec[h][0] = test[0][0];
            frame_eigenValVec[h][1] = test[0][1];
            frame_eigenValVec[h][2] = test[0][2];
            frame_eigenValVec[h][3] = test[0][3];
            frame_eigenValVec[h][4] = test[1][0];
            frame_eigenValVec[h][5] = test[1][1];
            frame_eigenValVec[h][6] = test[1][2];
            frame_eigenValVec[h][7] = test[1][3];
            frame_eigenValVec[h][8] = test[2][0];
            frame_eigenValVec[h][9] = test[2][1];
            frame_eigenValVec[h][10] = test[2][2];
            frame_eigenValVec[h][11] = test[2][3];

#ifdef DEBUG
            std::cout << "\neigenvalues\n" <<
                "---------------------------------\n" << 
                std::fixed << 
                std::setw(10) << std::setprecision(4) << frame_eigenValVec[h][0] << " " <<
                std::setw(10) << std::setprecision(4) << frame_eigenValVec[h][1] << " " <<
                std::setw(10) << std::setprecision(4) << frame_eigenValVec[h][2] << "\n" <<
                "\neigenvectors\n" <<
                "---------------------------------\n" << 
                std::setw(10) << std::setprecision(4) << frame_eigenValVec[h][3] << " " <<
                std::setw(10) << std::setprecision(4) << frame_eigenValVec[h][4] << " " <<
                std::setw(10) << std::setprecision(4) << frame_eigenValVec[h][5] << "\n" <<
                std::setw(10) << std::setprecision(4) << frame_eigenValVec[h][6] << " " <<
                std::setw(10) << std::setprecision(4) << frame_eigenValVec[h][7] << " " <<
                std::setw(10) << std::setprecision(4) << frame_eigenValVec[h][8] << "\n" <<
                std::setw(10) << std::setprecision(4) << frame_eigenValVec[h][9] << " " <<
                std::setw(10) << std::setprecision(4) << frame_eigenValVec[h][10] << " " <<
                std::setw(10) << std::setprecision(4) << frame_eigenValVec[h][11] << std::endl;  
#endif // DEBUG

        }
        //http://www.baptiste-wicht.com/2012/03/cp11-concurrency-tutorial-part-2-protect-shared-data/
        *eigenValVec = frame_eigenValVec;
        *cog = frame_cog;
}

//correction for the sign change of vectors in principal component result
void CorrectionSignEigenValuesVectors(params *me, std::vector<std::vector<std::vector<double>>> *eigenValVec){ 
    //define some constants
    const double pi = 3.1415926535;
    const double cut_off = me->jacobi_vector_correction_cutoff;

    //eigenvalvec is a three dimensional std::vector frame<atom_group<eigenvalues + eigenvectors>>
    // f {e1, v_x1, v_y1, v_z1, .....}
    //define a reference when trying to process the first frame of the entire trajectory
    if (me->dmov_calc_previous_last_eigenvalvec.empty())
    { 
        me->dmov_calc_previous_last_eigenvalvec = (*eigenValVec)[0];
    }

    //handle the first frame of a batch with respect to the last frame fo the previous batch
    //for every define atom group
    for (int j = 0; j < (*eigenValVec)[0].size(); j++)
    {
        //for each calculated eigenvector of a requested atom group
        for (int k = 0; k < 12; k+=4)
        {
            //detect a sign change with respect to a preceding frame
            int change = 
                ((*eigenValVec)[0][j][k+1] * me->dmov_calc_previous_last_eigenvalvec[j][k+1] < 0 && std::abs((*eigenValVec)[0][j][k+1] - me->dmov_calc_previous_last_eigenvalvec[j][k+1]) > cut_off) +
                ((*eigenValVec)[0][j][k+2] * me->dmov_calc_previous_last_eigenvalvec[j][k+2] < 0 && std::abs((*eigenValVec)[0][j][k+2] - me->dmov_calc_previous_last_eigenvalvec[j][k+2]) > cut_off) +
                ((*eigenValVec)[0][j][k+3] * me->dmov_calc_previous_last_eigenvalvec[j][k+3] < 0 && std::abs((*eigenValVec)[0][j][k+3] - me->dmov_calc_previous_last_eigenvalvec[j][k+3]) > cut_off);
            //correct for the sign change by *-1
            //the vectors i.e. the main component of each principle component can change direction -180 to +0 or 0 to +180
            if (change >= 1)
            {
                (*eigenValVec)[0][j][k+1] = -1 * (*eigenValVec)[0][j][k+1];
                (*eigenValVec)[0][j][k+2] = -1 * (*eigenValVec)[0][j][k+2];
                (*eigenValVec)[0][j][k+3] = -1 * (*eigenValVec)[0][j][k+3];
            }
        }
    }

    //apply the correction for the first frame of the batch too all frames in the batch
    //for each frame
    for (int i = 1; i < (*eigenValVec).size(); i++)
    {
        //for each defined atom group
        for (int j = 0; j < (*eigenValVec)[i].size(); j++)
        {
            //for each calculated eigenvector of a requested atom group
            for (int k = 0; k < 12; k+=4)
            {
                //detect a sign change with respect to a preceding frame
                int change = 
                    ((*eigenValVec)[i][j][k+1] * (*eigenValVec)[i-1][j][k+1] < 0 && std::abs((*eigenValVec)[i][j][k+1] - (*eigenValVec)[i-1][j][k+1]) > cut_off) +
                    ((*eigenValVec)[i][j][k+2] * (*eigenValVec)[i-1][j][k+2] < 0 && std::abs((*eigenValVec)[i][j][k+2] - (*eigenValVec)[i-1][j][k+2]) > cut_off) +
                    ((*eigenValVec)[i][j][k+3] * (*eigenValVec)[i-1][j][k+3] < 0 && std::abs((*eigenValVec)[i][j][k+3] - (*eigenValVec)[i-1][j][k+3]) > cut_off);
                //correct for the sign change by *-1
                //the vectors i.e. the main component of each principle component can change direction -180 to +0 or 0 to +180
                if (change >= 1)
                {
                    (*eigenValVec)[i][j][k+1] = -1 * (*eigenValVec)[i][j][k+1];
                    (*eigenValVec)[i][j][k+2] = -1 * (*eigenValVec)[i][j][k+2];
                    (*eigenValVec)[i][j][k+3] = -1 * (*eigenValVec)[i][j][k+3];
                }
            }
        }
    }

    //pass the last set of eigenvectors and eigenvalues of the last frame in this batch to the next batch processing
    me->dmov_calc_previous_last_eigenvalvec = (*eigenValVec)[(*eigenValVec).size()-1];
}

//correct for rotation
//fits the protein onto the z-axis.
void CorrectRotationSolute(frame *framedata, params *me, std::vector<std::vector<double>> *eigenValVec) {
    std::vector<double> soluteEigenValVec = (*eigenValVec)[0];
    double referenceAxis[] = { 0, 0, 1 };
    double cp[] = { 0, 0, 0 };
    double sign = 0;
    double a_r = 0;

    //cross-product
    cp[0] = soluteEigenValVec[2] * referenceAxis[2] - soluteEigenValVec[3] * referenceAxis[1];
    cp[1] = soluteEigenValVec[3] * referenceAxis[0] - soluteEigenValVec[1] * referenceAxis[2];
    cp[2] = soluteEigenValVec[1] * referenceAxis[1] - soluteEigenValVec[2] * referenceAxis[0];

    //calculated the angle
    a_r = ((sign > 0) - (sign < 0)) * 
        acos(
        (cp[0] * referenceAxis[0] + cp[1] * referenceAxis[1] + cp[2] * referenceAxis[2]) / (
        sqrt(cp[0] * cp[0] + cp[1] * cp[1] + cp[2] * cp[2]) *
        sqrt(referenceAxis[0] * referenceAxis[0] + referenceAxis[1] * referenceAxis[1] + referenceAxis[2] * referenceAxis[2])
        ));

    //Euclidian-space formulation from 
    //http://www.gamedev.net/page/resources/_/technical/math-and-physics/do-we-really-need-quaternions-r1199
    double c = cos(a_r);
    double s = sin(a_r);
    double t = 1-c;

    //calculate rotation matrix
    double RotationMatrix[3][3] = {}; 
    RotationMatrix[0][0] = t*cp[0]*cp[0] + c;
    RotationMatrix[0][1] = t*cp[0]*cp[1] + s*cp[2];
    RotationMatrix[0][2] = t*cp[0]*cp[2] - s*cp[1];
    RotationMatrix[1][0] = t*cp[0]*cp[1] - s *cp[2];
    RotationMatrix[1][1] = t*cp[1]*cp[1] + c;
    RotationMatrix[1][2] = t*cp[1]*cp[2] + s*cp[0];
    RotationMatrix[2][0] = t*cp[0]*cp[2] + s*cp[1];
    RotationMatrix[2][1] = t*cp[1]*cp[2] + s*cp[0];
    RotationMatrix[2][2] = t*cp[2]*cp[2] + c;

    //transform all x coordinates    
    double *x = framedata->x.data();
    double *y = framedata->y.data();
    double *z = framedata->z.data();

    std::vector<std::vector<double>> dim_temp (3, std::vector<double> (framedata->x.size()));
    for (int i = 0; i < framedata->x.size(); i++)
    {
        dim_temp[0][i] = RotationMatrix[0][0] * x[i] + 
            RotationMatrix[0][1] * y[i] + 
            RotationMatrix[0][2] * z[i];
    }

    //transform all y coordinates
    for (int i = 0; i < framedata->y.size(); i++)
    {
        dim_temp[1][i] = RotationMatrix[1][0] * x[i] + 
            RotationMatrix[1][1] * y[i] + 
            RotationMatrix[1][2] * z[i];
    }

    //transform all z coordinates
    for (int i = 0; i < framedata->z.size(); i++)
    {
        dim_temp[2][i] = RotationMatrix[2][0] * x[i] + 
            RotationMatrix[2][1] * y[i] + 
            RotationMatrix[2][2] * z[i];
    }

    framedata->x = dim_temp[0];
    framedata->y = dim_temp[1];
    framedata->z = dim_temp[2];
}

//code checked 20130122: CORRECT!
//write out the data in either CNF or PDB compatible format
//this is not for production... code has to be written differently 
void WriteOutFrame(frame *framedata, gz::ogzstream &outfile, params *me) {
    //check if it is possible to read file
    if (!outfile)
    {
        std::cerr << "cannot open otuput file" << "\n";
    }
    else {


        //quickfix
        std::vector<int> temp_fix; 
        for (int i = 0; i < me->solute_count; i++)
        {
            temp_fix.push_back(me->solute_molecules[i][0]-1);
            temp_fix.push_back(me->solute_molecules[i][1]-1);
        }
        for (int i = 0; i < me->solutes_cog_molecules.size(); i++)
        {
            temp_fix.push_back(me->solutes_cog_molecules[i][0]-1);
            temp_fix.push_back(me->solutes_cog_molecules[i][1]-1);
        }
        //quickfix

        if (me->outformat=="cnf" || me->outformat=="trc")
        {
            outfile << "TIMESTEP" << "\n";
            outfile << " " << std::setw(17) << std::setprecision(0) << framedata->timestep << " " << std::setw(19) << std::fixed << std::setprecision(9) << framedata->time << "\n";
            outfile << "END" << "\n";
            outfile << "POSITION" << "\n";	
            outfile << std::fixed << std::setprecision(9);

            //quickfix
            for (int i = 0; i < me->atomrecords; i++)
            {
                if (me->solvent_skip)
                {
                    for (int ii = 0; ii < temp_fix.size(); ii+=2)
                    {
                        if (i>=temp_fix[ii] && i<=temp_fix[ii+1])
                        {
                            outfile << framedata->prefix[i] << " " << std::setw(14) << framedata->x[i] << " " << std::setw(14) << framedata->y[i] << " " << std::setw(14) << framedata->z[i] << "\n";
                        }
                    }
                }
                else
                {                    
                    outfile << framedata->prefix[i] << " " << std::setw(14) << framedata->x[i] << " " << std::setw(14) << framedata->y[i] << " " << std::setw(14) << framedata->z[i] << "\n";
                }
            }
            //quickfix

            outfile << "END" << "\n";
            outfile << "GENBOX" << "\n";
            outfile << " " << framedata->boxtype << "\n";
            outfile << std::fixed << std::setprecision(9);
            outfile << " " << std::setw(14) << framedata->box_length_x << " " << std::setw(14) << framedata->box_length_y << " " << std::setw(14) << framedata->box_length_z << "\n";
            outfile << " " << std::setw(14) << framedata->box_angle_x << " " << std::setw(14) << framedata->box_angle_y << " " << std::setw(14) << framedata->box_angle_z << "\n";
            outfile << " " << std::setw(14) << framedata->box_3_x << " " << std::setw(14) << framedata->box_3_y << " " << std::setw(14) << framedata->box_3_z << "\n";
            outfile << " " << std::setw(14) << framedata->box_4_x << " " << std::setw(14) << framedata->box_4_y << " " << std::setw(14) << framedata->box_4_z << "\n";
            outfile << "END" << "\n";
        }
        if (me->outformat=="pdb")
        {
            //COLUMNS       DATA TYPE     FIELD         DEFINITION
            //--------------------------------------------------------------------------------------
            // 1 -  6       Record name   "REMARK"
            // 8 - 10       Integer       remarkNum     Remark  number. It is not an error for
            //                                          remark n to exist in an entry when
            //                                          remark n-1 does not.
            //12 - 79       LString       empty         Left  as white space in first line
            //                                          of each  new remark.
            outfile << "REMARK    1 " << std::setw(17) << std::setprecision(0) << framedata->timestep << " " << std::setw(19) << std::fixed << std::setprecision(9) << framedata->time << "\n";  

            //outfile << "MODEL " << std::setw(8) << "\n";
            outfile << "MODEL \n";
            //COLUMNS       DATA  TYPE    FIELD          DEFINITION
            //-------------------------------------------------------------
            // 1 -  6       Record name   "CRYST1"
            // 7 - 15       Real(9.3)     a              a (Angstroms).
            //16 - 24       Real(9.3)     b              b (Angstroms).
            //25 - 33       Real(9.3)     c              c (Angstroms).
            //34 - 40       Real(7.2)     alpha          alpha (degrees).
            //41 - 47       Real(7.2)     beta           beta (degrees).
            //48 - 54       Real(7.2)     gamma          gamma (degrees).
            //56 - 66       LString       sGroup         Space  group.
            //67 - 70       Integer       z              Z value.
            //CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1
            outfile << "CRYST " 
                << std::setw(9) << std::setprecision(3) << framedata->box_length_x
                << std::setw(9) << std::setprecision(3) << framedata->box_length_y
                << std::setw(9) << std::setprecision(3) << framedata->box_length_z
                << std::setw(7) << std::setprecision(2) << framedata->box_angle_x
                << std::setw(7) << std::setprecision(2) << framedata->box_angle_y
                << std::setw(7) << std::setprecision(2) << framedata->box_angle_z
                << " P 1           1"
                << "\n";
            outfile << std::fixed << std::setprecision(9);
            for (int i = 0; i < me->atomrecords; i++)
            {
                //quickfix
                if (me->solvent_skip)
                {
                    for (int ii = 0; ii < temp_fix.size(); ii+=2)
                    {
                        if (i>=temp_fix[ii] && i<=temp_fix[ii+1])
                        {
                            //quickfix
                            //    1 ASN   H1         1    1.021435895    2.079909498    0.623854235
                            //ATOM      1  N   THR A   1      -0.313  18.726  33.523  1.00 21.00           N
                            // 1 -  6        Record name   "ATOM  "
                            // 7 - 11        Integer       serial       Atom  serial number.
                            //13 - 16        Atom          name         Atom name.
                            //17             Character     altLoc       Alternate location indicator.
                            //18 - 20        Residue name  resName      Residue name.
                            //22             Character     chainID      Chain identifier.
                            //23 - 26        Integer       resSeq       Residue sequence number.
                            //27             AChar         iCode        Code for insertion of residues.
                            //31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
                            //39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
                            //47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
                            //55 - 60        Real(6.2)     occupancy    Occupancy.
                            //61 - 66        Real(6.2)     tempFactor   Temperature  factor.
                            //77 - 78        LString(2)    element      Element symbol, right-justified.
                            //79 - 80        LString(2)    charge       Charge  on the atom.
                            outfile << "ATOM  "
                                << std::setw(5) << framedata->prefix[i].substr(19, 5) << " "
                                << " " << std::setw(3) << framedata->prefix[i].substr(12, 3)
                                << " "
                                << std::setw(3) << framedata->prefix[i].substr(6, 3)
                                << " A"
                                << std::setw(4) << framedata->prefix[i].substr(1, 4)
                                << "    "
                                << std::fixed
                                << std::setw(8) << std::setprecision(3) << framedata->x[i] * 10
                                << std::setw(8) << std::setprecision(3) << framedata->y[i] * 10
                                << std::setw(8) << std::setprecision(3) << framedata->z[i] * 10
                                << std::setw(6) << std::setprecision(2) << 1.0
                                << std::setw(6) << std::setprecision(2) << 1.0 << "\n"                
                                ;
                        }
                    }
                }
                else
                {
                    //quickfix
                    //    1 ASN   H1         1    1.021435895    2.079909498    0.623854235
                    //ATOM      1  N   THR A   1      -0.313  18.726  33.523  1.00 21.00           N
                    // 1 -  6        Record name   "ATOM  "
                    // 7 - 11        Integer       serial       Atom  serial number.
                    //13 - 16        Atom          name         Atom name.
                    //17             Character     altLoc       Alternate location indicator.
                    //18 - 20        Residue name  resName      Residue name.
                    //22             Character     chainID      Chain identifier.
                    //23 - 26        Integer       resSeq       Residue sequence number.
                    //27             AChar         iCode        Code for insertion of residues.
                    //31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
                    //39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
                    //47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
                    //55 - 60        Real(6.2)     occupancy    Occupancy.
                    //61 - 66        Real(6.2)     tempFactor   Temperature  factor.
                    //77 - 78        LString(2)    element      Element symbol, right-justified.
                    //79 - 80        LString(2)    charge       Charge  on the atom.
                    outfile << "ATOM  "
                        << std::setw(5) << framedata->prefix[i].substr(19, 5) << " "
                        << " " << std::setw(3) << framedata->prefix[i].substr(12, 3)
                        << " "
                        << std::setw(3) << framedata->prefix[i].substr(6, 3)
                        << " A"
                        << std::setw(4) << framedata->prefix[i].substr(1, 4)
                        << "    "
                        << std::fixed
                        << std::setw(8) << std::setprecision(3) << framedata->x[i] * 10
                        << std::setw(8) << std::setprecision(3) << framedata->y[i] * 10
                        << std::setw(8) << std::setprecision(3) << framedata->z[i] * 10
                        << std::setw(6) << std::setprecision(2) << 1.0
                        << std::setw(6) << std::setprecision(2) << 1.0 << "\n"                
                        ;
                }
            }
            if (me->cog_write)
            {
                outfile << "HETATM" << std::setw(5) << me->atomrecords+2 << "  ZN   ZN A9999    "
                    << std::fixed
                    << std::setw(8) << std::setprecision(3) << framedata->solute_cog_x * 10
                    << std::setw(8) << std::setprecision(3) << framedata->solute_cog_y * 10
                    << std::setw(8) << std::setprecision(3) << framedata->solute_cog_z * 10
                    << "  1.00  0.00          ZN  \n";
            }
            outfile << "ENDMDL" << "\n";
        }
    }
}

//code checked 20130122: CORRECT!
char* XPathGetText(std::string xpath_query, xmlXPathContextPtr xpathCtx) 
{
    xmlXPathObjectPtr xpathObj = xmlXPathEvalExpression(BAD_CAST xpath_query.c_str(), xpathCtx);
    return (char*)xmlNodeGetContent(xpathObj->nodesetval->nodeTab[0]);	
}

//code checked 20130122: CORRECT!
//reads a references file for topology information that is needed for writting out
//files in CNF and PDB compatible format
void TrcReferenceFrame(std::vector<std::string> *prefix, std::string trc_reference) {

    std::ifstream infile(trc_reference);

    //check if it is possible to read file
    if (!infile)
    {
        std::cerr << "cannot open otuput file" << "\n";
    }
    else {
        //read line by line as string while not end-of-file
        int positionBlock_counter = 0;
        bool isPositionBlock = false;
        while (!infile.eof()) {
            std::string line;
            std::getline(infile, line);

            //ignore comments and empty lines
            if (line[0] != '#' && line.length() > 0)
            {
                if (line.substr(0,8) == "POSITION") 
                {
                    isPositionBlock = true;
                }
                else if (line.substr(0,3) == "END") 
                {
                    isPositionBlock = false;
                }
                else if (isPositionBlock)
                {                    
                    (*prefix)[positionBlock_counter] = line.substr(0,24);
                    positionBlock_counter += 1;
                }
            }
        }
    }
    // std::cout << "using reference file"  << std::endl;
}

//code checked 20130122: CORRECT!
//parse parameter input file
void ParseParamsInput(params *this_params, std::string job_id, std::string param_file) {
    //initialize libxml
    std::string xpath;

    xmlDocPtr doc;
    xmlXPathContextPtr xpathCtx; 
    xmlXPathObjectPtr xpathObj; 

    //load XML document
    doc = xmlParseFile(param_file.c_str());
    if (doc == NULL)
    {
        std::cout << "Unable to parse " << param_file << std::endl;
        exit(-1);
    }

    //create xpath context
    xpathCtx = xmlXPathNewContext(doc);
    if(xpathCtx == NULL) {
        std::cout << "Unable to create new XPath context " << param_file << std::endl;
        xmlFreeDoc(doc);
        exit(-1);
    }

    //evaluate xpath expression
    xpath = "/me/job[@id=\"" + job_id + "\"]";
    xpathObj = xmlXPathEvalExpression(BAD_CAST xpath.c_str(), xpathCtx);
    if(xpathObj == NULL) {
        std::cout << "Error: unable to evaluate xpath expression\n  " << xpath << std::endl;
        xmlXPathFreeContext(xpathCtx); 
        xmlFreeDoc(doc); 
        exit(-1);
    }

    for (int i = 0; i < xpathObj->nodesetval->nodeNr; i++)
    {
        //read current node;
        std::cout << xpathObj->nodesetval->nodeTab[i]->name << std::endl;
        std::cout << xmlGetProp(xpathObj->nodesetval->nodeTab[i], BAD_CAST ((std::string)"id").c_str())  << std::endl;

        //read the topology section
        xpathCtx->node = xpathObj->nodesetval->nodeTab[i];
        //get number of atoms records in the topology section
        this_params->atomrecords = atoi(XPathGetText("./topology/atom_records_count", xpathCtx));
        if (this_params->atomrecords < 0)
        {
            this_params->atomrecords = 0;
        }
        //get solvent definitions
        //get first atom
        this_params->solvent_atoms_first = atoi(XPathGetText("./topology/solvent/@first_atom", xpathCtx));
        //get last atom
        this_params->solvent_atoms_last = atoi(XPathGetText("./topology/solvent/@last_atom", xpathCtx));
        //skip or keep solvent
        this_params->solvent_skip = atoi(XPathGetText("./topology/solvent/@skip", xpathCtx));
        //get number of atoms
        this_params->solvent_size = atoi(XPathGetText("./topology/solvent/number_of_atoms", xpathCtx));
        //get dimension expansion factor
        this_params->solvent_dimension_expansion = atoi(XPathGetText("./topology/solvent/dimension_expansion", xpathCtx));

#ifdef DEBUG
        std::cout << "got solvent parameters" << std::endl;
#endif // DEBUG

        //get solute (cog) definitions
        xpath = "./topology/solutes_cog/solute";
        xmlXPathObjectPtr COGSolute = xmlXPathEvalExpression(BAD_CAST xpath.c_str(), xpathCtx);
        //loops through solvents
        //std::cout << xmlNodeGetContent(xpathObj->nodesetval->nodeTab[0]) << std::endl;
        //get all COG solutes
        for (int j = 0; j < COGSolute->nodesetval->nodeNr; j++)
        {
            std::vector<int> temp;
            xpathCtx->node = COGSolute->nodesetval->nodeTab[j];
            //get first atom of COG solute
            temp.push_back(atoi(XPathGetText("./@first_atom", xpathCtx)));
            //get last atom of COG solute
            temp.push_back(atoi(XPathGetText("./@last_atom", xpathCtx)));
            //get number of atom of COG solute
            temp.push_back(atoi(XPathGetText("./number_of_atoms", xpathCtx)));	
            //get dimension expansion factor
            temp.push_back(atoi(XPathGetText("./dimension_expansion", xpathCtx)));
            //skip or keep COG solute
            temp.push_back(atoi(XPathGetText("./@skip", xpathCtx)));	
            this_params->solutes_cog_molecules.push_back(temp);
            temp.clear();

        }
        xmlXPathFreeObject(COGSolute);

#ifdef DEBUG
        std::cout << "got cog_solutes parameters" << std::endl;
#endif // DEBUG

        xpathCtx->node = xpathObj->nodesetval->nodeTab[i];
        //get solute definitions
        xpath = "./topology/solutes/solute";
        xmlXPathObjectPtr solutes = xmlXPathEvalExpression(BAD_CAST xpath.c_str(), xpathCtx);
        //get all solutes
        for (int j = 0; j < solutes->nodesetval->nodeNr; j++)
        {
            std::vector<int> temp;
            xpathCtx->node = solutes->nodesetval->nodeTab[j];
            //get first atom of solute
            temp.push_back(atoi(XPathGetText("./@first_atom", xpathCtx)));
            //get last atom of solute
            temp.push_back(atoi(XPathGetText("./@last_atom", xpathCtx)));
            //get dimension expansion factor
            temp.push_back(atoi(XPathGetText("./dimension_expansion", xpathCtx)));
            //skip or keep solute
            temp.push_back(atoi(XPathGetText("./@skip", xpathCtx)));
            this_params->solute_molecules.push_back(temp);
            temp.clear();
        };
        this_params->solute_count = this_params->solute_molecules.size();
        xmlXPathFreeObject(solutes);

        //read the variables section
        xpathCtx->node = xpathObj->nodesetval->nodeTab[i];

        //first get the frameoutX specific options
        //write out original cog or not
        this_params->cog_write = atoi(XPathGetText("./analysis/frameout/cog_write", xpathCtx));
        //correct by shifting cog or not
        this_params->cog_correction = atoi(XPathGetText("./analysis/frameout/cog_correction", xpathCtx));
        //get distance cut
        this_params->distance_cut_off = atof(XPathGetText("./analysis/frameout/distance_cut_off", xpathCtx)) * atof(XPathGetText("./analysis/frameout/distance_cut_off", xpathCtx));
        //do rotation correction or not
        this_params->rotation_correction = atoi(XPathGetText("./analysis/frameout/rotation_correction", xpathCtx));
        //Jacobi maximum number of iterations
        this_params->jacobi_max_iteration = atoi(XPathGetText("./analysis/frameout/jacobi/@max_iteration", xpathCtx));
        //Jacobi eigenvector direction correction cut-off value
        this_params->jacobi_vector_correction_cutoff = atof(XPathGetText("./analysis/frameout/jacobi/@eigenvector_correction_cut_off", xpathCtx));

        //get number of frames per thread
        this_params->num_frame_per_thread = atoi(XPathGetText("./variables/frames_per_thread", xpathCtx));
        //get number of thread
        this_params->num_thread = atoi(XPathGetText("./variables/number_of_threads", xpathCtx));
        if (this_params->num_thread <= 0)
        {
            this_params->num_thread = std::thread::hardware_concurrency();
        }
        //get multiplier for thread
        if (atoi(XPathGetText("./variables/number_of_threads_multiplier", xpathCtx)) <= 0)
        {            
            this_params->num_thread_real = this_params->num_thread * 1;
        }
        else 
        {        
            this_params->num_thread_real = this_params->num_thread * atoi(XPathGetText("./variables/number_of_threads_multiplier", xpathCtx));
        }
        //get dimension expansion factor
        this_params->trc_reference = XPathGetText("./variables/reference_file", xpathCtx);

        //read input section
        //get input data format
        this_params->informat = XPathGetText("./input/format", xpathCtx);

        xpathCtx->node = xpathObj->nodesetval->nodeTab[i];
        //get input file paths
        xpath = "./input/files/file";
        xmlXPathObjectPtr inputFiles = xmlXPathEvalExpression(BAD_CAST xpath.c_str(), xpathCtx);
        for (int j = 0; j < inputFiles->nodesetval->nodeNr; j++)
        {
            xpathCtx->node = inputFiles->nodesetval->nodeTab[j];
            //get a file path
            this_params->input_files.push_back(XPathGetText(".", xpathCtx));
        };
        xmlXPathFreeObject(inputFiles);

        //read output block
        xpathCtx->node = xpathObj->nodesetval->nodeTab[i];
        this_params->outfilename = XPathGetText("./output/filename_prefix", xpathCtx);
        this_params->outformat = XPathGetText("./output/format", xpathCtx);
        this_params->output_fragment_size = atoi(XPathGetText("./output/frames_per_file", xpathCtx));
        this_params->output_fragment_skipframes = atoi(XPathGetText("./output/frame_interval", xpathCtx));
        this_params->output_fragment_skiptime = atof(XPathGetText("./output/time_interval", xpathCtx));
    }

    //clean up
    xmlXPathFreeObject(xpathObj);
    xmlXPathFreeContext(xpathCtx); 
    xmlFreeDoc(doc);
}

//code checked 20130122: CORRECT!
void GenboxParser(frame *currentFrame, int genBox_counter, std::string line) 
{
    switch (genBox_counter)
    {
    case 0:
        currentFrame->boxtype = std::stoi(line);
        break;
    case 1:
        currentFrame->box_length_x = std::stod(line.substr(0,15));
        currentFrame->box_length_y = std::stod(line.substr(15,15));
        currentFrame->box_length_z = std::stod(line.substr(30,15));
        break;
    case 2:
        currentFrame->box_angle_x = std::stod(line.substr(0,15));
        currentFrame->box_angle_y = std::stod(line.substr(15,15));
        currentFrame->box_angle_z = std::stod(line.substr(30,15));
        break;
    case 3:
        currentFrame->box_3_x = std::stod(line.substr(0,15));
        currentFrame->box_3_y = std::stod(line.substr(15,15));
        currentFrame->box_3_z = std::stod(line.substr(30,15));
        break;
    case 4:
        currentFrame->box_4_x = std::stod(line.substr(0,15));
        currentFrame->box_4_y = std::stod(line.substr(15,15));
        currentFrame->box_4_z = std::stod(line.substr(30,15));
        break;
    }
}

//code checked 20130122: CORRECT!
void PositionBlockParser(params &me, std::string &line, int positionBlock_counter, std::vector<std::string> *prefix, std::vector<double> *x, std::vector<double> *y, std::vector<double> *z) 
{
    if (me.informat == "trc")
    {
        (*x)[positionBlock_counter] = std::stod(line.substr(0,15));
        (*y)[positionBlock_counter] = std::stod(line.substr(15,15));
        (*z)[positionBlock_counter] = std::stod(line.substr(30,15));
    }
    if (me.informat == "cnf")
    {
        (*prefix)[positionBlock_counter] = line.substr(0,24);
        (*x)[positionBlock_counter] = std::stod(line.substr(25,15));
        (*y)[positionBlock_counter] = std::stod(line.substr(40,15));
        (*z)[positionBlock_counter] = std::stod(line.substr(55,15));
    }
}

//main program
int main(int argc, char* argv[])
{
    std::string job_id = argv[1];
    std::string param_file = argv[2];

    params me;
    ParseParamsInput(&me, job_id, param_file);

    //decide based on input, how many lines to write out in final output file
    int writeAtomRecordsCount = 0;
    if (me.solvent_skip)
    {
        int solutes_cog_count = me.solutes_cog_molecules.size();
        if (solutes_cog_count > 0)
        {                            
            writeAtomRecordsCount= me.solutes_cog_molecules[solutes_cog_count-1][1];
        }
        else
        {
            writeAtomRecordsCount = me.solute_molecules[me.solute_count-1][1];
        }
    }
    else {
        writeAtomRecordsCount = me.atomrecords;
    }

    //rotation parameters
    std::vector<std::vector<int>> ag (1, std::vector<int> (2));
    ag[0][0] = 0;
    ag[0][1] = me.solute_molecules[me.solute_molecules.size()-1][1];
    me.dmov_atomgroup = ag;

#ifdef DEBUG
    std::cout << "input parameters defined" << std::endl;
    std::cout << std::left;
    std::cout << "  number of atom records          : " << me.atomrecords << std::endl;
    std::cout << "  distance cut-off                : " << me.distance_cut_off << std::endl;
    std::cout << "  input format                    : " << me.informat << std::endl;
    std::cout << "  input files                     : " << std::endl;
    for (unsigned int i = 0; i < me.input_files.size(); i++)
    {
        std::cout << "                                    " << me.input_files[i] << std::endl;
    }
    std::cout << "  write-out COG                   : " << me.cog_write << std::endl;
    std::cout << "  correct COG                     : " << me.cog_correction << std::endl;
    std::cout << "  frames per thread               : " << me.num_frame_per_thread << std::endl;
    std::cout << "  number of threads               : " << me.num_thread << std::endl;
    std::cout << "  real number of threads          : " << me.num_thread_real << std::endl;
    std::cout << "  output filename prefix          : " << me.outfilename << std::endl;
    std::cout << "  output format                   : " << me.outformat << std::endl;
    std::cout << "  output frames per file          : " << me.output_fragment_size << std::endl;
    std::cout << "  output every n frame(s)         : " << me.output_fragment_skipframes << std::endl;
    std::cout << "  output every n picoseconds      : " << me.output_fragment_skiptime << std::endl;
    std::cout << "  reference coordinates           : " << me.ref_coords[0] << " " << me.ref_coords[1] << " " << me.ref_coords[2] << std::endl;
    std::cout << "  skip solvent                    : " << me.solvent_skip << std::endl;
    std::cout << "  solutes cog                     : " << std::endl;
    for (unsigned int i = 0; i < me.solutes_cog_molecules.size(); i++)
    {
        std::cout << "                                    " << me.solutes_cog_molecules[i][0] << " " << me.solutes_cog_molecules[i][1] << " " << me.solutes_cog_molecules[i][2] << " " << me.solutes_cog_molecules[i][3] <<  " " << me.solutes_cog_molecules[i][4] << std::endl;
    }
    std::cout << "  number of solutes               : " << me.solute_count << std::endl;
    std::cout << "  solutes atom numbering          : " << std::endl;
    for (unsigned int i = 0; i < me.solute_molecules.size(); i++)
    {
        std::cout << "                                    " << me.solute_molecules[i][0] << " " << me.solute_molecules[i][1] << " " << me.solute_molecules[i][2] << " " << me.solute_molecules[i][3] << std::endl;
    }
    std::cout << "  solvent cog dimension expansion : " << me.solvent_dimension_expansion << std::endl;
    std::cout << "  solvent size                    : " << me.solvent_size << std::endl;
    std::cout << "  trc refernce file               : " << me.trc_reference << std::endl;
#endif // DEBUG

    //holds the active frames
    int activeFrame_count = me.num_thread_real * me.num_frame_per_thread;
    std::vector<frame> activeFrames (activeFrame_count);
    std::vector<frame> activeFramesCopy (activeFrame_count);

    //eigenvalues and eigenvectors
    std::vector<std::vector<std::vector<double>>> eigenValVec (activeFrame_count);
    std::vector<std::vector<std::vector<double>>> centerOfGeometry (activeFrame_count);

    //performance log
    auto start = std::chrono::system_clock::now();

    //only use the first TITLE block found
    bool firstPass = true;

    //frame counter which can be used for skipping frames
    int frame_counter = 0;
    double frame_time = 0;
    bool processThisFrame = false;

    //position block array size
    std::vector<std::string> prefix (me.atomrecords);
    std::vector<double> x (me.atomrecords);
    std::vector<double> y (me.atomrecords);
    std::vector<double> z (me.atomrecords);

#ifdef DEBUG
    std::cout << "--->atom records c++ vectors initialized" << std::endl;
#endif // DEBUG

    //parse reference file for prefix if it is a trajectory file in TRC format
    if (me.informat=="trc")
    {
        TrcReferenceFrame(&prefix, me.trc_reference);
    }

#ifdef DEBUG
    std::cout << "--->reference frame collected" << std::endl;
#endif // DEBUG

    //determine output filename
    gz::ogzstream outfile;
    if (me.outformat == "pdb")
    {
        outfile.open((me.outfilename + me.outformat + ".gz").c_str(), std::ios_base::out);
    }
    if (me.outformat == "cnf")
    {
        outfile.open((me.outfilename + me.outformat + ".gz").c_str(), std::ios_base::out);
    }

    //output thread
    std::thread outfileThread = std::thread([](){return 0;});

#ifdef DEBUG
    std::cout << "--->output file thread initialized" << std::endl;
#endif // DEBUG

    //process the trajectory files sequentially
    std::cout << "processing trajectory files: " << std::endl;

    for (unsigned int i = 0; i < me.input_files.size(); i++)
    {
        std::cout << "  " << me.input_files[i] << std::endl;

        //remember if box has been shift to gather first atom of a frame
        int init_shift[3] = { 0, 0, 0 };

        //define file to read
        gz::igzstream file(me.input_files[i].c_str());

        //boolean specifying current active block
        bool isTitleBlock = false;
        bool isTimestepBlock = false;
        bool isPositionBlock = false;
        bool isGenboxBlock = false;

        //content holder for the blocks
        std::string titleBlock, timestepBlock, positionBlock, genboxBlock;

        //counters
        int positionBlock_counter = 0;
        int activeFrame_counter = 0;
        int genBox_counter = 0;

        //check if it is possible to read file
        if (!file)
        {
            std::cerr << "cannot open otuput file" << "\n";
        }
        else {

#ifdef DEBUG
            std::cout << "--->starting to read file" << std::endl;
#endif // DEBUG

            //define holder for current frame
            frame currentFrame;

            //read line by line as string while not end-of-file
            while (!file.eof()) {
                std::string line;
                std::getline(file, line);

                //ignore comments and empty lines
                if (line[0] != '#' && line.length() > 0)
                { //detect the start of a block and remember till the end of the block is found
                    if (line.substr(0,6) == "TITLE")
                    {
                        isTitleBlock = true;

#ifdef DEBUG
                        std::cout << "--->    found TITLE" << std::endl;
#endif // DEBUG

                    }
                    else if (line.substr(0,8) == "TIMESTEP")
                    {
                        isTimestepBlock = true;
#ifdef DEBUG
                        std::cout << "--->    found TIMESTEP" << std::endl;
#endif // DEBUG

                    }
                    else if (line.substr(0,8) == "POSITION")
                    {
                        positionBlock_counter = 0;
                        isPositionBlock = true;

#ifdef DEBUG
                        std::cout << "--->    found POSITION" << std::endl;
#endif // DEBUG

                    }
                    else if (line.substr(0,6) == "GENBOX")
                    {
                        genBox_counter = 0;
                        isGenboxBlock = true;

#ifdef DEBUG
                        std::cout << "--->    found GENBOX" << std::endl;
#endif // DEBUG

                    }
                    //reset if end of the block is found
                    else if (line.substr(0,3) == "END")
                    {    
                        //what to do if end of a TITLE block
                        if (isTitleBlock)
                        {
                            //TITLE block
                            if (me.outformat == "cnf" && firstPass)
                            {
                                outfile << "TITLE" << "\n";
                                outfile << titleBlock;
                                outfile << "END" << "\n";
                            }
                        }
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
                            if (frame_counter == 0)
                            {
                                processThisFrame = true;                                
                                if (me.informat == "trc")
                                {
                                    frame_time = std::stod(timestepBlock.substr(15,15));
                                }
                                else if (me.informat == "cnf")
                                {
                                    frame_time = std::stod(timestepBlock.substr(19,19));                                    
                                }
                            }
                            //if skip by frame counting
                            else if (me.output_fragment_skipframes > 0 && (frame_counter % me.output_fragment_skipframes) == 0)
                            {
                                processThisFrame = true;
                            } //if skip by time interval
                            else if (me.output_fragment_skiptime > 0)
                            {   
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
                            if (me.informat == "trc")
                            {
                                currentFrame.time = std::stod(timestepBlock.substr(15,15));
                                currentFrame.timestep = std::stol(timestepBlock.substr(0,15));
                            }
                            else if (me.informat == "cnf")
                            {
                                currentFrame.time = std::stod(timestepBlock.substr(19,38));
                                currentFrame.timestep = std::stol(timestepBlock.substr(0,19));
                            }

#ifdef DEBUG
                            if (processThisFrame)
                            {
                                std::cout << "      got the time and step ( " << currentFrame.time << " / " << frame_counter << " )" << std::endl;
                            }
#endif // DEBUG

                        }
                        if (isPositionBlock)
                        {
                            currentFrame.prefix = prefix;
                            currentFrame.x = x;
                            currentFrame.y = y;
                            currentFrame.z = z;

#ifdef DEBUG
                            if (processThisFrame)
                            {
                                std::cout << "      got the coordinates " << std::endl;
                            }
#endif // DEBUG

                        } //what to do if end of a GENBOX block
                        if (isGenboxBlock)
                        {
                            //the first frame has been collected, it's no longer the first pass through the loops
                            firstPass = false;

                            int first_atom_solute = 0;
                            if (me.solute_count >  0)
                            {
                                first_atom_solute = me.solute_molecules[0][0] - 1;
                            }

                            //skip frame[0] from the reference process, because it is the first frame 
                            //and has no reference coordinates for the first atom
                            //also get the time of the first frame as reference for time-based skipping
                            if (frame_counter==0)
                            {
                                currentFrame.init_shift_x = 0;
                                currentFrame.init_shift_y = 0;
                                currentFrame.init_shift_z = 0;
                                me.ref_coords[0] = x[first_atom_solute];
                                me.ref_coords[1] = y[first_atom_solute];
                                me.ref_coords[2] = z[first_atom_solute];
                            }
                            else 
                            {
                                //do not combine the following lines
                                //gather should be kept separate from reference coordinates assignment
                                FirstAtomBasedBoxShifter(&currentFrame, first_atom_solute, &me);
                                me.ref_coords[0] = currentFrame.x[first_atom_solute];
                                me.ref_coords[1] = currentFrame.y[first_atom_solute];
                                me.ref_coords[2] = currentFrame.z[first_atom_solute];
                            }

                            //until n frames have been read in
                            if (processThisFrame)
                            {
                                activeFrames[activeFrame_counter] = currentFrame;
                                activeFrame_counter += 1;
                            }

                            //if after processing the GENBOX block, n frames are stored than do something
                            if (activeFrame_counter == activeFrame_count)
                            {
#ifdef DEBUG
                                for (int x = 0; x < activeFrames.size(); x++)
                                {
                                    std::cout << activeFrames[x].time << " " << activeFrames[x].timestep << std::endl;  
                                }

                                std::cout << "--->processing the batch ( " << activeFrame_counter << " / " << activeFrame_count << " )" << std::endl;
#endif // DEBUG

                                //STEP 2
                                //C++11 threading (batched)
                                std::thread *myThreads = new std::thread[me.num_thread_real];

#ifdef DEBUG
                                std::cout << "      threads c++ vector initialized" << std::endl;
#endif // DEBUG

                                try
                                {
                                    for (int g = 0; g < me.num_thread_real; g++)
                                    {
                                        myThreads[g] = std::thread([&activeFrames, &me, g, &eigenValVec, &centerOfGeometry]()
                                        {
                                            for (int f = (g * me.num_frame_per_thread); f < ((g + 1) * me.num_frame_per_thread); f++)
                                            {
                                                //STEP 2a: //solute gathering
                                                SoluteGatherer(&activeFrames[f], f, &me);

#ifdef DEBUG
                                                std::cout << "      threads " << g << ": solutes gathered" << std::endl;
#endif // DEBUG

                                                for (unsigned int j = 0; j < me.solutes_cog_molecules.size(); j++)
                                                {
                                                    //STEP 2b: gather the non-solvent with respect to the COG of solutes
                                                    COGGatherer(&activeFrames[f], me.solutes_cog_molecules[j][0],  me.solutes_cog_molecules[j][1], me.solutes_cog_molecules[j][2], me.solutes_cog_molecules[j][3]);
                                                }

#ifdef DEBUG
                                                std::cout << "      threads " << g << ": solutes (COG) gathered" << std::endl;
#endif // DEBUG

                                                CalculateEigenValuesVectors(&activeFrames[f], &me, &eigenValVec[f], &centerOfGeometry[f]);

#ifdef DEBUG
                                                std::cout << "      threads " << g << ": COG correction performed" << std::endl;
#endif // DEBUG

                                            }
                                        });					
                                    }
                                    for (int g = 0; g < me.num_thread_real; g++)
                                    {
                                        myThreads[g].join();
#ifdef DEBUG
                                        std::cout << "--->joining thread: " << g << std::endl;
#endif // DEBUG
                                    }
#ifdef DEBUG
                                    std::cout << "--->all thread finished successfully" << std::endl;
#endif // DEBUG

                                    CorrectionSignEigenValuesVectors(&me, &eigenValVec);

                                    for (int g = 0; g < me.num_thread_real; g++)
                                    {
                                        myThreads[g] = std::thread([&activeFrames, &me, g, &eigenValVec, &centerOfGeometry]()
                                        {
                                            for (int f = (g * me.num_frame_per_thread); f < ((g + 1) * me.num_frame_per_thread); f++)
                                            {
                                                CorrectRotationSolute(&activeFrames[f], &me, &eigenValVec[f]);

                                                //STEP 2c: gather the solvent with respect to the COG of solutes
                                                COGGatherer(&activeFrames[f], me.solvent_atoms_first, me.solvent_atoms_last, me.solvent_size, me.solvent_dimension_expansion);

#ifdef DEBUG
                                                std::cout << "      threads " << g << ": solvent (COG) gathered" << std::endl;
#endif // DEBUG

                                                //STEP 2d: correct COG to (0,0,0)
                                                if (me.cog_correction)
                                                {
                                                    for (int j = 0; j < me.atomrecords; j++)
                                                    {                                                
                                                        activeFrames[f].x[j] = activeFrames[f].x[j] - activeFrames[f].solute_cog_x;
                                                    }
                                                    for (int j = 0; j < me.atomrecords; j++)
                                                    {                                                
                                                        activeFrames[f].y[j] = activeFrames[f].y[j] - activeFrames[f].solute_cog_y;
                                                    }
                                                    for (int j = 0; j < me.atomrecords; j++)
                                                    {                                                
                                                        activeFrames[f].z[j] = activeFrames[f].z[j] - activeFrames[f].solute_cog_z;
                                                    }
                                                }

#ifdef DEBUG
                                                std::cout << "      threads " << g << ": COG correction performed" << std::endl;
#endif // DEBUG

                                            }
                                        });					
                                    }
                                    for (int g = 0; g < me.num_thread_real; g++)
                                    {
                                        myThreads[g].join();
#ifdef DEBUG
                                        std::cout << "--->joining thread: " << g << std::endl;
#endif // DEBUG
                                    }

                                }
                                catch (const std::exception &e) {
                                    std::wcout << "\nEXCEPTION: " << e.what() << std::endl;
                                }

                                delete[] myThreads;

                                //STEP 3
                                //write out all frames sequentially
                                try
                                {
                                    if (outfileThread.joinable())
                                    {
                                        outfileThread.join();
                                    }
                                }
                                catch (const std::exception &e) {
                                    std::wcout << "\nEXCEPTION (join): " << e.what() << std::endl;
                                }

                                try                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
                                {
                                    activeFramesCopy = activeFrames;
                                    outfileThread = std::thread([&activeFramesCopy, activeFrame_count, &me, &outfile, writeAtomRecordsCount]()
                                    {
                                        for (int g = 0; g < activeFrame_count; g++)
                                        {
                                            //std::cout << activeFramesCopy[g].time << " " << activeFramesCopy[g].timestep << std::endl;

                                            //WriteOutFrame(&activeFramesCopy[g], outfile, writeAtomRecordsCount, me.outformat, me.cog_write);
                                            WriteOutFrame(&activeFramesCopy[g], outfile, &me);
                                        }
                                    });
                                }
                                catch (const std::exception &e) {
                                    std::wcout << "\nEXCEPTION (write): " << e.what() << std::endl;
                                }
                                //reset activeFrame counter to 0
                                //i.e. get the next n frames and process those
                                activeFrame_counter = 0;
                            }
                            //clear all existing block content
                            titleBlock = "";
                            timestepBlock = "";
                            //increment frame_counter by 1, for keeping track which frame is being read in.
                            frame_counter += 1;
                        }

                        //reset booleans and increment frame count by 1
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
                        }
                        //what to do if currently in a POSITION block
                        if (isPositionBlock)
                        {
                            if (processThisFrame)
                            {
                                PositionBlockParser(me, line, positionBlock_counter, &prefix, &x, &y, &z);
                            }
                            //a frame is not process the data of the previous frame is retained,
                            //but the data of the first atom is changed to match that of the current frame
                            else if (!processThisFrame && positionBlock_counter == 0)
                            {
                                PositionBlockParser(me, line, positionBlock_counter, &prefix, &x, &y, &z);
                            }
                            positionBlock_counter += 1;
                        }
                        //what to do if currently in a GENBOX block
                        if (isGenboxBlock)
                        {
                            GenboxParser(&currentFrame, genBox_counter, line);
                            genBox_counter += 1;
                        }
                    }
                }
            }

            //what to do with remaining frames if end-of-file is reach
            if (activeFrame_counter > 0)
            {

#ifdef DEBUG
                std::cout << activeFrame_counter << " left-overs found" << std::endl;  
#endif // DEBUG

                std::vector<std::thread> workers;
                int remainder = activeFrame_count % me.num_thread_real;
                int perThread = activeFrame_count / me.num_thread_real;   

                //divide as equally
                for (int g = 0; g < me.num_thread_real; g++)
                {
                    workers.push_back(std::thread([&activeFrames, &me, g, perThread, &eigenValVec, &centerOfGeometry]()
                    {
                        for (int f = (g * perThread); f < ( (g + 1) * perThread); f++)
                        {
                            //STEP 2a: //solute gathering
                            SoluteGatherer(&activeFrames[f], f, &me);
                            
                            for (unsigned int j = 0; j < me.solutes_cog_molecules.size(); j++)
                            {
                                //STEP 2b: gather the non-solvent with respect to the COG of solutes
                                COGGatherer(&activeFrames[f], me.solutes_cog_molecules[j][0],  me.solutes_cog_molecules[j][1], me.solutes_cog_molecules[j][2], me.solutes_cog_molecules[j][3]);
                            }

                            CalculateEigenValuesVectors(&activeFrames[f], &me, &eigenValVec[f], &centerOfGeometry[f]);
                        }
                    }));
                }
                for (int g = 0; g < workers.size(); g++)
                {
                    workers[g].join();
                }
                workers.clear();
                CorrectionSignEigenValuesVectors(&me, &eigenValVec);
                for (int g = 0; g < me.num_thread_real; g++)
                {
                    workers.push_back(std::thread([&activeFrames, &me, g, perThread, &eigenValVec, &centerOfGeometry]()
                    {
                        for (int f = (g * perThread); f < ( (g + 1) * perThread); f++)
                        {

                            CorrectRotationSolute(&activeFrames[f], &me, &eigenValVec[f]);

                            if (!me.solvent_skip)
                            {
                                //STEP 2c: gather the solvent with respect to the COG of solutes
                                COGGatherer(&activeFrames[f], me.solvent_atoms_first, me.solvent_atoms_last, me.solvent_size, me.solvent_dimension_expansion);
                            }

                            //STEP 2d: correct COG to (0,0,0)
                            if (me.cog_correction)
                            {
                                for (int j = 0; j < me.atomrecords; j++)
                                {                                                
                                    activeFrames[f].x[j] = activeFrames[f].x[j] - activeFrames[f].solute_cog_x;
                                }
                                for (int j = 0; j < me.atomrecords; j++)
                                {                                                
                                    activeFrames[f].y[j] = activeFrames[f].y[j] - activeFrames[f].solute_cog_y;
                                }
                                for (int j = 0; j < me.atomrecords; j++)
                                {                                                
                                    activeFrames[f].z[j] = activeFrames[f].z[j] - activeFrames[f].solute_cog_z;
                                }
                            }
                        }
                    }));
                }
                for (int g = 0; g < workers.size(); g++)
                {
                    workers[g].join();
                }
                workers.clear();

                //chew on whatever remains
                for (int g = (perThread * me.num_thread_real); g < remainder; g++)
                {
                    workers.push_back(std::thread([&activeFrames, &me, g, &eigenValVec, &centerOfGeometry]()
                    {
                        //STEP 2a: //solute gathering
                        SoluteGatherer(&activeFrames[g], g, &me);

                        for (unsigned int j = 0; j < me.solutes_cog_molecules.size(); j++)
                        {
                            //STEP 2b: gather the non-solvent with respect to the COG of solutes
                            COGGatherer(&activeFrames[g], me.solutes_cog_molecules[j][0],  me.solutes_cog_molecules[j][1], me.solutes_cog_molecules[j][2], me.solutes_cog_molecules[j][3]);
                        }

                        CalculateEigenValuesVectors(&activeFrames[g], &me, &eigenValVec[g], &centerOfGeometry[g]);
                    }));
                }
                for (int g = 0; g < workers.size(); g++)
                {
                    workers[g].join();
                }
                workers.clear(); 
                CorrectionSignEigenValuesVectors(&me, &eigenValVec);
                for (int g = (perThread * me.num_thread_real); g < remainder; g++)
                {
                    workers.push_back(std::thread([&activeFrames, &me, g, &eigenValVec, &centerOfGeometry]()
                    {
                        CorrectRotationSolute(&activeFrames[g], &me, &eigenValVec[g]);

                        if (!me.solvent_skip)
                        {
                            //STEP 2c: gather the solvent with respect to the COG of solutes
                            COGGatherer(&activeFrames[g], me.solvent_atoms_first, me.solvent_atoms_last, me.solvent_size, me.solvent_dimension_expansion);
                        }

                        //STEP 2d: correct COG to (0,0,0)
                        if (me.cog_correction)
                        {
                            for (int j = 0; j < me.atomrecords; j++)
                            {                                                
                                activeFrames[g].x[j] = activeFrames[g].x[j] - activeFrames[g].solute_cog_x;
                            }
                            for (int j = 0; j < me.atomrecords; j++)
                            {                                                
                                activeFrames[g].y[j] = activeFrames[g].y[j] - activeFrames[g].solute_cog_y;
                            }
                            for (int j = 0; j < me.atomrecords; j++)
                            {                                                
                                activeFrames[g].z[j] = activeFrames[g].z[j] - activeFrames[g].solute_cog_z;
                            }
                        }
                    }));
                }
                for (int g = 0; g < workers.size(); g++)
                {
                    workers[g].join();
                }
                workers.clear();


                if (outfileThread.joinable())
                {
                    outfileThread.join();
                }

                //STEP 3
                //write out all frames sequentially
                for (int g = 0; g < activeFrame_counter; g++)
                {
                    //WriteOutFrame(&activeFrames[g], outfile, writeAtomRecordsCount, me.outformat, me.cog_write);
                    WriteOutFrame(&activeFrames[g], outfile, &me);
                }
            }
            else
            {
                if (outfileThread.joinable())
                {
                    outfileThread.join();
                }
            }
        }
        //clean up

#ifdef DEBUG
        std::cout << "left-overs processed" << std::endl;  
#endif // DEBUG

        file.close();
    }

    //clean up and report
    if (outfileThread.joinable())
    {
        outfileThread.join();
    }

    outfile.close();
    auto end = std::chrono::system_clock::now();
    auto diff = end - start;
    std::cout << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
    return 0;
}