//#define DEBUG
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
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm> 
#include <string>
#include <chrono>
#include <thread>
#include <vector>
#include <future>
#include <numeric>
#include <complex>
#include <cmath>
#include "gzstream.h"
#include <libxml/parser.h>
#include <libxml/tree.h>                     
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

//a single frame containing coordinates data
struct frame {
    int frame_id;
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

    //shared parameters
    int atomrecords;    
    double distance_cut_off;
    int verbosity;

    //dmovX parameters
    std::vector<std::vector<int>> dmov_atomgroup;
    std::vector<std::vector<int>> dmov_angles;
    std::vector<std::vector<int>> dmov_dihedral_angles;
    std::vector<std::vector<int>> dmov_atomgroup_inclusion_type;
};

//Jacobi
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

//core algorithm 
void FrameCalculation(frame *framedata, params *me, std::vector<std::vector<double>> *eigenValVec, 
    std::vector<double> *angles, std::vector<double> *dihedrals, std::vector<std::vector<double>> *cog){ 
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
            Jacobi(matrix, eigenvalues, eigenvectors, 100000);

            //output matrix

            frame_eigenValVec[h][0] = eigenvalues[0];
            frame_eigenValVec[h][1] = eigenvalues[1];
            frame_eigenValVec[h][2] = eigenvalues[2];
            frame_eigenValVec[h][3] = eigenvectors[0];
            frame_eigenValVec[h][4] = eigenvectors[1];
            frame_eigenValVec[h][5] = eigenvectors[2];
            frame_eigenValVec[h][6] = eigenvectors[3];
            frame_eigenValVec[h][7] = eigenvectors[4];
            frame_eigenValVec[h][8] = eigenvectors[5];
            frame_eigenValVec[h][9] = eigenvectors[6];
            frame_eigenValVec[h][10] = eigenvectors[7];
            frame_eigenValVec[h][11] = eigenvectors[8];

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
        //if one would prefer to use push_back
        //than use mutex.lock() and mutex.unlock()
        //http://www.baptiste-wicht.com/2012/03/cp11-concurrency-tutorial-part-2-protect-shared-data/
        *eigenValVec = frame_eigenValVec;
        *cog = frame_cog;

        //std::cout << "test" << std::endl;
        //for (int h = 0; h < me->dmov_atomgroup.size(); h++)
        //{ 
        //    std::cout << std::setw(10) << std::setprecision(4) << (*eigenValVec)[h][2] << " " << std::setw(10) << std::setprecision(4) << (*eigenValVec)[h][2] << " " << std::setw(10) << std::setprecision(4) << (*eigenValVec)[h][2] << " " << std::endl;
        //}
        //std::cout << std::endl;

        //calculated the angles or dihedrals
        //declare some shared variables
        std::vector<double> a_1 (4);
        std::vector<double> a_2 (4);
        std::vector<double> a_3 (4);
        std::vector<double> a_4 (4);
        std::vector<double> d_1 (4);
        std::vector<double> d_2 (4);
        std::vector<double> d_s (4);
        double a_r = 0;
        double sign = 0;
        std::vector<double> frame_angles;

        //calculated the angles
        for (int i = 0; i < me->dmov_angles.size(); i++)
        {
            a_1 = (*eigenValVec)[me->dmov_angles[i][0]];
            a_2 = (*eigenValVec)[me->dmov_angles[i][1]];
            a_r = acos((a_1[0]*a_2[0] +  a_1[1]*a_2[1] +  a_1[2]*a_2[2])
                /
                sqrt((a_1[0]*a_1[0] + a_1[1]*a_1[1] + a_1[2]*a_1[2])*(a_2[0]*a_2[0] + a_2[1]*a_2[1] + a_2[2]*a_2[2])));
            frame_angles.push_back(a_r);
        }
        *angles = frame_angles;
        frame_angles.clear();

        //calculated dihedrals angles
        for (int i = 0; i < me->dmov_dihedral_angles.size(); i++)
        {
            if (me->dmov_dihedral_angles[i][2] == -1 && me->dmov_dihedral_angles[i][3] == -1)
            {
                //first vector
                a_1 = (*eigenValVec)[me->dmov_dihedral_angles[i][0]];

                //second vector calculated from center of geometries of two groups
                a_2 = (*cog)[me->dmov_dihedral_angles[i][0]];            
                a_4 = (*cog)[me->dmov_dihedral_angles[i][1]];
                a_2[0] = a_4[0] - a_2[0];
                a_2[1] = a_4[1] - a_2[1];
                a_2[2] = a_4[2] - a_2[2];

                //third vector
                a_3 = (*eigenValVec)[me->dmov_dihedral_angles[i][1]];            

                //cross-product
                d_1[0] = a_1[1] * a_2[2] - a_1[2] * a_2[1];
                d_1[1] = a_1[2] * a_2[0] - a_1[0] * a_2[2];
                d_1[2] = a_1[0] * a_2[1] - a_1[1] * a_2[0];

                d_2[0] = a_2[1] * a_3[2] - a_2[2] * a_3[1];
                d_2[1] = a_2[2] * a_3[0] - a_2[0] * a_3[2];
                d_2[2] = a_2[0] * a_3[1] - a_2[1] * a_3[0];

                //determine the sign for the angle
                d_s[0] = d_1[1] * d_2[2] - d_1[2] * d_2[1];
                d_s[1] = d_1[2] * d_2[0] - d_1[0] * d_2[2];
                d_s[2] = d_1[0] * d_2[1] - d_1[1] * d_2[0];
                sign = d_s[0] * a_3[0] +  d_s[1] * a_3[1] +  d_s[2] * a_3[2];

                //calculate angle in degrees with correct sign
                a_r = ((sign > 0) - (sign < 0)) * acos(
                    (d_1[0] * d_2[0] + d_1[1] * d_2[1] + d_1[2] * d_2[2]) / (
                    sqrt(d_1[0] * d_1[0] + d_1[1] * d_1[1] + d_1[2] * d_1[2]) *
                    sqrt(d_2[0] * d_2[0] + d_2[1] * d_2[1] + d_2[2] * d_2[2])
                    )) * (180/pi);
            }
            else if (me->dmov_dihedral_angles[i][2] != -1 && me->dmov_dihedral_angles[i][3] != -1)
            {
                //calculate everything using center of geometry
                //first vector calculated from center of geometries of two groups
                a_1 = (*cog)[me->dmov_dihedral_angles[i][0]];            
                a_2 = (*cog)[me->dmov_dihedral_angles[i][1]];
                a_1[0] = a_2[0] - a_1[0];
                a_1[1] = a_2[1] - a_1[1];
                a_1[2] = a_2[2] - a_1[2];

                //second vector calculated from center of geometries of two groups
                a_2 = (*cog)[me->dmov_dihedral_angles[i][1]];            
                a_4 = (*cog)[me->dmov_dihedral_angles[i][2]];
                a_2[0] = a_4[0] - a_2[0];
                a_2[1] = a_4[1] - a_2[1];
                a_2[2] = a_4[2] - a_2[2];

                //third vector calculated from center of geometries of two groups
                a_3 = (*cog)[me->dmov_dihedral_angles[i][2]];            
                a_4 = (*cog)[me->dmov_dihedral_angles[i][3]];
                a_3[0] = a_4[0] - a_3[0];
                a_3[1] = a_4[1] - a_3[1];
                a_3[2] = a_4[2] - a_3[2];         

                //cross-product
                d_1[0] = a_1[1] * a_2[2] - a_1[2] * a_2[1];
                d_1[1] = a_1[2] * a_2[0] - a_1[0] * a_2[2];
                d_1[2] = a_1[0] * a_2[1] - a_1[1] * a_2[0];

                d_2[0] = a_2[1] * a_3[2] - a_2[2] * a_3[1];
                d_2[1] = a_2[2] * a_3[0] - a_2[0] * a_3[2];
                d_2[2] = a_2[0] * a_3[1] - a_2[1] * a_3[0];

                //determine the sign for the angle
                d_s[0] = d_1[1] * d_2[2] - d_1[2] * d_2[1];
                d_s[1] = d_1[2] * d_2[0] - d_1[0] * d_2[2];
                d_s[2] = d_1[0] * d_2[1] - d_1[1] * d_2[0];
                sign = d_s[0] * a_3[0] +  d_s[1] * a_3[1] +  d_s[2] * a_3[2];

                //calculate angle in degrees with correct sign
                a_r = ((sign > 0) - (sign < 0)) * acos(
                    (d_1[0] * d_2[0] + d_1[1] * d_2[1] + d_1[2] * d_2[2]) / (
                    sqrt(d_1[0] * d_1[0] + d_1[1] * d_1[1] + d_1[2] * d_1[2]) *
                    sqrt(d_2[0] * d_2[0] + d_2[1] * d_2[1] + d_2[2] * d_2[2])
                    )) * (180/pi);
            }
            else
            {
                a_r = 0;
            }
            frame_angles.push_back(a_r);
        }
        *dihedrals = frame_angles;
}

//write out collected data
void WriteOut(frame *framedata, gz::ogzstream &outfile, params *me, 
              std::vector<std::vector<double>> *eigenValVec, std::vector<std::vector<double>> *centerOfGeometry, 
              std::vector<double> *angles, std::vector<double> *dihedralAngles, bool verbose) {
                  if (!outfile)
                  {
                      std::cerr << "cannot open output file" << "\n" << std::flush;
                  }
                  else 
                  {      
                      //std::cout << std::right << std::setw(16) << std::setprecision(0) << framedata->timestep << " " << std::setw(16) << std::fixed << std::setprecision(4) << framedata->time << " " 
                      //    << std::setprecision(9)
                      //    << framedata->prefix[0] << " " << std::setw(14) << framedata->x[0] << " " << std::setw(14) << framedata->y[0] << " " << std::setw(14) << framedata->z[0] << std::endl;
                      outfile << std::setw(16) << std::setprecision(0) << framedata->timestep << " " 
                          << std::setw(16) << std::fixed << std::setprecision(4) << framedata->time << " "
                          << std::setprecision(9);
                      for (int h = 0; h < me->dmov_atomgroup.size(); h++)
                      {
                          outfile << std::setw(16) << (*eigenValVec)[h][2] << " "
                              << std::setw(16) << (*eigenValVec)[h][9] << " "
                              << std::setw(16) << (*eigenValVec)[h][10] << " "
                              << std::setw(16) << (*eigenValVec)[h][11] << " ";
                      }
                      for (int h = 0; h < me->dmov_atomgroup.size(); h++)
                      {
                          outfile << std::setw(16) << (*centerOfGeometry)[h][0] << " ";
                          outfile << std::setw(16) << (*centerOfGeometry)[h][1] << " ";
                          outfile << std::setw(16) << (*centerOfGeometry)[h][2] << " ";
                      }
                      for (int h = 0; h < me->dmov_angles.size(); h++)
                      {
                          outfile << std::setw(16) << (*angles)[h] << " ";
                      }
                      for (int h = 0; h < me->dmov_dihedral_angles.size(); h++)
                      {
                          outfile << std::setw(16) << (*dihedralAngles)[h] << " ";
                      }
                      outfile << "\n" << std::flush;
                  }
}

//get value in xpath
char* XPathGetText(std::string xpath_query, xmlXPathContextPtr xpathCtx) 
{
    xmlXPathObjectPtr xpathObj = xmlXPathEvalExpression(BAD_CAST xpath_query.c_str(), xpathCtx);
    return (char*)xmlNodeGetContent(xpathObj->nodesetval->nodeTab[0]);	
}

//reads a references file for topology information that is needed for writting out
void TrcReferenceFrame(std::vector<std::string> *prefix, std::string trc_reference) {

    std::ifstream infile(trc_reference);

    //check if it is possible to read file
    if (!infile)
    {
        std::cerr << "cannot open output file" << "\n" << std::flush;
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
}

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

        xpathCtx->node = xpathObj->nodesetval->nodeTab[i];
        //get topology info
        this_params->atomrecords = atoi(XPathGetText("./topology/atom_records_count", xpathCtx));
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
        //get number of frames per thread
        this_params->num_frame_per_thread = atoi(XPathGetText("./variables/frames_per_thread", xpathCtx));
        //get number of thread
        this_params->num_thread = atoi(XPathGetText("./variables/number_of_threads", xpathCtx));
        if (this_params->num_thread <= 0)
        {
            this_params->num_thread = std::thread::hardware_concurrency();
        }
        //get multiplier for thread
        if (atoi(XPathGetText("./variables/number_of_threads_multiplier", xpathCtx)) <= 1)
        {            
            this_params->num_thread_real = this_params->num_thread;
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

        //get dmov analysis atom groups
        xpathCtx->node = xpathObj->nodesetval->nodeTab[i];
        xpath = "./analysis/dmov/atom_groups/atom_group";
        xmlXPathObjectPtr atom_groups = xmlXPathEvalExpression(BAD_CAST xpath.c_str(), xpathCtx);
        for (int j = 0; j < atom_groups->nodesetval->nodeNr; j++)
        {
            std::vector<int> temp;
            xpathCtx->node = atom_groups->nodesetval->nodeTab[j];
            xpath = "./selection";
            xmlXPathObjectPtr selections = xmlXPathEvalExpression(BAD_CAST xpath.c_str(), xpathCtx);
            for (int k = 0; k < selections->nodesetval->nodeNr; k++)
            {
                xpathCtx->node = selections->nodesetval->nodeTab[k];
                int first_atom = atoi(XPathGetText("./@first_atom", xpathCtx));
                int last_atom = atoi(XPathGetText("./@last_atom", xpathCtx));
                for (int l = first_atom - 1; l < last_atom; l++)
                {
                    temp.push_back(l);
                }
            };
            std::sort(temp.begin(), temp.end());
            this_params->dmov_atomgroup.push_back(temp);
            xmlXPathFreeObject(selections);
        };
        xmlXPathFreeObject(atom_groups);

        //get dmov analysis angles 
        xpathCtx->node = xpathObj->nodesetval->nodeTab[i];
        xpath = "./analysis/dmov/angles/angle";
        xmlXPathObjectPtr angles = xmlXPathEvalExpression(BAD_CAST xpath.c_str(), xpathCtx);
        for (int j = 0; j < angles->nodesetval->nodeNr; j++)
        {
            std::vector<int> temp;
            xpathCtx->node = angles->nodesetval->nodeTab[j];
            temp.push_back(atoi(XPathGetText("./@first_group", xpathCtx))  -1);
            temp.push_back(atoi(XPathGetText("./@second_group", xpathCtx)) -1);
            this_params->dmov_angles.push_back(temp);
        };
        xmlXPathFreeObject(angles);

        //get dmov analysis dihedrals 
        xpathCtx->node = xpathObj->nodesetval->nodeTab[i];
        xpath = "./analysis/dmov/dihedrals/dihedral";
        xmlXPathObjectPtr dihedrals = xmlXPathEvalExpression(BAD_CAST xpath.c_str(), xpathCtx);
        for (int j = 0; j < dihedrals->nodesetval->nodeNr; j++)
        {
            std::vector<int> temp;
            xpathCtx->node = dihedrals->nodesetval->nodeTab[j];
            temp.push_back(atoi(XPathGetText("./@first_group", xpathCtx)) -1);
            temp.push_back(atoi(XPathGetText("./@second_group", xpathCtx))-1);
            temp.push_back(atoi(XPathGetText("./@third_group", xpathCtx)) -1);
            temp.push_back(atoi(XPathGetText("./@fourth_group", xpathCtx))-1);
            this_params->dmov_dihedral_angles.push_back(temp);
        };
        xmlXPathFreeObject(dihedrals);

        //read output block
        xpathCtx->node = xpathObj->nodesetval->nodeTab[i];
        this_params->outfilename = XPathGetText("./output/filename_prefix", xpathCtx);
        this_params->output_fragment_skipframes = atoi(XPathGetText("./output/frame_interval", xpathCtx));
        this_params->output_fragment_skiptime = atof(XPathGetText("./output/time_interval", xpathCtx));
        this_params->verbosity = atoi(XPathGetText("./output/verbosity", xpathCtx));        
    }

    this_params->ref_coords[0] = 0;
    this_params->ref_coords[1] = 0;
    this_params->ref_coords[2] = 0;

    //clean up
    xmlXPathFreeObject(xpathObj);
    xmlXPathFreeContext(xpathCtx); 
    xmlFreeDoc(doc);
}

//parses the genbox block
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

//parses the position block
void PositionBlockParser(params &me, std::string &line, int positionBlock_counter, std::vector<std::string> *prefix, std::vector<double> *x, std::vector<double> *y, std::vector<double> *z) 
{
    if (me.informat == "trc")
    {
        (*x)[positionBlock_counter] = std::stod(line.substr(0,15));
        (*y)[positionBlock_counter] = std::stod(line.substr(15,15));
        (*z)[positionBlock_counter] = std::stod(line.substr(30,15));
    }
    else if (me.informat == "cnf")
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

    //determine output filename
    gz::ogzstream outfile;
    outfile.open((me.outfilename + ".gz").c_str(), std::ios_base::out);
    outfile << std::fixed;

    //some streamstring holder
    std::ostringstream iss ("");

    //output thread
    std::thread outfileThread = std::thread([](){return 0;});

#ifdef DEBUG
    std::cout << "--->output file and thread defined" << std::endl;
    std::cout << "#dmov: domain movement" << "\n" << std::flush;
    std::cout << "#------------------------------------------------------------------------------------------------------------------------" << "\n" << std::flush;
    std::cout << "# defined input parameters\n" << std::flush;
    std::cout << std::left;
    std::cout << "#  input format                    : " << me.informat << "\n" << std::flush;
    std::cout << "#  input files                     : " << "\n" << std::flush;
    for (unsigned int i = 0; i < me.input_files.size(); i++)
    {
        std::cout << "#                                    " << me.input_files[i] << "\n" << std::flush;
    }
    std::cout << "#  frames per thread               : " << me.num_frame_per_thread << "\n" << std::flush;
    std::cout << "#  number of threads               : " << me.num_thread << "\n" << std::flush;
    std::cout << "#  real number of threads          : " << me.num_thread_real << "\n" << std::flush;
    std::cout << "#  output filename prefix          : " << me.outfilename << "\n" << std::flush;
    std::cout << "#  output every n frame(s)         : " << me.output_fragment_skipframes << "\n" << std::flush;
    std::cout << "#  output every n picoseconds      : " << me.output_fragment_skiptime << "\n" << std::flush;
    std::cout << "#  reference coordinates           : " << me.ref_coords[0] << " " << me.ref_coords[1] << " " << me.ref_coords[2] << "\n" << std::flush;
    std::cout << "#  number of solutes               : " << me.solute_count << "\n" << std::flush;
    std::cout << "#  solutes atom numbering          : " << "\n" << std::flush;
    for (unsigned int i = 0; i < me.solute_molecules.size(); i++)
    {
        std::cout << "#                                    " << me.solute_molecules[i][0] << " " << me.solute_molecules[i][1] << " " << me.solute_molecules[i][2] << " " << me.solute_molecules[i][3] << "\n" << std::flush;
    }

    std::cout << "#  analysis atoms selection        : " << "\n" << std::flush;
    for (unsigned int i = 0; i < me.dmov_atomgroup.size(); i++)
    {
        std::cout << "#      atom group " << std::setw(6) << i + 1 << "           :\n" << std::flush;
        iss.str("");
        for (unsigned int k = 0; k < me.dmov_atomgroup[i].size(); k++)
        {
            iss << me.dmov_atomgroup[i][k] + 1 << " ";
            if (iss.str().length() > 80)
            {
                std::cout << "#                                    "<< iss.str() << "\n" << std::flush;
                iss.str("");
            }
        }
        if (iss.str().length() > 0)
        {
            std::cout << "#                                    "<< iss.str() << "\n" << std::flush;
            iss.str("");
        }
    }

    std::cout << "#  analysis angles                 : " << "\n" << std::flush;
    for (unsigned int i = 0; i < me.dmov_angles.size(); i++)
    {
        std::cout << "#      angle " << std::setw(6) << i + 1 << "                : " << me.dmov_angles[i][0] + 1 << " " << me.dmov_angles[i][1] + 1 << "\n" << std::flush;
    }

    std::cout << "#  analysis dihedral angles        : " << "\n" << std::flush;
    for (unsigned int i = 0; i < me.dmov_dihedral_angles.size(); i++)
    {
        std::cout << "#      dihedral angle " << std::setw(6) << i + 1 << "       : " << me.dmov_dihedral_angles[i][0] + 1 << " " << me.dmov_dihedral_angles[i][1] + 1
            << " " << me.dmov_dihedral_angles[i][2] + 1 << " " << me.dmov_dihedral_angles[i][3] + 1 << "\n" << std::flush;
    }

    std::cout << "#  trc refernce file               : " << me.trc_reference << "\n" << std::flush;
#endif // DEBUG

    outfile << "#dmov: domain movement" << "\n" << std::flush;
    outfile << "#------------------------------------------------------------------------------------------------------------------------" << "\n" << std::flush;
    outfile << "# defined input parameters\n" << std::flush;
    outfile << std::left;
    outfile << "#  input format                    : " << me.informat << "\n" << std::flush;
    outfile << "#  input files                     : " << "\n" << std::flush;
    for (unsigned int i = 0; i < me.input_files.size(); i++)
    {
        outfile << "#                                    " << me.input_files[i] << "\n" << std::flush;
    }
    outfile << "#  frames per thread               : " << me.num_frame_per_thread << "\n" << std::flush;
    outfile << "#  number of threads               : " << me.num_thread << "\n" << std::flush;
    outfile << "#  real number of threads          : " << me.num_thread_real << "\n" << std::flush;
    outfile << "#  output filename prefix          : " << me.outfilename << "\n" << std::flush;
    outfile << "#  output every n frame(s)         : " << me.output_fragment_skipframes << "\n" << std::flush;
    outfile << "#  output every n picoseconds      : " << me.output_fragment_skiptime << "\n" << std::flush;
    outfile << "#  reference coordinates           : " << me.ref_coords[0] << " " << me.ref_coords[1] << " " << me.ref_coords[2] << "\n" << std::flush;
    outfile << "#  number of solutes               : " << me.solute_count << "\n" << std::flush;
    outfile << "#  solutes atom numbering          : " << "\n" << std::flush;
    for (unsigned int i = 0; i < me.solute_molecules.size(); i++)
    {
        outfile << "#                                    " << me.solute_molecules[i][0] << " " << me.solute_molecules[i][1] << " " << me.solute_molecules[i][2] << " " << me.solute_molecules[i][3] << "\n" << std::flush;
    }

    outfile << "#  analysis atoms selection        : " << "\n" << std::flush;
    for (unsigned int i = 0; i < me.dmov_atomgroup.size(); i++)
    {
        outfile << "#      atom group " << std::setw(6) << i + 1 << "           :\n" << std::flush;
        iss.str("");
        for (unsigned int k = 0; k < me.dmov_atomgroup[i].size(); k++)
        {
            iss << me.dmov_atomgroup[i][k] + 1 << " ";
            if (iss.str().length() > 80)
            {
                outfile << "#                                    "<< iss.str() << "\n" << std::flush;
                iss.str("");
            }
        }
        if (iss.str().length() > 0)
        {
            outfile << "#                                    "<< iss.str() << "\n" << std::flush;
            iss.str("");
        }
    }

    outfile << "#  analysis angles                 : " << "\n" << std::flush;
    for (unsigned int i = 0; i < me.dmov_angles.size(); i++)
    {
        outfile << "#      angle " << std::setw(6) << i + 1 << "                : " << me.dmov_angles[i][0] << " " << me.dmov_angles[i][1] + 1 << "\n" << std::flush;
    }

    outfile << "#  analysis dihedral angles        : " << "\n" << std::flush;
    for (unsigned int i = 0; i < me.dmov_dihedral_angles.size(); i++)
    {
        outfile << "#      dihedral angle " << std::setw(6) << i + 1 << "       : " << me.dmov_dihedral_angles[i][0] + 1 << " " << me.dmov_dihedral_angles[i][1] + 1
            << " " << me.dmov_dihedral_angles[i][2] + 1 << " " << me.dmov_dihedral_angles[i][3] + 1 << "\n" << std::flush;
    }

    outfile << "#  trc refernce file               : " << me.trc_reference << "\n" << std::flush;
    outfile << "#------------------------------------------------------------------------------------------------------------------------" << "\n" << std::flush;
    outfile << std::right;

    //write the column headers
    outfile << "#" << std::setw(9) << "timestep" << " " 
        << std::setw(16) << "time" << " "; 
    for (int i = 0; i < me.dmov_atomgroup.size(); i++)
    {
        iss.str("");
        iss << "e" << i+1 << "_v";
        outfile << std::setw(16) << iss.str() << " ";
        iss.str("");
        iss << "e" << i+1 << "_x";
        outfile << std::setw(16) << iss.str() << " ";
        iss.str("");
        iss << "e" << i+1 << "_y";
        outfile << std::setw(16) << iss.str() << " ";
        iss.str("");
        iss << "e" << i+1 << "_z";
        outfile << std::setw(16) << iss.str() << " ";
    }

    for (int h = 0; h < me.dmov_atomgroup.size(); h++)
    {
        iss.str("");
        iss << "cog_" << h+1 << "_x";;
        outfile << std::setw(16) << iss.str() << " ";
        iss.str("");
        iss << "cog_" << h+1 << "_y";;
        outfile << std::setw(16) << iss.str() << " ";
        iss.str("");
        iss << "cog_" << h+1 << "_z";;
        outfile << std::setw(16) << iss.str() << " ";
    }
    for (int i = 0; i < me.dmov_angles.size(); i++)
    {
        iss.str("");
        iss << "angle_" << i + 1;
        outfile << std::setw(16) << iss.str() << " ";        
    }
    for (int i = 0; i < me.dmov_dihedral_angles.size(); i++)
    {
        iss.str("");
        iss  << "dihedral_" << i + 1;
        outfile << std::setw(16) << iss.str() << " ";        
    }
    outfile << "\n" << std::flush;

    //holds the active frames
    int activeFrame_count = me.num_thread_real * me.num_frame_per_thread;
    std::vector<frame> activeFrames (activeFrame_count);
    std::vector<frame> activeFramesCopy (activeFrame_count);

    //hold the dmov frame eigenvalues and eigenvectors data
    std::vector<std::vector<std::vector<double>>> eigenValVec (activeFrame_count);
    std::vector<std::vector<std::vector<double>>> eigenValVecCopy (activeFrame_count);
    std::vector<std::vector<double>> angles (activeFrame_count);
    std::vector<std::vector<double>> anglesCopy (activeFrame_count);
    std::vector<std::vector<double>> dihedralAngles (activeFrame_count);
    std::vector<std::vector<double>> dihedralAnglesCopy (activeFrame_count);
    std::vector<std::vector<std::vector<double>>> centerOfGeometry (activeFrame_count);
    std::vector<std::vector<std::vector<double>>> centerOfGeometryCopy (activeFrame_count);

    //performance log - starting time
    auto start = std::chrono::system_clock::now();

    //only use the first TITLE block found
    bool firstPass = true;

    //frame counter which can be used for skipping frames
    int frame_counter = 0, s = 0, molecule_start = 0, molecule_end = 0, px = 0, py = 0, pz = 0;
    double frame_time = 0;
    bool processThisFrame = false;

    //dimension expansion factor
    int periodic_copies = 1;

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

    //process the trajectory files sequentially
    std::cout << "processing trajectory files: " << std::endl;
    for (unsigned int i = 0; i < me.input_files.size(); i++)
    {
        std::cout << "  " << me.input_files[i] << std::endl;

        //define file to read
        gz::igzstream file(me.input_files[i].c_str());

        //boolean specifying current active block
        bool isTitleBlock = false , isTimestepBlock = false, isPositionBlock = false, isGenboxBlock = false;

        //content holder for the blocks
        std::string titleBlock(""), timestepBlock(""), positionBlock(""), genboxBlock("");

        //counters
        int positionBlock_counter = 0, activeFrame_counter = 0, genBox_counter = 0;

        //check if it is possible to read file
        if (!file)
        {
            std::cerr << "cannot open output file" << "\n" << std::flush;
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
                    //reset if end of the block is found
                    else if (line.substr(0,3) == "END")
                    {    
                        //what to do if end of a TITLE block
                        if (isTitleBlock)
                        {
                            //TITLE block
                            if (me.outformat == "cnf" && firstPass)
                            {
                                outfile << "TITLE" << "\n" << std::flush;
                                outfile << titleBlock;
                                outfile << "END" << "\n" << std::flush;
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
                                currentFrame.time = std::stod(timestepBlock.substr(19,19));
                                currentFrame.timestep = std::stol(timestepBlock.substr(0,18));
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

                            //until n frames have been read in
                            if (processThisFrame)
                            {
                                currentFrame.frame_id = frame_counter;
                                activeFrames[activeFrame_counter] = currentFrame;
                                activeFrame_counter += 1;
                            }

                            //if after processing the GENBOX block, n frames are stored than do something
                            if (activeFrame_counter == activeFrame_count)
                            {
#ifdef DEBUG
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
                                        //PAST THE RIGHT VARIABLES, OBJECTS, WHATEVER....
                                        myThreads[g] = std::thread([&activeFrames, &me, g, &eigenValVec, &angles, &dihedralAngles, &centerOfGeometry]() {
                                            for (int f = (g * me.num_frame_per_thread); f < ((g + 1) * me.num_frame_per_thread); f++)
                                            {
                                                //DO SOMETHING
                                                FrameCalculation(&activeFrames[f], &me, &eigenValVec[f], &angles[f], &dihedralAngles[f], &centerOfGeometry[f]);
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
                                    eigenValVecCopy = eigenValVec;
                                    anglesCopy = angles;
                                    dihedralAnglesCopy = dihedralAngles;
                                    centerOfGeometryCopy = centerOfGeometry;

                                    //PAST THE RIGHT VARIABLES, OBJECTS, WHATEVER....
                                    outfileThread = std::thread([&activeFramesCopy, &me, &outfile, &eigenValVecCopy, &anglesCopy, &dihedralAnglesCopy, &centerOfGeometryCopy]()
                                    {
                                        //DO SOMETHING
                                        for (int g = 0; g < activeFramesCopy.size(); g++)
                                        {
                                            WriteOut(&activeFramesCopy[g], outfile, &me, &eigenValVecCopy[g], &centerOfGeometryCopy[g], &anglesCopy[g], &dihedralAnglesCopy[g], false); 
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
                    //PAST THE RIGHT VARIABLES, OBJECTS, WHATEVER....
                    workers.push_back(std::thread([&activeFrames, &me, g, &eigenValVec, &angles, &dihedralAngles, &centerOfGeometry]() {
                        for (int f = (g * me.num_frame_per_thread); f < ((g + 1) * me.num_frame_per_thread); f++)
                        {
                            //DO SOMETHING
                            FrameCalculation(&activeFrames[f], &me, &eigenValVec[f], &angles[f], &dihedralAngles[f], &centerOfGeometry[f]);
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
                    //PAST THE RIGHT VARIABLES, OBJECTS, WHATEVER....
                    workers.push_back(std::thread([&activeFrames, &me, g, &eigenValVec, &angles, &dihedralAngles, &centerOfGeometry]() {
                        for (int f = (g * me.num_frame_per_thread); f < ((g + 1) * me.num_frame_per_thread); f++)
                        {
                            //DO SOMETHING
                            FrameCalculation(&activeFrames[f], &me, &eigenValVec[f], &angles[f], &dihedralAngles[f], &centerOfGeometry[f]);
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

                for (int g = 0; g < activeFrame_counter; g++)
                {                    
                    //DO SOMETHING
                    WriteOut(&activeFrames[g], outfile, &me, &eigenValVec[g], &centerOfGeometry[g], &angles[g], &dihedralAngles[g], false); 
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