#define EIGEN_NO_DEBUG  
#define EIGEN_MPL2_ONLY
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include "FrameGeometry.h"
#include "Structs.h"
#include "gzstream.h"
#include <iomanip>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SVD>

FrameGeometry::FrameGeometry(void)
{
}


FrameGeometry::~FrameGeometry(void)
{
}


void FrameGeometry::Gather(void)
{
}

//parse a given line in the GENBOX block
void FrameGeometry::TrajectoryGenboxBlockLineParser(Structs::FrameGeometric *currentFrame, int &genBox_counter, std::string &line) 
{
    switch (genBox_counter)
    {
    case 0:
        currentFrame->boxtype = std::stoi(line);
        break;
    case 1:
        currentFrame->box_length.x() = std::stod(line.substr(0,15));
        currentFrame->box_length.y() = std::stod(line.substr(15,15));
        currentFrame->box_length.z() = std::stod(line.substr(30,15));
        break;
    case 2:
        currentFrame->box_angle.x() = std::stod(line.substr(0,15));
        currentFrame->box_angle.y() = std::stod(line.substr(15,15));
        currentFrame->box_angle.z() = std::stod(line.substr(30,15));
        break;
    case 3:
        currentFrame->box_3.x() = std::stod(line.substr(0,15));
        currentFrame->box_3.y() = std::stod(line.substr(15,15));
        currentFrame->box_3.z() = std::stod(line.substr(30,15));
        break;
    case 4:
        currentFrame->box_4.x() = std::stod(line.substr(0,15));
        currentFrame->box_4.y() = std::stod(line.substr(15,15));
        currentFrame->box_4.z() = std::stod(line.substr(30,15));
        break;
    }
}

//parse a given line in the POSITION block
void FrameGeometry::TrajectoryPositionBlockLineParser(Structs::InputParametersFrameout &params, std::string &line, int &positionBlock_counter, std::vector<std::string> *prefix, Eigen::MatrixXd *coordinates) 
{
    if (params.informat == "trc")
    {
        (*coordinates)(0,positionBlock_counter) = std::stod(line.substr(0,15));
        (*coordinates)(1,positionBlock_counter) = std::stod(line.substr(15,15));
        (*coordinates)(2,positionBlock_counter) = std::stod(line.substr(30,15));
    }
    if (params.informat == "cnf")
    {
        (*prefix)[positionBlock_counter] = line.substr(0,24);
        (*coordinates)(0,positionBlock_counter) = std::stod(line.substr(25,15));
        (*coordinates)(1,positionBlock_counter) = std::stod(line.substr(40,15));
        (*coordinates)(2,positionBlock_counter) = std::stod(line.substr(55,15));
    }
}

//reads a references file for topology information that is needed for writting out
//files in CNF and PDB compatible format
void FrameGeometry::TrcReferenceFrame(std::vector<std::string> *prefix, std::string trc_reference) {
    std::ifstream infile(trc_reference);

    //check if it is possible to read file
    if (!infile)
    {
        std::cerr << "cannot open reference input file" << "\n";
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

//write out the data in either CNF or PDB compatible format
//this is not for production... code has to be written differently 
//FIGURE OUT WHAT THE QUICKFIX MEANS
void FrameGeometry::WriteOutFrame(Structs::FrameGeometric *framedata, gz::ogzstream &outfile, Structs::InputParametersFrameout *me) {
    //check if it is possible to read file
    if (!outfile)
    {
        std::cerr << "cannot open output file" << "\n";
    }
    else 
    {
        //THE QUICK FIX IS FOR WRITING OUT CORRECTLY WHEN NO SOLVENT IS PRESENT OR GATHERED
        //it just determines the range of everything that is not solvent and writes that out
        std::vector<int> temp_fix; 
        //this get the first and last atom of each solute
        for (int i = 0; i < me->solute_count; i++)
        {
            temp_fix.push_back(me->solute_molecules(0,i)-1);
            temp_fix.push_back(me->solute_molecules(1,i)-1);
        }

        for (int i = 0; i < me->solute_cog_molecules.cols(); i++)
        {
            temp_fix.push_back(me->solute_cog_molecules(0,i)-1);
            temp_fix.push_back(me->solute_cog_molecules(1,i)-1);
        }

        for (int i = 0; i < me->ion_molecules.cols(); i++)
        {
            temp_fix.push_back(me->ion_molecules(0,i)-1);
            temp_fix.push_back(me->ion_molecules(1,i)-1);
        }
        //quickfix

        if (me->outformat=="cnf" || me->outformat=="trc")
        {
            outfile << "TIMESTEP" << "\n";
            outfile << " " << std::setw(17) << std::setprecision(0) << framedata->timestep << " " << std::setw(19) << std::fixed << std::setprecision(9) << framedata->time << "\n";
            outfile << "END" << "\n";
            outfile << "POSITION" << "\n";	
            outfile << std::fixed << std::setprecision(9);

            //quickfix for what?
            //write out the coordinates based on skipping solvent or not
            for (int i = 0; i < me->atomrecords; i++)
            {
                if (me->solvent_skip)
                {
                    //if the current atom record is within the range of the molecule atom numbers
                    for (unsigned int ii = 0; ii < temp_fix.size(); ii+=2)
                    {
                        if (i>=temp_fix[ii] && i<=temp_fix[ii+1])
                        {
                            outfile << framedata->prefix[i] << " " << std::setw(14) << framedata->coordinates(0,i) << " " << std::setw(14) << framedata->coordinates(1,i) << " " << std::setw(14) << framedata->coordinates(2,i) << "\n";
                        }
                    }
                }
                else
                {                    
                    outfile << framedata->prefix[i] << " " << std::setw(14) << framedata->coordinates(0,i) << " " << std::setw(14) << framedata->coordinates(1,i) << " " << std::setw(14) << framedata->coordinates(2,i) << "\n";
                }
            }
            //quickfix

            outfile << "END" << "\n";
            outfile << "GENBOX" << "\n";
            outfile << " " << framedata->boxtype << "\n";
            outfile << std::fixed << std::setprecision(9);
            outfile << " " << std::setw(14) << framedata->box_length.x() << " " << std::setw(14) << framedata->box_length.y() << " " << std::setw(14) << framedata->box_length.z() << "\n";
            outfile << " " << std::setw(14) << framedata->box_angle.x() << " " << std::setw(14) << framedata->box_angle.y() << " " << std::setw(14) << framedata->box_angle.z() << "\n";
            outfile << " " << std::setw(14) << framedata->box_3.x() << " " << std::setw(14) << framedata->box_3.y() << " " << std::setw(14) << framedata->box_3.z() << "\n";
            outfile << " " << std::setw(14) << framedata->box_4.x() << " " << std::setw(14) << framedata->box_4.y() << " " << std::setw(14) << framedata->box_4.z() << "\n";
            outfile << "END" << "\n" << std::flush;;
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
            outfile << std::fixed << std::setprecision(9);
            outfile << "CRYST1" 
                << std::setw(9) << std::setprecision(3) << framedata->box_length.x()*10
                << std::setw(9) << std::setprecision(3) << framedata->box_length.y()*10
                << std::setw(9) << std::setprecision(3) << framedata->box_length.z()*10
                << std::setw(7) << std::setprecision(2) << framedata->box_angle.x()
                << std::setw(7) << std::setprecision(2) << framedata->box_angle.y()
                << std::setw(7) << std::setprecision(2) << framedata->box_angle.z()
                << " P 1           1"
                << "\n";
            outfile << std::fixed << std::setprecision(9);
            for (int i = 0; i < me->atomrecords; i++)
            {
                //quickfix
                if (me->solvent_skip)
                {
                    for (unsigned int ii = 0; ii < temp_fix.size(); ii+=2)
                    {
                        if (i>=temp_fix[ii] && i<=temp_fix[ii+1])
                        {
                            //20130819                                                     

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
                                << std::right << std::fixed 
                                << std::setw(8) << std::setprecision(3) << framedata->coordinates(0,i) * 10
                                << std::setw(8) << std::setprecision(3) << framedata->coordinates(1,i) * 10
                                << std::setw(8) << std::setprecision(3) << framedata->coordinates(2,i) * 10
                                << std::setw(6) << std::setprecision(2) << 1.0
                                << std::setw(6) << std::setprecision(2) << 1.0 << "\n"                
                                << std::flush;
                            //20130819   
                        }
                    }
                }
                else
                {
                    //20130819                            
                    //if (!(framedata->coordinates.col(i).array() == 0).all()) {
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
                            << std::right << std::fixed 
                            << std::setw(8) << std::setprecision(3) << framedata->coordinates(0,i) * 10
                            << std::setw(8) << std::setprecision(3) << framedata->coordinates(1,i) * 10
                            << std::setw(8) << std::setprecision(3) << framedata->coordinates(2,i) * 10
                            << std::setw(6) << std::setprecision(2) << 1.0
                            << std::setw(6) << std::setprecision(2) << 1.0 << "\n"                
                            << std::flush;
                        //20130819   
                    //}
                }
            }
            if (me->cog_write)
            {
                outfile << "HETATM" << std::setw(5) << me->atomrecords+2 << "  ZN   ZN A9999    "
                    << std::fixed
                    << std::setw(8) << std::setprecision(3) << framedata->solute_cog.x() * 10
                    << std::setw(8) << std::setprecision(3) << framedata->solute_cog.y() * 10
                    << std::setw(8) << std::setprecision(3) << framedata->solute_cog.z() * 10
                    << "  1.00  0.00          ZN  \n";
            }
            outfile << "ENDMDL" << "\n";
        }
    }
}

//WARNING: never include any ions or solvent in the rotation
void FrameGeometry::CorrectionRotational(Structs::FrameGeometric *framedata, Structs::FrameGeometric *ref_framedata, Structs::InputParametersFrameout *me) {

    //where to stop applying correction
    int last_atom = 0;

    //apply correction to the last of the solutes
    if ((me->solute_molecules.cols() > 0 && me->solute_cog_molecules.cols() > 0))
    {
        last_atom = me->solute_cog_molecules(1,me->solute_cog_molecules.cols()-1);
    }
    //apply correction to the essentail solutes
    else         if ((me->solute_molecules.cols() > 0))
    {
        last_atom = me->solute_molecules(1,me->solute_molecules.cols()-1);
    }
    //if ions are not present, apply correction based on atoms up to the first solvent
    else
    {
        return;
        last_atom = me->solvent_molecules(0,0) - 1;
    }

    for (int i = 0; i < 10; i++)
    {
        Eigen::JacobiSVD<Eigen::Matrix<double, 3, Eigen::Dynamic>> svd (framedata->coordinates.leftCols(last_atom) * ref_framedata->coordinates.leftCols(last_atom).transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
        auto R  = svd.matrixV()*svd.matrixU().transpose();
        for (int j = 0; j < framedata->coordinates.cols(); j++)
        {
            //framedata->coordinates.col(i) = ref_framedata->solute_cog + (R * (framedata->coordinates.col(i) - framedata->solute_cog));
            framedata->coordinates.col(j) = (R * framedata->coordinates.col(j));
        }   
    } 
}