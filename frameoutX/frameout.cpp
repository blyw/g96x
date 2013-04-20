//#define DEBUG
#define POSIX

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
#include <Eigen/Dense>

//a single frame containing coordinates data
struct frame {
    long timestep;
    double time;
    std::vector<std::string> prefix;
    Eigen::Matrix<double, 3, Eigen::Dynamic> coordinates;
    Eigen::Vector3d solute_cog;
    int boxtype;
    Eigen::Vector3d box_length;
    Eigen::Vector3d box_angle;
    Eigen::Vector3d box_3;
    Eigen::Vector3d box_4;
    Eigen::Vector3i init_shift;
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
    Eigen::Matrix<int, 4, Eigen::Dynamic> solute_molecules;
    Eigen::Matrix<int, 5, Eigen::Dynamic> solutes_cog_molecules;

    //topology parameters - ions
    Eigen::Matrix<int, 4, Eigen::Dynamic>  ions_molecules;
    bool ions_skip;

    //topology parameters - solvent
    Eigen::Matrix<int, 4, 1>  solvent_molecules;
    bool solvent_skip;

    //input paramters
    std::string informat;
    std::string trc_reference;
    std::vector<std::string> input_files;

    //shared parameters
    int atomrecords;    
    double distance_cut_off;
    int verbosity;

    //frameoutX parameters
    Eigen::Vector3d ref_coords;
    bool cog_write;
    bool cog_correction;
    bool gather;
};

//simple gathering of a specified atom in frame with respect to a given coordinate
void FirstAtomBasedBoxShifter(frame *framedata, int atom_number, params *me) {
    if (me->gather)
    {
        double distance_shortest = 1E20;

        Eigen::Vector3d coords(0,0,0);
        Eigen::Vector3i min_shift(0,0,0);

        //dimension expansion factor
        int periodic_copies = 2;
        if (me->solute_count > 0)
        {
            periodic_copies = (me->solute_molecules(2,0)+1)*2;
        }
        //shift the coordinates of the specified atom in the original framedata based on the previous frame shift
        if (!(framedata->init_shift(0) == 0 && framedata->init_shift(1) == 0 && framedata->init_shift(2) == 0))
        {
            framedata->coordinates.colwise() += (framedata->init_shift).cast<double>().cwiseProduct(framedata->box_length);
        }

        //build search grid
        Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> grid;
        grid.resize(3, (periodic_copies*2+1)*(periodic_copies*2+1)*(periodic_copies*2+1));
        int grid_counter = 0;

        for (int x = 0-periodic_copies; x <= periodic_copies; x++)
        {
            for (int y = 0-periodic_copies; y <= periodic_copies; y++)
            {
                for (int z = 0-periodic_copies; z <= periodic_copies; z++)
                {
                    grid(0,grid_counter) = x;
                    grid(1,grid_counter) = y;
                    grid(2,grid_counter) = z;
                    grid_counter += 1;
                }
            }
        }

        //gather the first atom fo the frame
        for (int i = 0; i < grid.cols(); i++)
        {
            double distance = 
                (framedata->coordinates.col(atom_number) + framedata->box_length.cwiseProduct((grid.col(i)).cast<double>()) - me->ref_coords).dot(framedata->coordinates.col(atom_number) + framedata->box_length.cwiseProduct((grid.col(i)).cast<double>()) - me->ref_coords);
            if (distance < distance_shortest)
            {
                distance_shortest = distance;
                min_shift = grid.col(i);
            }
        }

        //shift again if first atom was shifted
        if (!( min_shift(0) == 0 && min_shift(1) == 0 && min_shift(2) == 0))
        {
            framedata->coordinates.colwise() += framedata->box_length.cwiseProduct(min_shift.cast<double>());
        }

        //maybe useful for debugging?
        framedata->init_shift = min_shift;
    }
}

//solute gathering i.e. cluster algorithm (nearest neighbour)
void SoluteGatherer(frame *framedata, int frameId, params *me){ // std::vector<std::vector<int>> solute_molecules, int solute_count, int frameId, int periodic_copies) {	
    //init variables
    double distance_shortest;
    Eigen::Vector3i min_shift(0,0,0);

    //cut-off is derived from (longest bond)^2
    double cut_off = me->distance_cut_off;
    int molecule_start, molecule_end, molecule_start_previous, periodic_copies;

    //for center of geometry of all solute molecules calculation
    Eigen::Vector3d coordinates_sum(0,0,0);

    //build search grid
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> grid;
    grid.resize(3, (periodic_copies*2+1)*(periodic_copies*2+1)*(periodic_copies*2+1));
    int grid_counter = 0;
    for (int x = 0-periodic_copies; x <= periodic_copies; x++)
    {
        for (int y = 0-periodic_copies; y <= periodic_copies; y++)
        {
            for (int z = 0-periodic_copies; z <= periodic_copies; z++)
            {
                grid(0,grid_counter) = x;
                grid(1,grid_counter) = y;
                grid(2,grid_counter) = z;
                grid_counter += 1;
            }
        }
    }

    //cross-reference gathering
    for (int s = 0; s < me->solute_count; s++)
    {
        //define first atom, last atom and number of periodic copies of a solute molecule
        molecule_start = me->solute_molecules(0,s);
        molecule_end = me->solute_molecules(1,s);
        periodic_copies = me->solute_molecules(2,s);

        //the first atom in a solute molecule which is not the first solute molecule is 
        //gather with respect to all atoms of the previous solute molecule
        if (s > 0)
        { 
            distance_shortest = 1E20;
            min_shift.Zero();
            //molecule_start_previous is re-defined each time a solute molecule is gathered
            //molecule_start - 1 --> array index of the first atom of the current solute molecule
            //molecule_start - 2 --> array index of the last atom of the previous solute molecule
            for (int i = (molecule_start - 2); i >= (molecule_start_previous - 1); i--)
            {	
                //gather
                for (int g = 0; g < grid.cols(); g++)
                {
                    //defined cutt-off use to prevent unnecessary looping through the full list
                    if (distance_shortest <= cut_off)
                    {
                        break;
                    }

                    double distance = 
                        (framedata->coordinates.col(molecule_start - 1) + framedata->box_length.cwiseProduct(grid.col(g).cast<double>()) - framedata->coordinates.col(i)).dot(
                        framedata->coordinates.col(molecule_start - 1) + framedata->box_length.cwiseProduct(grid.col(g).cast<double>()) - framedata->coordinates.col(i));
                    if (distance < distance_shortest)
                    {
                        distance_shortest = distance;
                        min_shift = grid.col(g);
                    }
                }
            }
            //shift the coordinates of the atom in the original frame
            framedata->coordinates.col(molecule_start - 1) += framedata->box_length.cwiseProduct(min_shift.cast<double>());
        }

        //dimension sum for cog 
        coordinates_sum += framedata->coordinates.col(molecule_start - 1);

        //double xval, yval, zval;
        //gather solute molecule
        for (int i = molecule_start; i < molecule_end; i++)
        {				
            //atom dependent variables that should be reset for each atom
            distance_shortest = 1E20;
            min_shift.Zero();
            //reverse search from coordinate closest to already-in-list atoms of current solute molecule
            for (int j = (i - 1); j >= (molecule_start - 1); j--)
            {
                //gather
                for (int g = 0; g < grid.cols(); g++)
                {
                    //defined cutt-off use to prevent unnecessary looping through the full list
                    if (distance_shortest <= cut_off)
                    {
                        break;
                    }

                    double distance = 
                        (framedata->coordinates.col(molecule_start - 1) + framedata->box_length.cwiseProduct(grid.col(g).cast<double>()) - framedata->coordinates.col(i)).dot(
                        framedata->coordinates.col(molecule_start - 1) + framedata->box_length.cwiseProduct(grid.col(g).cast<double>()) - framedata->coordinates.col(i));
                    if (distance < distance_shortest)
                    {
                        distance_shortest = distance;
                        min_shift = grid.col(g);
                    }
                }
            }

            //shift the coordinates of the atom in the original frame	
            framedata->coordinates.col(i) += framedata->box_length.cwiseProduct(min_shift.cast<double>());

            //dimension sum for cog 
            coordinates_sum += framedata->coordinates.col(i);
        }
        //done gathering for a solute molecule
        //set reference for gathering next solute molecule
        molecule_start_previous = molecule_start;
    }

    if (me->solute_count > 0)
    {
        //use this center of geometry for the entire frame
        framedata->solute_cog = coordinates_sum / molecule_end;
    }
    else
    {
        //use this center of geometry for the entire frame        
        framedata->solute_cog = framedata->coordinates.col(0);
    }
}

//gather solvent molecules with respect to COG of solutes
void COGGatherer(frame *framedata, int first_atom, int last_atom, int molecule_size, int periodic_copies) {
    double distance_shortest = 1E20;
    Eigen::Vector3i min_shift(0,0,0);
    //cut-off is derived from (longest bond)^2
    //double cut_off = me->distance_cut_off;
    double cut_off = 0.3 * 0.3;

    //build search grid
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> grid;
    grid.resize(3, (periodic_copies*2+1)*(periodic_copies*2+1)*(periodic_copies*2+1));
    int grid_counter = 0;

    for (int x = 0-periodic_copies; x <= periodic_copies; x++)
    {
        for (int y = 0-periodic_copies; y <= periodic_copies; y++)
        {
            for (int z = 0-periodic_copies; z <= periodic_copies; z++)
            {
                grid(0, grid_counter) = x;
                grid(1, grid_counter) = y;
                grid(2, grid_counter) = z;
                grid_counter += 1;
            }
        }
    }

    for (int i = first_atom - 1; i < last_atom; i+=molecule_size)
    {
        distance_shortest = 1E20;
        min_shift.Zero();

        //gather the first atom of the solvent molecule
        for (int g = 0; g < grid.cols(); g++)
        {
            //defined cutt-off use to prevent unnecessary looping through the full list
            double distance = 
                (framedata->coordinates.col(i) + framedata->box_length.cwiseProduct(grid.col(g).cast<double>()) - framedata->solute_cog).dot(
                framedata->coordinates.col(i) + framedata->box_length.cwiseProduct(grid.col(g).cast<double>()) - framedata->solute_cog);
            if (distance < distance_shortest)
            {
                distance_shortest = distance;
                min_shift = grid.col(g);
            }
        }

        //shift the coordinates of the first atom in the original framedata
        framedata->coordinates.col(i) += framedata->box_length.cwiseProduct(min_shift.cast<double>());

        ////gather other atoms of the solvent molecule
        for (int j = (i + 1); j < (i + molecule_size); j++)
        {
            //reset for the each additional atom in solvent molecule
            distance_shortest = 1E20;
            min_shift.Zero();

            for (int k = i; k < j; k++)
            {   

                for (int g = 0; g < grid.cols(); g++)
                {                   
                    //defined cutt-off use to prevent unnecessary looping through the full list
                    if (distance_shortest <= cut_off)
                    {
                        break;
                    }

                    double distance = 
                        (framedata->coordinates.col(j) + framedata->box_length.cwiseProduct(grid.col(g).cast<double>()) - framedata->coordinates.col(k)).dot(
                        framedata->coordinates.col(j) + framedata->box_length.cwiseProduct(grid.col(g).cast<double>()) - framedata->coordinates.col(k));
                    if (distance < distance_shortest)
                    {
                        distance_shortest = distance;
                        min_shift = grid.col(g);
                    }
                }
            }
            //	//shift atom correctly
            framedata->coordinates.col(j) += framedata->box_length.cwiseProduct(min_shift.cast<double>());
        }
    }
}

//write out the data in either CNF or PDB compatible format
//this is not for production... code has to be written differently 
void WriteOutFrame(frame *framedata, gz::ogzstream &outfile, params *me) {
    //check if it is possible to read file
    if (!outfile)
    {
        std::cerr << "cannot open output file" << "\n";
    }
    else 
    {
        //quickfix for what???
        std::vector<int> temp_fix; 
        //this get the first and last atom of each solute
        for (int i = 0; i < me->solute_count; i++)
        {
            temp_fix.push_back(me->solute_molecules(0,i)-1);
            temp_fix.push_back(me->solute_molecules(1,i)-1);
        }
        //for (int i = 0; i < me->solutes_cog_molecules.size(); i++)
        for (int i = 0; i < me->solutes_cog_molecules.cols(); i++)
        {
            temp_fix.push_back(me->solutes_cog_molecules(0,i)-1);
            temp_fix.push_back(me->solutes_cog_molecules(1,i)-1);
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
                    for (int ii = 0; ii < temp_fix.size(); ii+=2)
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
                                << std::setw(8) << std::setprecision(3) << framedata->coordinates(0,i) * 10
                                << std::setw(8) << std::setprecision(3) << framedata->coordinates(1,i) * 10
                                << std::setw(8) << std::setprecision(3) << framedata->coordinates(2,i) * 10
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
                        << std::setw(8) << std::setprecision(3) << framedata->coordinates(0,i) * 10
                        << std::setw(8) << std::setprecision(3) << framedata->coordinates(1,i) * 10
                        << std::setw(8) << std::setprecision(3) << framedata->coordinates(2,i) * 10
                        << std::setw(6) << std::setprecision(2) << 1.0
                        << std::setw(6) << std::setprecision(2) << 1.0 << "\n"                
                        ;
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

char* XPathGetText(std::string xpath_query, xmlXPathContextPtr xpathCtx) 
{
    xmlXPathObjectPtr xpathObj = xmlXPathEvalExpression(BAD_CAST xpath_query.c_str(), xpathCtx);
    return (char*)xmlNodeGetContent(xpathObj->nodesetval->nodeTab[0]);	
}

//reads a references file for topology information that is needed for writting out
//files in CNF and PDB compatible format
void TrcReferenceFrame(std::vector<std::string> *prefix, std::string trc_reference) {
    std::ifstream infile(trc_reference);

    //check if it is possible to read file
    if (!infile)
    {
        std::cerr << "cannot open input file" << "\n";
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
        this_params->solvent_molecules(0,0) = atoi(XPathGetText("./topology/solvent/@first_atom", xpathCtx));
        //get last atom
        this_params->solvent_molecules(1,0) = atoi(XPathGetText("./topology/solvent/@last_atom", xpathCtx));
        //skip or keep solvent
        this_params->solvent_skip = atoi(XPathGetText("./topology/solvent/@skip", xpathCtx));
        //get number of atoms
        this_params->solvent_molecules(2,0) = atoi(XPathGetText("./topology/solvent/number_of_atoms", xpathCtx));
        //get dimension expansion factor
        this_params->solvent_molecules(3,0) = atoi(XPathGetText("./topology/solvent/dimension_expansion", xpathCtx));

#ifdef DEBUG
        std::cout << "got solvent parameters" << std::endl;
#endif // DEBUG

        //get solute (cog) definitions
        xpath = "./topology/solutes_cog/solute";
        xmlXPathObjectPtr COGSolute = xmlXPathEvalExpression(BAD_CAST xpath.c_str(), xpathCtx);
        //loops through solvents
        //std::cout << xmlNodeGetContent(xpathObj->nodesetval->nodeTab[0]) << std::endl;
        //get all COG solutes
        this_params->solutes_cog_molecules.resize(5, COGSolute->nodesetval->nodeNr); 
        for (int j = 0; j < COGSolute->nodesetval->nodeNr; j++)
        {
            xpathCtx->node = COGSolute->nodesetval->nodeTab[j];
            //get first atom of COG solute
            this_params->solutes_cog_molecules(0,i) = atoi(XPathGetText("./@first_atom", xpathCtx));
            //get last atom of COG solute
            this_params->solutes_cog_molecules(1,i) = atoi(XPathGetText("./@last_atom", xpathCtx));
            //get number of atom of COG solute
            this_params->solutes_cog_molecules(2,i) = atoi(XPathGetText("./number_of_atoms", xpathCtx));
            //get dimension expansion factor
            this_params->solutes_cog_molecules(3,i) = atoi(XPathGetText("./dimension_expansion", xpathCtx));
            //skip or keep COG solute
            this_params->solutes_cog_molecules(4,i) = atoi(XPathGetText("./@skip", xpathCtx));
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
        this_params->solute_molecules.resize(4, solutes->nodesetval->nodeNr);

        std::cout << solutes->nodesetval->nodeNr << std::endl;
        std::cout << this_params->solute_molecules.cols() << std::endl;
        std::cout << this_params->solute_molecules.rows() << std::endl;
        std::cout << this_params->solute_molecules.size() << std::endl;

        for (int j = 0; j < solutes->nodesetval->nodeNr; j++)
        {
            xpathCtx->node = solutes->nodesetval->nodeTab[j];
            //get first atom of solute
            this_params->solute_molecules(0,j) = atoi(XPathGetText("./@first_atom", xpathCtx));
            //get last atom of solute
            this_params->solute_molecules(1,j) = atoi(XPathGetText("./@last_atom", xpathCtx));
            //get dimension expansion factor
            this_params->solute_molecules(2,j) = atoi(XPathGetText("./dimension_expansion", xpathCtx));
            //skip or keep solute
            this_params->solute_molecules(3,j) = atoi(XPathGetText("./@skip", xpathCtx));
        };
        this_params->solute_count = this_params->solute_molecules.cols();
        xmlXPathFreeObject(solutes);

#ifdef DEBUG
        std::cout << "got solutes parameters" << std::endl;
#endif // DEBUG

        //read the variables section
        xpathCtx->node = xpathObj->nodesetval->nodeTab[i];
        //write out cog or not
        this_params->cog_write = atoi(XPathGetText("./analysis/frameout/cog_write", xpathCtx));
        //correct by shifting cog or not
        this_params->cog_correction = atoi(XPathGetText("./analysis/frameout/cog_correction", xpathCtx));
        //get distance cut
        this_params->distance_cut_off = atof(XPathGetText("./analysis/frameout/distance_cut_off", xpathCtx))*atof(XPathGetText("./analysis/frameout/distance_cut_off", xpathCtx));

#ifdef DEBUG
        std::cout << "got analysis parameters" << std::endl;
#endif // DEBUG

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

#ifdef DEBUG
        std::cout << "got hardware parameters" << std::endl;
#endif // DEBUG

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

void GenboxParser(frame *currentFrame, int genBox_counter, std::string line) 
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

//undecided
void PositionBlockParser(params &me, std::string &line, int positionBlock_counter, std::vector<std::string> *prefix, Eigen::MatrixXd *coordinates) 
{
    if (me.informat == "trc")
    {
        (*coordinates)(0,positionBlock_counter) = std::stod(line.substr(0,15));
        (*coordinates)(1,positionBlock_counter) = std::stod(line.substr(15,15));
        (*coordinates)(2,positionBlock_counter) = std::stod(line.substr(30,15));
    }
    if (me.informat == "cnf")
    {
        (*prefix)[positionBlock_counter] = line.substr(0,24);
        (*coordinates)(0,positionBlock_counter) = std::stod(line.substr(25,15));
        (*coordinates)(1,positionBlock_counter) = std::stod(line.substr(40,15));
        (*coordinates)(2,positionBlock_counter) = std::stod(line.substr(55,15));
    }
}

//main program
int main(int argc, char* argv[])
{
    std::string job_id = argv[1];
    std::string param_file = argv[2];

    params me;
    ParseParamsInput(&me, job_id, param_file);

    me.gather = true;

    //decide based on input, how many lines to write out in final output file
    int writeAtomRecordsCount = 0;
    if (me.solvent_skip)
    {
        int solutes_cog_count = me.solutes_cog_molecules.cols();
        if (solutes_cog_count > 0)
        {                            
            writeAtomRecordsCount= me.solutes_cog_molecules(1,solutes_cog_count-1);
        }
        else
        {
            writeAtomRecordsCount = me.solute_molecules(1,me.solute_count-1);
        }
    }
    else {
        writeAtomRecordsCount = me.atomrecords;
    }

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
    std::cout << "  reference coordinates           : " << me.ref_coords.x() << " " << me.ref_coords.y() << " " << me.ref_coords.z() << std::endl;
    std::cout << "  skip solvent                    : " << me.solvent_skip << std::endl;
    std::cout << "  solutes cog                     : " << std::endl;
    for (unsigned int i = 0; i < me.solutes_cog_molecules.cols(); i++)
    {
        std::cout << "                                    " << me.solutes_cog_molecules(0,i) << " " << me.solutes_cog_molecules(1,i) << " " << me.solutes_cog_molecules(2,i) << " " << me.solutes_cog_molecules(3,i) <<  " " << me.solutes_cog_molecules(4,i) << std::endl;
    }
    std::cout << "  number of solutes               : " << me.solute_count << std::endl;
    std::cout << "  solutes atom numbering          : " << std::endl;
    for (unsigned int i = 0; i < me.solute_molecules.cols(); i++)
    {
        std::cout << "                                    " << me.solute_molecules(0,i) << " " << me.solute_molecules(1,i) << " " << me.solute_molecules(2,i) << " " << me.solute_molecules(3,i) << std::endl;
    }
    std::cout << "  solvent cog dimension expansion : " << me.solvent_molecules(3,0) << std::endl;
    std::cout << "  solvent size                    : " << me.solvent_molecules(2,0) << std::endl;
    std::cout << "  trc refernce file               : " << me.trc_reference << std::endl;  

    //holds the active frames
    int activeFrame_count = me.num_thread_real * me.num_frame_per_thread;
    std::vector<frame> activeFrames (activeFrame_count);
    std::vector<frame> activeFramesCopy (activeFrame_count);

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
    Eigen::MatrixXd coordinates(0,0);
    coordinates.resize(3,me.atomrecords);

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
        Eigen::Vector3d init_shift(0,0,0);

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
            std::cerr << "cannot open input file" << "\n";
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
                                currentFrame.time = std::stod(timestepBlock.substr(18,20));
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
                            currentFrame.coordinates = coordinates;

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

                            //first atom which is used to define reference between frames in the code below
                            //is defined as the first atom of the first solute molecule specified in the 
                            //input file
                            int first_atom_solute = 0;
                            if (me.solute_count >  0)
                            {
                                first_atom_solute = me.solute_molecules(0,0) - 1;
                            }

                            //skip frame[0] from the reference process, because it is the first frame 
                            //and has no reference coordinates for the first atom
                            //also get the time of the first frame as reference for time-based skipping
                            if (frame_counter==0)
                            {
                                currentFrame.init_shift(0) = 0;
                                currentFrame.init_shift(1) = 0;
                                currentFrame.init_shift(2) = 0;
                                me.ref_coords.x() = coordinates(0,first_atom_solute);
                                me.ref_coords.y() = coordinates(1,first_atom_solute);
                                me.ref_coords.z() = coordinates(2,first_atom_solute);
                            }
                            else 
                            {
                                //do not combine the following lines
                                //gather should be kept separate from reference coordinates assignment
                                FirstAtomBasedBoxShifter(&currentFrame, first_atom_solute, &me);
                                me.ref_coords.x() = coordinates(0,first_atom_solute);
                                me.ref_coords.y() = coordinates(1,first_atom_solute);
                                me.ref_coords.z() = coordinates(2,first_atom_solute);
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
                                    if (me.gather)
                                    {

                                        for (int g = 0; g < me.num_thread_real; g++)
                                        {
                                            myThreads[g] = std::thread([&activeFrames, &me, g]()
                                            {
                                                for (int f = (g * me.num_frame_per_thread); f < ((g + 1) * me.num_frame_per_thread); f++)
                                                {
                                                    //STEP 2a: //solute gathering
                                                    SoluteGatherer(&activeFrames[f], f, &me);

#ifdef DEBUG
                                                    std::cout << "      threads " << g << ": solutes gathered" << std::endl;
#endif // DEBUG

                                                    for (unsigned int j = 0; j < me.solutes_cog_molecules.cols(); j++)
                                                    {
                                                        //STEP 2b: gather the non-solvent with respect to the COG of solutes
                                                        COGGatherer(&activeFrames[f], me.solutes_cog_molecules(0,j),  me.solutes_cog_molecules(1,j), me.solutes_cog_molecules(2,j), me.solutes_cog_molecules(3,j));
                                                    }

#ifdef DEBUG
                                                    std::cout << "      threads " << g << ": solutes (COG) gathered" << std::endl;
#endif // DEBUG

                                                    //STEP 2c: gather the solvent with respect to the COG of solutes
                                                    COGGatherer(&activeFrames[f], me.solvent_molecules(0,0), me.solvent_molecules(1,0), me.solvent_molecules(2,0), me.solvent_molecules(3,0));

#ifdef DEBUG
                                                    std::cout << "      threads " << g << ": solvent (COG) gathered" << std::endl;
#endif // DEBUG

                                                    //STEP 2d: correct COG to (0,0,0)
                                                    if (me.cog_correction)
                                                    {
                                                        activeFrames[f].coordinates.colwise() -= activeFrames[f].solute_cog;  
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
#ifdef DEBUG
                                        std::cout << "--->all thread finished successfully" << std::endl;
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
                                PositionBlockParser(me, line, positionBlock_counter, &prefix, &coordinates);
                            }
                            //a frame is not process the data of the previous frame is retained,
                            //but the data of the first atom is changed to match that of the current frame
                            else if (!processThisFrame && positionBlock_counter == 0)
                            {
                                PositionBlockParser(me, line, positionBlock_counter, &prefix, &coordinates);
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

                if (me.gather)
                {

                    //divide as equally
                    for (int g = 0; g < me.num_thread_real; g++)
                    {
                        //DO SOMETHING
                        workers.push_back(std::thread([&activeFrames, &me, g, perThread]()
                        {
                            for (int f = (g * perThread); f < ( (g + 1) * perThread); f++)
                            {
                                //STEP 2a: //solute gathering
                                SoluteGatherer(&activeFrames[f], f, &me);

                                for (unsigned int j = 0; j < me.solutes_cog_molecules.size(); j++)
                                {
                                    //STEP 2b: gather the non-solvent with respect to the COG of solutes
                                    COGGatherer(&activeFrames[f], me.solutes_cog_molecules(0,j),  me.solutes_cog_molecules(1,j), me.solutes_cog_molecules(2,j), me.solutes_cog_molecules(3,j));
                                }

                                if (!me.solvent_skip)
                                {
                                    //STEP 2c: gather the solvent with respect to the COG of solutes
                                    COGGatherer(&activeFrames[f], me.solvent_molecules(0,0), me.solvent_molecules(1,0), me.solvent_molecules(2,0), me.solvent_molecules(3,0));
                                }

                                //STEP 2d: correct COG to (0,0,0)
                                if (me.cog_correction)
                                {
                                    activeFrames[f].coordinates.colwise() -= activeFrames[f].solute_cog;
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
                        //DO SOMETHING
                        workers.push_back(std::thread([&activeFrames, &me, g]()
                        {
                            //STEP 2a: //solute gathering
                            SoluteGatherer(&activeFrames[g], g, &me);

                            for (unsigned int j = 0; j < me.solutes_cog_molecules.size(); j++)
                            {
                                //STEP 2b: gather the non-solvent with respect to the COG of solutes
                                COGGatherer(&activeFrames[g], me.solutes_cog_molecules(0,j),  me.solutes_cog_molecules(1,j), me.solutes_cog_molecules(2,j), me.solutes_cog_molecules(3,j));
                            }

                            if (!me.solvent_skip)
                            {
                                //STEP 2c: gather the solvent with respect to the COG of solutes
                                COGGatherer(&activeFrames[g], me.solvent_molecules(0,0), me.solvent_molecules(1,0), me.solvent_molecules(2,0), me.solvent_molecules(3,0));
                            }

                            //STEP 2d: correct COG to (0,0,0)
                            if (me.cog_correction)
                            {
                                activeFrames[g].coordinates.colwise() -= activeFrames[g].solute_cog;
                            }
                        }));
                    }
                    for (int g = 0; g < workers.size(); g++)
                    {
                        workers[g].join();
                    }
                    workers.clear();

                }

                if (outfileThread.joinable())
                {
                    outfileThread.join();
                }

                //STEP 3
                //write out all frames sequentially
                for (int g = 0; g < activeFrame_counter; g++)
                {
                    //DO SOMETHING
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