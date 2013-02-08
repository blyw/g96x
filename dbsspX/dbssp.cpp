#define DEBUG
#define POSIX

#include <iostream>
#include <fstream>
#include <sstream>
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
#include <numeric>
#include <algorithm>

//a single frame containing coordinates data
struct frame {
    int frame_id;
    long timestep;
    double time;
    //15944 SOLV  HW2    52493    1.966742392    4.000388119    1.796537431
    std::vector<std::string> prefix;
    std::vector<int> residue_number;
    std::vector<std::string> residue_name;
    std::vector<std::string> atom_type;
    std::vector<int> atom_number;
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
    int num_thread;
    int num_frame_per_thread;
    int num_thread_real;
    int output_fragment_size;
    int output_fragment_skipframes;
    double output_fragment_skiptime;
    int solute_count;
    std::vector<std::vector<int>> solute_molecules;
    std::vector<std::vector<int>> solutes_cog_molecules;
    int protein_count;
    std::vector<std::vector<int>> proteins;
    int solvent_atoms_first;
    int solvent_atoms_last;
    int solvent_dimension_expansion;
    int solvent_size;
    bool solvent_skip;
    int atomrecords;
    double distance_cut_off;
    std::string outformat;
    std::string informat;
    std::string trc_reference;
    std::string outfilename;
    std::vector<std::string> input_files;
    double ref_coords[3];
    bool cog_write;
    bool cog_correction;
    int verbosity;

};

//code checked 20130207: CORRECT, not tested?
//calculate the dihedral angle or torsion for atoms in the backbone
void CalculateBackboneDihedralsaAtan2(frame *framedata, int frame_id, params *me, std::vector<std::vector<std::vector<double>>> *backbone_dihedrals, std::vector<std::vector<int>> *backbone_dihedral_atoms)
{ 
    std::vector<double> phi;
    std::vector<double> psi;
    std::vector<std::vector<double>> proteins;

    for (int i = 0; i < me->protein_count; i++)
    { 
        int *bbda = (*backbone_dihedral_atoms)[i].data();
        double *x = framedata->x.data();
        double *y = framedata->y.data();
        double *z = framedata->z.data();

        //vector for calculating dihedral angle
        double vector_a[3] = { 0, 0, 0 };
        double vector_b[3] = { 0, 0, 0 };
        double vector_c[3] = { 0, 0, 0 };
        double vector_ab[3] = { 0, 0, 0 };
        double vector_bc[3] = { 0, 0, 0 };
        double vector_bm[3] = { 0, 0, 0 };
        //for determining the sign
        double vector_s[3] = { 0, 0, 0 };

        //loop through list of atoms
        int j = 0;
        int dihedral_atom_count = (*backbone_dihedral_atoms)[i].size();
        double b_mag = 0;

        //psi angles
        for (j = 2; j < dihedral_atom_count - 3; j+=3)
        {
            //std::cout << "phi: " << framedata->prefix[bbda[j+0]] << " " 
            //    << framedata->prefix[bbda[j+1]] << " " 
            //    << framedata->prefix[bbda[j+2]] << " " 
            //    << framedata->prefix[bbda[j+3]] << "\n";

            vector_a[0] = x[bbda[j+1]] - x[bbda[j+0]];
            vector_b[0] = x[bbda[j+2]] - x[bbda[j+1]]; 
            vector_c[0] = x[bbda[j+3]] - x[bbda[j+2]]; 

            vector_a[1] = y[bbda[j+1]] - y[bbda[j+0]];
            vector_b[1] = y[bbda[j+2]] - y[bbda[j+1]]; 
            vector_c[1] = y[bbda[j+3]] - y[bbda[j+2]]; 

            vector_a[2] = z[bbda[j+1]] - z[bbda[j+0]];
            vector_b[2] = z[bbda[j+2]] - z[bbda[j+1]]; 
            vector_c[2] = z[bbda[j+3]] - z[bbda[j+2]]; 

            b_mag = (sqrt(vector_b[0]*vector_b[0] + vector_b[1]*vector_b[1] + vector_b[2]*vector_b[2]));

            vector_bm[0] = b_mag * vector_a[0];
            vector_bm[1] = b_mag * vector_a[1];
            vector_bm[2] = b_mag * vector_a[2];

            vector_ab[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
            vector_ab[1] = vector_a[2] * vector_b[0] - vector_a[0] * vector_b[2];
            vector_ab[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];

            vector_bc[0] = vector_b[1] * vector_c[2] - vector_b[2] * vector_c[1];
            vector_bc[1] = vector_b[2] * vector_c[0] - vector_b[0] * vector_c[2];
            vector_bc[2] = vector_b[0] * vector_c[1] - vector_b[1] * vector_c[0];

            phi.push_back(atan2(
                (vector_bm[0] * vector_bc[0] + vector_bm[1] * vector_bc[1] + vector_bm[2] * vector_bc[2]),
                (vector_ab[0] * vector_bc[0] + vector_ab[1] * vector_bc[1] + vector_ab[2] * vector_bc[2]))
                * (180/3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679));

        }

        //phi angles
        for (j = 3; j < dihedral_atom_count; j+=3)
        {
            //std::cout << "psi: " << framedata->prefix[bbda[j+0]] << " " 
            //    << framedata->prefix[bbda[j+1]] << " " 
            //    << framedata->prefix[bbda[j+2]] << " " 
            //    << framedata->prefix[bbda[j+3]] << "\n";

            vector_a[0] = x[bbda[j+1]] - x[bbda[j+0]];
            vector_b[0] = x[bbda[j+2]] - x[bbda[j+1]]; 
            vector_c[0] = x[bbda[j+3]] - x[bbda[j+2]]; 

            vector_a[1] = y[bbda[j+1]] - y[bbda[j+0]];
            vector_b[1] = y[bbda[j+2]] - y[bbda[j+1]]; 
            vector_c[1] = y[bbda[j+3]] - y[bbda[j+2]]; 

            vector_a[2] = z[bbda[j+1]] - z[bbda[j+0]];
            vector_b[2] = z[bbda[j+2]] - z[bbda[j+1]]; 
            vector_c[2] = z[bbda[j+3]] - z[bbda[j+2]]; 

            b_mag = (sqrt(vector_b[0]*vector_b[0] + vector_b[1]*vector_b[1] + vector_b[2]*vector_b[2]));

            vector_bm[0] = b_mag * vector_a[0];
            vector_bm[1] = b_mag * vector_a[1];
            vector_bm[2] = b_mag * vector_a[2];

            vector_ab[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
            vector_ab[1] = vector_a[2] * vector_b[0] - vector_a[0] * vector_b[2];
            vector_ab[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];

            vector_bc[0] = vector_b[1] * vector_c[2] - vector_b[2] * vector_c[1];
            vector_bc[1] = vector_b[2] * vector_c[0] - vector_b[0] * vector_c[2];
            vector_bc[2] = vector_b[0] * vector_c[1] - vector_b[1] * vector_c[0];

            psi.push_back(atan2(
                (vector_bm[0] * vector_bc[0] + vector_bm[1] * vector_bc[1] + vector_bm[2] * vector_bc[2]),
                (vector_ab[0] * vector_bc[0] + vector_ab[1] * vector_bc[1] + vector_ab[2] * vector_bc[2]))
                * (180/3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679));
        }

        proteins.push_back(phi);
        proteins.push_back(psi);

        phi.clear();
        psi.clear();
    }
    (*backbone_dihedrals)[frame_id] = (proteins);
}


//code checked 20130208: CORRECT!
//calculate the dihedral angle or torsion for atoms that are part of the backbone
void CalculateBackboneDihedrals(frame *framedata, int frame_id, params *me, std::vector<std::vector<std::vector<double>>> *backbone_dihedrals, std::vector<std::vector<int>> *backbone_dihedral_atoms)
{ 
    //some data holders
    std::vector<double> phi;
    std::vector<double> psi;
    std::vector<double> atomnumber;
    std::vector<std::vector<double>> proteins;

    //vector for calculating dihedral angle
    double vector_a[3] = { 0, 0, 0 };
    double vector_b[3] = { 0, 0, 0 };
    double vector_c[3] = { 0, 0, 0 };
    double vector_ab[3] = { 0, 0, 0 };
    double vector_bc[3] = { 0, 0, 0 };
    //for determining the sign
    double vector_s[3] = { 0, 0, 0 };
    double sign = 0;
    //for loops
    unsigned int i, j, dihedral_atom_count = 0;
    int *bbda;
    double *x;
    double *y;
    double *z;

    //loops through each defined protein
    for (i = 0; i < me->protein_count; i++)
    { 
        bbda = (*backbone_dihedral_atoms)[i].data();
        x = framedata->x.data();
        y = framedata->y.data();
        z = framedata->z.data();

        //loop through list of atoms
        j = 0;
        dihedral_atom_count = (*backbone_dihedral_atoms)[i].size();
        sign = 0;

        //psi angles
        for (j = 2; j < dihedral_atom_count - 3; j+=4)
        {
            //std::cout << "phi: " << framedata->prefix[bbda[j+0]] << " " 
            //    << framedata->prefix[bbda[j+1]] << " " 
            //    << framedata->prefix[bbda[j+2]] << " " 
            //    << framedata->prefix[bbda[j+3]] << "\n";

            vector_a[0] = x[bbda[j+1]] - x[bbda[j+0]];
            vector_b[0] = x[bbda[j+2]] - x[bbda[j+1]]; 
            vector_c[0] = x[bbda[j+3]] - x[bbda[j+2]]; 

            vector_a[1] = y[bbda[j+1]] - y[bbda[j+0]];
            vector_b[1] = y[bbda[j+2]] - y[bbda[j+1]]; 
            vector_c[1] = y[bbda[j+3]] - y[bbda[j+2]]; 

            vector_a[2] = z[bbda[j+1]] - z[bbda[j+0]];
            vector_b[2] = z[bbda[j+2]] - z[bbda[j+1]]; 
            vector_c[2] = z[bbda[j+3]] - z[bbda[j+2]]; 

            vector_ab[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
            vector_ab[1] = vector_a[2] * vector_b[0] - vector_a[0] * vector_b[2];
            vector_ab[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];

            vector_bc[0] = vector_b[1] * vector_c[2] - vector_b[2] * vector_c[1];
            vector_bc[1] = vector_b[2] * vector_c[0] - vector_b[0] * vector_c[2];
            vector_bc[2] = vector_b[0] * vector_c[1] - vector_b[1] * vector_c[0];

            //sign
            vector_s[0] = vector_ab[1] * vector_bc[2] - vector_ab[2] * vector_bc[1];
            vector_s[1] = vector_ab[2] * vector_bc[0] - vector_ab[0] * vector_bc[2];
            vector_s[2] = vector_ab[0] * vector_bc[1] - vector_ab[1] * vector_bc[0];

            sign = vector_s[0] * vector_c[0] +  vector_s[1] * vector_c[1] +  vector_s[2] * vector_c[2];

            phi.push_back(((sign > 0) - (sign < 0)) * acos(
                (vector_ab[0] * vector_bc[0] + vector_ab[1] * vector_bc[1] + vector_ab[2] * vector_bc[2]) / (
                sqrt(vector_ab[0] * vector_ab[0] + vector_ab[1] * vector_ab[1] + vector_ab[2] * vector_ab[2]) *
                sqrt(vector_bc[0] * vector_bc[0] + vector_bc[1] * vector_bc[1] + vector_bc[2] * vector_bc[2])
                )) * (180/3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679));

            atomnumber.push_back(bbda[j+1]);
            //}

            ////phi angles
            //for (j = 3; j < dihedral_atom_count; j+=3)
            //{
            //std::cout << "psi: " << framedata->prefix[bbda[j+0]] << " " 
            //    << framedata->prefix[bbda[j+1]] << " " 
            //    << framedata->prefix[bbda[j+2]] << " " 
            //    << framedata->prefix[bbda[j+3]] << "\n";

            vector_a[0] = x[bbda[j+2]] - x[bbda[j+1]];
            vector_b[0] = x[bbda[j+3]] - x[bbda[j+2]]; 
            vector_c[0] = x[bbda[j+4]] - x[bbda[j+3]]; 

            vector_a[1] = y[bbda[j+2]] - y[bbda[j+1]];
            vector_b[1] = y[bbda[j+3]] - y[bbda[j+2]]; 
            vector_c[1] = y[bbda[j+4]] - y[bbda[j+3]]; 

            vector_a[2] = z[bbda[j+2]] - z[bbda[j+1]];
            vector_b[2] = z[bbda[j+3]] - z[bbda[j+2]]; 
            vector_c[2] = z[bbda[j+4]] - z[bbda[j+3]]; 

            vector_ab[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
            vector_ab[1] = vector_a[2] * vector_b[0] - vector_a[0] * vector_b[2];
            vector_ab[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];

            vector_bc[0] = vector_b[1] * vector_c[2] - vector_b[2] * vector_c[1];
            vector_bc[1] = vector_b[2] * vector_c[0] - vector_b[0] * vector_c[2];
            vector_bc[2] = vector_b[0] * vector_c[1] - vector_b[1] * vector_c[0];

            //sign
            vector_s[0] = vector_ab[1] * vector_bc[2] - vector_ab[2] * vector_bc[1];
            vector_s[1] = vector_ab[2] * vector_bc[0] - vector_ab[0] * vector_bc[2];
            vector_s[2] = vector_ab[0] * vector_bc[1] - vector_ab[1] * vector_bc[0];

            sign = vector_s[0] * vector_c[0] +  vector_s[1] * vector_c[1] +  vector_s[2] * vector_c[2];

            psi.push_back(((sign > 0) - (sign < 0)) * acos(
                (vector_ab[0] * vector_bc[0] + vector_ab[1] * vector_bc[1] + vector_ab[2] * vector_bc[2]) / (
                sqrt(vector_ab[0] * vector_ab[0] + vector_ab[1] * vector_ab[1] + vector_ab[2] * vector_ab[2]) *
                sqrt(vector_bc[0] * vector_bc[0] + vector_bc[1] * vector_bc[1] + vector_bc[2] * vector_bc[2])
                )) * (180/3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679));
        }

        proteins.push_back(phi);
        proteins.push_back(psi);
        proteins.push_back(atomnumber);

        phi.clear();
        psi.clear();
        atomnumber.clear();
    }
    (*backbone_dihedrals)[frame_id] = (proteins);
}

//void PrintDihedrals(frame *framedata, int frame_id, std::vector<std::vector<std::vector<double>>> *backbone_dihedrals, gz::ogzstream &outfile) 
//{
//    int k, l, m = 0;
//    std::vector<std::vector<double>> bdf = (*backbone_dihedrals)[frame_id];
//    outfile << std::fixed;
//    //std::cout << std::endl << (*backbone_dihedrals)[j].size() << std::endl;
//    for (k = 0; k < bdf.size()/3; k++)
//    {
//        //process phi-psi data pairs for each defined protein
//        for (l = k*3; l < (k+1)*3; l++)
//        {
//            outfile << std::setw(16) << framedata->time << " ";
//            for (m = 0; m < bdf[l].size(); m++)
//            {
//                outfile << std::setw(16) << std::setprecision(9) << bdf[l][m] << " ";
//            }
//            outfile << "\n" <<std::flush;
//            //std::cout << std::endl << (*backbone_dihedrals)[j][l].size() << std::endl;
//        }
//    }
//
//}

void PrintDihedrals(frame *framedata, int frame_id, std::vector<std::vector<std::vector<double>>> *backbone_dihedrals, gz::ogzstream &outfile) 
{
    unsigned int k, l, m = 0;
    std::vector<std::vector<double>> bdf = (*backbone_dihedrals)[frame_id];
    outfile << std::fixed;
    //std::cout << std::endl << (*backbone_dihedrals)[j].size() << std::endl;
    for (k = 0; k < bdf.size()/3; k++)
    {
        //process atomnumber + phi-psi data pairs for each defined protein
        l = k*3;
        for (m = 0; m < bdf[l].size(); m++)
        {
            outfile << std::setw(16) << std::setprecision(4) << framedata->time << " ";
            outfile << std::setw(16) << std::setprecision(9) << framedata->residue_number[bdf[l+2][m]] << " ";
            outfile << std::setw(16) << std::setprecision(9) << bdf[l][m] << " ";
            outfile << std::setw(16) << std::setprecision(9) << bdf[l+1][m] << " ";
            outfile << "\n" <<std::flush;
        }
        //std::cout << std::endl << (*backbone_dihedrals)[j][l].size() << std::endl;
    }

}

//void PrintDihedrals(frame *framedata, int frame_id, std::vector<std::vector<std::vector<double>>> *backbone_dihedrals, gz::ogzstream &outfile) 
//{
//    //std::cout << std::endl << (*backbone_dihedrals)[j].size() << std::endl;
//    for (int k = 0; k < (*backbone_dihedrals)[frame_id].size()/3; k++)
//    {
//        //process phi-psi data pairs for each defined protein
//        for (int l = k*3; l < (k+1)*3; l+=3)
//        {
//            for (int m = 0; m < (*backbone_dihedrals)[frame_id][l].size(); m++)
//            {
//                outfile << std::fixed;
//                outfile << std::setw(9) << std::setprecision(2) << framedata->time << " " 
//                    << std::setw(9) << framedata->residue_number[(*backbone_dihedrals)[frame_id][l+2][m+1]] 
//                    << std::setw(16) << std::setprecision(9) << (*backbone_dihedrals)[frame_id][l][m] << " " 
//                    << std::setw(16) << std::setprecision(9) <<  (*backbone_dihedrals)[frame_id][l+1][m]  
//                    << "\n" <<std::flush;;
//                //outfile << std::fixed;
//                //outfile << std::setprecision(6) << framedata->time << " " << std::setw(4) 
//                //    << framedata->residue_number[(*backbone_dihedrals)[frame_id][l+2][m+1]] << " "
//                //<< std::setw(12) << std::setprecision(2) << (*backbone_dihedrals)[frame_id][l][m]  
//                //<< " " << std::setw(12) << std::setprecision(2) << (*backbone_dihedrals)[frame_id][l+1][m] << "\n" <<std::flush;;
//            }
//            //std::cout << std::endl << (*backbone_dihedrals)[j][l].size() << std::endl;
//        }
//    }
//
//}

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

        //get the protein definition of the input file
        xpathCtx->node = xpathObj->nodesetval->nodeTab[i];
        //get protein definitions
        xpath = "./dbssp/proteins/protein";
        xmlXPathObjectPtr proteins = xmlXPathEvalExpression(BAD_CAST xpath.c_str(), xpathCtx);
        //get all solutes
        for (int j = 0; j < proteins->nodesetval->nodeNr; j++)
        {
            std::vector<int> temp;
            xpathCtx->node = proteins->nodesetval->nodeTab[j];
            //get first atom of solute
            temp.push_back(atoi(XPathGetText("./@first_atom", xpathCtx)));
            //get last atom of solute
            temp.push_back(atoi(XPathGetText("./@last_atom", xpathCtx)));
            this_params->proteins.push_back(temp);
            temp.clear();
        };
        this_params->protein_count = this_params->proteins.size();
        xmlXPathFreeObject(proteins);

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

#ifdef DEBUG
    std::cout << "input parameters defined" << std::endl;
    std::cout << std::left;
    std::cout << "  distance cut-off                : " << me.distance_cut_off << std::endl;
    std::cout << "  input format                    : " << me.informat << std::endl;
    std::cout << "  input files                     : " << std::endl;
    for (unsigned int i = 0; i < me.input_files.size(); i++)
    {
        std::cout << "                                    " << me.input_files[i] << std::endl;
    }
    std::cout << "  frames per thread               : " << me.num_frame_per_thread << std::endl;
    std::cout << "  number of threads               : " << me.num_thread << std::endl;
    std::cout << "  real number of threads          : " << me.num_thread_real << std::endl;
    std::cout << "  output filename prefix          : " << me.outfilename << std::endl;
    std::cout << "  output every n frame(s)         : " << me.output_fragment_skipframes << std::endl;
    std::cout << "  output every n picoseconds      : " << me.output_fragment_skiptime << std::endl;
    std::cout << "  reference coordinates           : " << me.ref_coords[0] << " " << me.ref_coords[1] << " " << me.ref_coords[2] << std::endl;
    std::cout << "  number of solutes               : " << me.solute_count << std::endl;
    std::cout << "  solutes atom numbering          : " << std::endl;
    for (unsigned int i = 0; i < me.solute_molecules.size(); i++)
    {
        std::cout << "                                    " << me.solute_molecules[i][0] << " " << me.solute_molecules[i][1] << " " << me.solute_molecules[i][2] << " " << me.solute_molecules[i][3] << std::endl;
    }
    std::cout << "  trc refernce file               : " << me.trc_reference << std::endl;
#endif // DEBUG

    //holds the active frames
    int activeFrame_count = me.num_thread_real * me.num_frame_per_thread;
    std::vector<frame> activeFrames (activeFrame_count);
    std::vector<std::vector<std::vector<double>>> backbone_dihedrals (activeFrame_count);
    std::vector<frame> activeFramesCopy (activeFrame_count);
    std::vector<std::vector<std::vector<double>>> backbone_dihedralsCopy (activeFrame_count);

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

    //backbone dihedrals: participating atoms
    std::vector<std::vector<int>> backbone_dihedral_atoms (me.protein_count); 

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
    outfile.open((me.outfilename + ".gz").c_str(), std::ios_base::out);

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

        //define file to read
        gz::igzstream file(me.input_files[i].c_str());

        //boolean specifying current active block
        bool isTitleBlock, isTimestepBlock, isPositionBlock, isGenboxBlock = false;

        //content holder for the blocks
        std::string titleBlock, timestepBlock, positionBlock, genboxBlock;

        //counters
        unsigned int positionBlock_counter, activeFrame_counter = 0, genBox_counter = 0;

        //check if it is possible to read file
        if (!file)
        {
            std::cerr << "cannot open otuput file" << "\n" << std::flush;
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
                            //backbone dihedrals
                            if (firstPass)
                            {
                                //parse frame prefix
                                for (int j = 0; j < currentFrame.prefix.size(); j++)
                                {
                                    //15944 SOLV  HW2    52493    1.966742392    4.000388119    1.796537431
                                    currentFrame.residue_number.push_back(stoi(currentFrame.prefix[j].substr(0,5)));
                                    currentFrame.residue_name.push_back(currentFrame.prefix[j].substr(6,4));
                                    currentFrame.atom_type.push_back(currentFrame.prefix[j].substr(12,4));
                                    currentFrame.atom_number.push_back(stoi(currentFrame.prefix[j].substr(16,8)));
                                }

                                //get the backbone atoms
                                for (int j = 0; j < me.protein_count; j++)
                                {                   
                                    for (int k = me.proteins[j][0]-1; k < me.proteins[j][1]; k++)
                                    {         
                                        if (currentFrame.atom_type[k] == "C   ")
                                        {
                                            backbone_dihedral_atoms[j].push_back(k);
                                        }
                                        if (currentFrame.atom_type[k] == "CA  ")
                                        {
                                            backbone_dihedral_atoms[j].push_back(k);                                        
                                        }
                                        if (currentFrame.atom_type[k] == "N   ")
                                        {                               
                                            backbone_dihedral_atoms[j].push_back(k);
                                        }
                                    }
                                }
                            }

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
                                        myThreads[g] = std::thread([&activeFrames, &me, &backbone_dihedrals, g, &backbone_dihedral_atoms]() {
                                            for (int f = (g * me.num_frame_per_thread); f < ((g + 1) * me.num_frame_per_thread); f++)
                                            {
                                                //DO SOMETHING
                                                CalculateBackboneDihedrals(&activeFrames[f], f, &me, &backbone_dihedrals, &backbone_dihedral_atoms);
                                                //CalculateBackboneDihedralsaAtan2(&activeFrames[f], f, &me, &backbone_dihedrals, &backbone_dihedral_atoms);
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
                                    backbone_dihedralsCopy = backbone_dihedrals;
                                    outfileThread = std::thread([&activeFramesCopy, &backbone_dihedralsCopy, &me, &outfile]()
                                    {
                                        for (int g = 0; g < activeFramesCopy.size(); g++)
                                        {
                                            //DO SOMETHING
                                            PrintDihedrals(&activeFramesCopy[g], g, &backbone_dihedralsCopy, outfile);
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
                    workers.push_back(std::thread([&activeFrames, &me, &backbone_dihedrals, g, &backbone_dihedral_atoms]() {
                        for (int f = (g * me.num_frame_per_thread); f < ((g + 1) * me.num_frame_per_thread); f++)
                        {
                            //DO SOMETHING
                            CalculateBackboneDihedrals(&activeFrames[f], f, &me, &backbone_dihedrals, &backbone_dihedral_atoms);
                            //CalculateBackboneDihedralsaAtan2(&activeFrames[f], f, &me, &backbone_dihedrals, &backbone_dihedral_atoms);
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
                    workers.push_back(std::thread([&activeFrames, &me, &backbone_dihedrals, g, &backbone_dihedral_atoms]() {
                        for (int f = (g * me.num_frame_per_thread); f < ((g + 1) * me.num_frame_per_thread); f++)
                        {
                            //DO SOMETHING
                            CalculateBackboneDihedrals(&activeFrames[f], f, &me, &backbone_dihedrals, &backbone_dihedral_atoms);
                            //CalculateBackboneDihedralsaAtan2(&activeFrames[f], f, &me, &backbone_dihedrals, &backbone_dihedral_atoms);
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
                    PrintDihedrals(&activeFrames[g], g, &backbone_dihedrals, outfile);
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