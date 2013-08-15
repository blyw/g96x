//#define DEBUG
#define EIGEN_NO_DEBUG  
#define EIGEN_MPL2_ONLY
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include "gzstream.h"
#include "InputParameters.h"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <thread>
#include <Eigen/Dense>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <Eigen/Dense>
#include "Structs.h"

InputParameters::InputParameters(std::string _inputFilePath)
{
    inputFilePath = _inputFilePath;
}

InputParameters::~InputParameters(void)
{
}

//get the contents of a node
char* InputParameters::XPathGetText(std::string xpath_query, xmlXPathContextPtr xpathCtx) 
{
    xmlXPathObjectPtr xpathObj = xmlXPathEvalExpression(BAD_CAST xpath_query.c_str(), xpathCtx);
    return (char*)xmlNodeGetContent(xpathObj->nodesetval->nodeTab[0]);	
}

//parse parameter input file
//void InputParameters::ParseInputFile(Structs::InputParametersFrameout *this_params, std::string job_id, std::string param_file)
void InputParameters::ParseInputFile(Structs::InputParametersFrameout *this_params, std::string job_id)
{
    this_params->ions_skip = false;
    this_params->cog_write = false;
    this_params->correction_translation = false;
    this_params->correction_rotation = false;

    //initialize variables for libxml2
    std::string xpath;
    xmlDocPtr doc;
    xmlXPathContextPtr xpathCtx; 
    xmlXPathObjectPtr xpathObj; 

    //load XML document
    doc = xmlParseFile(inputFilePath.c_str());
    //if the documuent did not parse properly
    if (doc == NULL)
    {
        std::cout << "Unable to parse " << inputFilePath << std::endl;
        exit(-1);
    }

    //create xpath context
    xpathCtx = xmlXPathNewContext(doc);
    if(xpathCtx == NULL) {
        std::cout << "Unable to create new XPath context " << inputFilePath << std::endl;
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
        //get images factor
        this_params->solvent_molecules(3,0) = atoi(XPathGetText("./topology/solvent/dimension_search", xpathCtx));

#ifdef DEBUG
        std::cout << "got solvent parameters" << std::endl;
#endif // DEBUG

        //get solute (cog) definitions
        xpath = "./topology/solutes_cog/solute";
        xmlXPathObjectPtr COGSolute = xmlXPathEvalExpression(BAD_CAST xpath.c_str(), xpathCtx);
        //loops through solvents
        //std::cout << xmlNodeGetContent(xpathObj->nodesetval->nodeTab[0]) << std::endl;
        //get all COG solutes
        this_params->solute_cog_molecules.resize(5, COGSolute->nodesetval->nodeNr); 
        for (int j = 0; j < COGSolute->nodesetval->nodeNr; j++)
        {
            xpathCtx->node = COGSolute->nodesetval->nodeTab[j];
            //get first atom of COG solute
            this_params->solute_cog_molecules(0,i) = atoi(XPathGetText("./@first_atom", xpathCtx));
            //get last atom of COG solute
            this_params->solute_cog_molecules(1,i) = atoi(XPathGetText("./@last_atom", xpathCtx));
            //get number of atom of COG solute
            this_params->solute_cog_molecules(2,i) = atoi(XPathGetText("./number_of_atoms", xpathCtx));
            //get images factor
            this_params->solute_cog_molecules(3,i) = atoi(XPathGetText("./dimension_search", xpathCtx));
            //gather this molecule with respect tot this atom
            this_params->solute_cog_molecules(4,i) = atoi(XPathGetText("./@init_atom", xpathCtx));
            //skip or keep COG solute
            this_params->solute_cog_molecules(5,i) = atoi(XPathGetText("./@skip", xpathCtx));
        }
        xmlXPathFreeObject(COGSolute);

#ifdef DEBUG
        std::cout << "got cog_solutes parameters" << std::endl;
#endif // DEBUG

        //get solute (cog) definitions
        xpath = "./topology/ions/ion";
        xmlXPathObjectPtr Ions = xmlXPathEvalExpression(BAD_CAST xpath.c_str(), xpathCtx);
        //loops through solvents
        //std::cout << xmlNodeGetContent(xpathObj->nodesetval->nodeTab[0]) << std::endl;
        //get all ions
        this_params->ion_molecules.resize(5, Ions->nodesetval->nodeNr); 
        for (int j = 0; j < Ions->nodesetval->nodeNr; j++)
        {
            std::cout << j << std::endl;
            xpathCtx->node = Ions->nodesetval->nodeTab[j];
            //get first atom of COG solute
            this_params->ion_molecules(0,i) = atoi(XPathGetText("./@first_atom", xpathCtx));
            //get last atom of COG solute
            this_params->ion_molecules(1,i) = atoi(XPathGetText("./@last_atom", xpathCtx));
            //get number of atom of COG solute
            this_params->ion_molecules(2,i) = atoi(XPathGetText("./number_of_atoms", xpathCtx));
            //get images factor
            this_params->ion_molecules(3,i) = atoi(XPathGetText("./dimension_search", xpathCtx));
            //gather this molecule with respect tot this atom
            this_params->ion_molecules(4,i) = atoi(XPathGetText("./@init_atom", xpathCtx));
            //skip or keep COG solute
            this_params->ion_molecules(5,i) = atoi(XPathGetText("./@skip", xpathCtx));
        }
        xmlXPathFreeObject(Ions);

#ifdef DEBUG
        std::cout << "got ions parameters" << std::endl;
#endif // DEBUG

        xpathCtx->node = xpathObj->nodesetval->nodeTab[i];
        //get solute definitions
        xpath = "./topology/solutes/solute";
        xmlXPathObjectPtr solutes = xmlXPathEvalExpression(BAD_CAST xpath.c_str(), xpathCtx);
        //get all solutes
        this_params->solute_molecules.resize(4, solutes->nodesetval->nodeNr);

        for (int j = 0; j < solutes->nodesetval->nodeNr; j++)
        {
            xpathCtx->node = solutes->nodesetval->nodeTab[j];
            //get first atom of solute
            this_params->solute_molecules(0,j) = atoi(XPathGetText("./@first_atom", xpathCtx));
            //get last atom of solute
            this_params->solute_molecules(1,j) = atoi(XPathGetText("./@last_atom", xpathCtx));
            //get images factor
            this_params->solute_molecules(2,j) = atoi(XPathGetText("./dimension_search", xpathCtx));
            //gather this molecule with respect tot this atom
            this_params->solute_molecules(3,j) = atoi(XPathGetText("./@init_atom", xpathCtx));
            //skip or keep solute
            this_params->solute_molecules(4,j) = atoi(XPathGetText("./@skip", xpathCtx));
        };
        this_params->solute_count = this_params->solute_molecules.cols();
        xmlXPathFreeObject(solutes);

#ifdef DEBUG
        std::cout << "got solutes parameters" << std::endl;
#endif // DEBUG

        //read the variables section
        xpathCtx->node = xpathObj->nodesetval->nodeTab[i];
        //to gather or not
        this_params->gather = atoi(XPathGetText("./analysis/frameout/@gather", xpathCtx));
        //write out cog or not
        this_params->cog_write = atoi(XPathGetText("./analysis/frameout/cog_write", xpathCtx));
        //correct for translation
        this_params->correction_translation = atoi(XPathGetText("./analysis/frameout/correction_translation", xpathCtx));
        //correct for rotation
        this_params->correction_rotation = atoi(XPathGetText("./analysis/frameout/correction_rotation", xpathCtx));
        //get distance cut
        this_params->distance_cut_off = atof(XPathGetText("./analysis/frameout/distance_cut_off", xpathCtx));

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

        //get images factor
        this_params->trc_reference = XPathGetText("./variables/reference_file", xpathCtx);
    }

    //clean up
    xmlXPathFreeObject(xpathObj);
    xmlXPathFreeContext(xpathCtx); 
    xmlFreeDoc(doc);
}

//print out the parameters gained from parsing the input file
void InputParameters::PrintInputParametersFrameOut(Structs::InputParametersFrameout &me) {
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
    std::cout << "  correct translation             : " << me.correction_translation << std::endl;
    std::cout << "  correct rotation                : " << me.correction_rotation << std::endl;
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
    for (int i = 0; i < me.solute_cog_molecules.cols(); i++)
    {
        std::cout << "                                    " << me.solute_cog_molecules(0,i) << " " << me.solute_cog_molecules(1,i) << " " << me.solute_cog_molecules(2,i) << " " << me.solute_cog_molecules(3,i) <<  " " << me.solute_cog_molecules(4,i) << std::endl;
    }
    std::cout << "  number of solutes               : " << me.solute_count << std::endl;
    std::cout << "  solutes atom numbering          : " << std::endl;
    for (int i = 0; i < me.solute_molecules.cols(); i++)
    {
        std::cout << "                                    " << me.solute_molecules(0,i) << " " << me.solute_molecules(1,i) << " " << me.solute_molecules(2,i) << " " << me.solute_molecules(3,i) << std::endl;
    }
    std::cout << "  ions cog                        : " << std::endl;
    for (int i = 0; i < me.ion_molecules.cols(); i++)
    {
        std::cout << "                                    " << me.ion_molecules(0,i) << " " << me.ion_molecules(1,i) << " " << me.ion_molecules(2,i) << " " << me.ion_molecules(3,i) <<  " " << me.ion_molecules(4,i) << std::endl;
    }
    std::cout << "  skip ions                       : " << me.ions_skip << std::endl;
    std::cout << "  solvent cog images              : " << me.solvent_molecules(3,0) << std::endl;
    std::cout << "  solvent size                    : " << me.solvent_molecules(2,0) << std::endl;
    std::cout << "  trc refernce file               : " << me.trc_reference << std::endl;  
    std::cout << "  gathering enabled               : " << me.gather << std::endl;  
}

//print out the parameters gained from parsing the input file
void InputParameters::PrintInputParametersFrameOut(Structs::InputParametersFrameout &me, gz::ogzstream &outfile) {
    outfile << "#frameout: input parameters defined by user" << std::endl;
    outfile << "#--------------------------------------------------------------------------------" << std::endl;
    outfile << std::left;
    outfile << "#   number of atom records          : " << me.atomrecords << std::endl;
    outfile << "#   distance cut-off                : " << me.distance_cut_off << std::endl;
    outfile << "#   input format                    : " << me.informat << std::endl;
    outfile << "#   input files                     : " << std::endl;
    for (unsigned int i = 0; i < me.input_files.size(); i++)
    {
        outfile << "#                                     " << me.input_files[i] << std::endl;
    }
    outfile << "#   write-out COG                   : " << me.cog_write << std::endl;
    outfile << "#   correct translation             : " << me.correction_translation << std::endl;
    outfile << "#   correct rotation                : " << me.correction_rotation << std::endl;
    outfile << "#   frames per thread               : " << me.num_frame_per_thread << std::endl;
    outfile << "#   number of threads               : " << me.num_thread << std::endl;
    outfile << "#   real number of threads          : " << me.num_thread_real << std::endl;
    outfile << "#   output filename prefix          : " << me.outfilename << std::endl;
    outfile << "#   output format                   : " << me.outformat << std::endl;
    outfile << "#   output frames per file          : " << me.output_fragment_size << std::endl;
    outfile << "#   output every n frame(s)         : " << me.output_fragment_skipframes << std::endl;
    outfile << "#   output every n picoseconds      : " << me.output_fragment_skiptime << std::endl;
    outfile << "#   initial reference coordinates   : " << me.ref_coords.x() << " " << me.ref_coords.y() << " " << me.ref_coords.z() << std::endl;
    outfile << "#   skip solvent                    : " << me.solvent_skip << std::endl;
    outfile << "#   solutes cog                     : " << std::endl;
    for (int i = 0; i < me.solute_cog_molecules.cols(); i++)
    {
        outfile << "#                                     " << me.solute_cog_molecules(0,i) << " " << me.solute_cog_molecules(1,i) << " " << me.solute_cog_molecules(2,i) << " " << me.solute_cog_molecules(3,i) <<  " " << me.solute_cog_molecules(4,i) << std::endl;
    }
    outfile << "#   number of solutes               : " << me.solute_count << std::endl;
    outfile << "#   solutes atom numbering          : " << std::endl;
    for (int i = 0; i < me.solute_molecules.cols(); i++)
    {
        outfile << "#                                     " << me.solute_molecules(0,i) << " " << me.solute_molecules(1,i) << " " << me.solute_molecules(2,i) << " " << me.solute_molecules(3,i) << std::endl;
    }
    outfile << "#   ions cog                        : " << std::endl;
    for (int i = 0; i < me.ion_molecules.cols(); i++)
    {
        outfile << "#                                     " << me.ion_molecules(0,i) << " " << me.ion_molecules(1,i) << " " << me.ion_molecules(2,i) << " " << me.ion_molecules(3,i) <<  " " << me.ion_molecules(4,i) << std::endl;
    }
    outfile << "#   skip ions                       : " << me.ions_skip << std::endl;
    outfile << "#   solvent cog images              : " << me.solvent_molecules(3,0) << std::endl;
    outfile << "#   solvent size                    : " << me.solvent_molecules(2,0) << std::endl;
    outfile << "#   trc refernce file (prefix)      : " << me.trc_reference << std::endl;  
    outfile << "#   gathering enabled               : " << me.gather << std::endl;  
}

//print out the parameters gained from parsing the input file
void InputParameters::PrintInputParametersFrameOut(Structs::InputParametersFrameout &me, std::stringstream &outfile) {
    outfile << "# input parameters defined" << std::endl;
    outfile << std::left;
    outfile << "#   number of atom records          : " << me.atomrecords << std::endl;
    outfile << "#   distance cut-off                : " << me.distance_cut_off << std::endl;
    outfile << "#   input format                    : " << me.informat << std::endl;
    outfile << "#   input files                     : " << std::endl;
    for (unsigned int i = 0; i < me.input_files.size(); i++)
    {
        outfile << "#                                     " << me.input_files[i] << std::endl;
    }
    outfile << "#   write-out COG                   : " << me.cog_write << std::endl;
    outfile << "#   correct translation             : " << me.correction_translation << std::endl;
    outfile << "#   correct rotation                : " << me.correction_rotation << std::endl;
    outfile << "#   frames per thread               : " << me.num_frame_per_thread << std::endl;
    outfile << "#   number of threads               : " << me.num_thread << std::endl;
    outfile << "#   real number of threads          : " << me.num_thread_real << std::endl;
    outfile << "#   output filename prefix          : " << me.outfilename << std::endl;
    outfile << "#   output format                   : " << me.outformat << std::endl;
    outfile << "#   output frames per file          : " << me.output_fragment_size << std::endl;
    outfile << "#   output every n frame(s)         : " << me.output_fragment_skipframes << std::endl;
    outfile << "#   output every n picoseconds      : " << me.output_fragment_skiptime << std::endl;
    outfile << "#   reference coordinates           : " << me.ref_coords.x() << " " << me.ref_coords.y() << " " << me.ref_coords.z() << std::endl;
    outfile << "#   skip solvent                    : " << me.solvent_skip << std::endl;
    outfile << "#   solutes cog                     : " << std::endl;
    for (int i = 0; i < me.solute_cog_molecules.cols(); i++)
    {
        outfile << "#                                     " << me.solute_cog_molecules(0,i) << " " << me.solute_cog_molecules(1,i) << " " << me.solute_cog_molecules(2,i) << " " << me.solute_cog_molecules(3,i) <<  " " << me.solute_cog_molecules(4,i) << std::endl;
    }
    outfile << "#   number of solutes               : " << me.solute_count << std::endl;
    outfile << "#   solutes atom numbering          : " << std::endl;
    for (int i = 0; i < me.solute_molecules.cols(); i++)
    {
        outfile << "#                                     " << me.solute_molecules(0,i) << " " << me.solute_molecules(1,i) << " " << me.solute_molecules(2,i) << " " << me.solute_molecules(3,i) << std::endl;
    }
    outfile << "#   ions cog                        : " << std::endl;
    for (int i = 0; i < me.ion_molecules.cols(); i++)
    {
        outfile << "#                                     " << me.ion_molecules(0,i) << " " << me.ion_molecules(1,i) << " " << me.ion_molecules(2,i) << " " << me.ion_molecules(3,i) <<  " " << me.ion_molecules(4,i) << std::endl;
    }
    outfile << "#   skip ions                       : " << me.ions_skip << std::endl;
    outfile << "#   solvent cog images              : " << me.solvent_molecules(3,0) << std::endl;
    outfile << "#   solvent size                    : " << me.solvent_molecules(2,0) << std::endl;
    outfile << "#   trc refernce file               : " << me.trc_reference << std::endl;  
    outfile << "#   gathering enabled               : " << me.gather << std::endl;  
}