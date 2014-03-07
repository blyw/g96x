#define DEBUG
#define EIGEN_NO_DEBUG  
#define EIGEN_MPL2_ONLY
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO


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
//void InputParameters::ParseInputFile(Structs::GenericParameters *this_params, std::string job_id, std::string param_file)
void InputParameters::ParseInputFile(Structs::GenericParameters *this_params, std::string job_id)
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
	if (xpathCtx == NULL) {
		std::cout << "Unable to create new XPath context " << inputFilePath << std::endl;
		xmlFreeDoc(doc);
		exit(-1);
	}

	//evaluate xpath expression
	xpath = "/me/job[@id=\"" + job_id + "\"]";
	xpathObj = xmlXPathEvalExpression(BAD_CAST xpath.c_str(), xpathCtx);
	if (xpathObj == NULL) {
		std::cout << "Error: unable to evaluate xpath expression\n  " << xpath << std::endl;
		xmlXPathFreeContext(xpathCtx);
		xmlFreeDoc(doc);
		exit(-1);
	}

	for (int i = 0; i < xpathObj->nodesetval->nodeNr; i++)
	{
		//read current node;
		std::cout << xpathObj->nodesetval->nodeTab[i]->name << " id " << xmlGetProp(xpathObj->nodesetval->nodeTab[i], BAD_CAST((std::string)"id").c_str()) << std::endl;

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
		this_params->solvent_molecules(0, 0) = atoi(XPathGetText("./topology/solvent/@first_atom", xpathCtx));
		//get last atom
		this_params->solvent_molecules(1, 0) = atoi(XPathGetText("./topology/solvent/@last_atom", xpathCtx));
		//skip or keep solvent
		this_params->solvent_skip = atoi(XPathGetText("./topology/solvent/@skip", xpathCtx));
		//get number of atoms
		this_params->solvent_molecules(2, 0) = atoi(XPathGetText("./topology/solvent/number_of_atoms", xpathCtx));
		//get images factor
		this_params->solvent_molecules(3, 0) = atoi(XPathGetText("./topology/solvent/dimension_search", xpathCtx));

#ifdef DEBUG
		std::cout << "got solvent parameters" << std::endl;
#endif // DEBUG

		//get solute (cog) definitions
		xpath = "./topology/solutes_cog/solute";
		xmlXPathObjectPtr COGSolute = xmlXPathEvalExpression(BAD_CAST xpath.c_str(), xpathCtx);
		//loops through solvents
		//std::cout << xmlNodeGetContent(xpathObj->nodesetval->nodeTab[0]) << std::endl;
		//get all COG solutes
		this_params->solute_cog_molecules.resize(6, COGSolute->nodesetval->nodeNr);
		for (int j = 0; j < COGSolute->nodesetval->nodeNr; j++)
		{
			xpathCtx->node = COGSolute->nodesetval->nodeTab[j];
			//get first atom of COG solute
			this_params->solute_cog_molecules(0, j) = atoi(XPathGetText("./@first_atom", xpathCtx));
			//get last atom of COG solute
			this_params->solute_cog_molecules(1, j) = atoi(XPathGetText("./@last_atom", xpathCtx));
			//get number of atom of COG solute
			this_params->solute_cog_molecules(2, j) = atoi(XPathGetText("./number_of_atoms", xpathCtx));
			//get images factor
			this_params->solute_cog_molecules(3, j) = atoi(XPathGetText("./dimension_search", xpathCtx));
			//gather this molecule with respect tot this atom
			this_params->solute_cog_molecules(4, j) = atoi(XPathGetText("./@init_atom", xpathCtx));
			//skip or keep COG solute
			this_params->solute_cog_molecules(5, j) = std::strcmp(XPathGetText("./@skip", xpathCtx), "true") == 0;
		}
		xmlXPathFreeObject(COGSolute);

#ifdef DEBUG
		std::cout << "got cog_solutes parameters" << std::endl;
#endif // DEBUG

		xpathCtx->node = xpathObj->nodesetval->nodeTab[i];
		//get solute (cog) definitions
		xpath = "./topology/ions/ion";
		xmlXPathObjectPtr Ions = xmlXPathEvalExpression(BAD_CAST xpath.c_str(), xpathCtx);
		//loops through solvents
		//std::cout << xmlNodeGetContent(xpathObj->nodesetval->nodeTab[0]) << std::endl;
		//get all ions
		this_params->ion_molecules.resize(6, Ions->nodesetval->nodeNr);
		//std::cout << Ions->nodesetval->nodeNr << std::endl;
		for (int j = 0; j < Ions->nodesetval->nodeNr; j++)
		{
			std::cout << j << std::endl;
			xpathCtx->node = Ions->nodesetval->nodeTab[j];
			//get first atom of COG solute
			this_params->ion_molecules(0, j) = atoi(XPathGetText("./@first_atom", xpathCtx));
			//get last atom of COG solute
			this_params->ion_molecules(1, j) = atoi(XPathGetText("./@last_atom", xpathCtx));
			//get number of atom of COG solute
			this_params->ion_molecules(2, j) = atoi(XPathGetText("./number_of_atoms", xpathCtx));
			//get images factor
			this_params->ion_molecules(3, j) = atoi(XPathGetText("./dimension_search", xpathCtx));
			//gather this molecule with respect tot this atom
			this_params->ion_molecules(4, j) = atoi(XPathGetText("./@init_atom", xpathCtx));
			//skip or keep COG solute
			this_params->ion_molecules(5, j) = std::strcmp(XPathGetText("./@skip", xpathCtx), "true") == 0;
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
		this_params->solute_molecules.resize(5, solutes->nodesetval->nodeNr);

		for (int j = 0; j < solutes->nodesetval->nodeNr; j++)
		{
			xpathCtx->node = solutes->nodesetval->nodeTab[j];
			//get first atom of solute
			this_params->solute_molecules(0, j) = atoi(XPathGetText("./@first_atom", xpathCtx));
			//get last atom of solute
			this_params->solute_molecules(1, j) = atoi(XPathGetText("./@last_atom", xpathCtx));
			//get images factor
			this_params->solute_molecules(2, j) = atoi(XPathGetText("./dimension_search", xpathCtx));
			//gather this molecule with respect tot this atom
			this_params->solute_molecules(3, j) = atoi(XPathGetText("./@init_atom", xpathCtx));
			//skip or keep solute
			this_params->solute_molecules(4, j) = std::strcmp(XPathGetText("./@skip", xpathCtx), "true") == 0;
		};
		this_params->solute_count = this_params->solute_molecules.cols();
		xmlXPathFreeObject(solutes);

#ifdef DEBUG
		std::cout << "got solutes parameters" << std::endl;
#endif // DEBUG

		xpathCtx->node = xpathObj->nodesetval->nodeTab[i];
		//get distance-matrix atom sets definitions
		xpath = "./analysis/distance_matrix/atoms_set";
		xmlXPathObjectPtr distanceMatrix = xmlXPathEvalExpression(BAD_CAST xpath.c_str(), xpathCtx);
		//get all atomsets
		this_params->distanceMatrixSets.resize(5, distanceMatrix->nodesetval->nodeNr);

		//<atoms_set setA_start = "" setA_end = "" setB_start = "" setB_end = "" cut_off = "1.4" write_matrix = "true" / >
		for (int j = 0; j < distanceMatrix->nodesetval->nodeNr; j++)
		{
			xpathCtx->node = distanceMatrix->nodesetval->nodeTab[j];
			//get first atom of set A
			this_params->distanceMatrixSets(0, j) = atoi(XPathGetText("./@setA_first", xpathCtx));
			//get last atom of set A
			this_params->distanceMatrixSets(1, j) = atoi(XPathGetText("./@setA_last", xpathCtx));
			//get first atom of set B
			this_params->distanceMatrixSets(2, j) = atoi(XPathGetText("./@setB_first", xpathCtx));
			//get last atom of set B
			this_params->distanceMatrixSets(3, j) = atoi(XPathGetText("./@setB_last", xpathCtx));
			//compute the number of elements expected of the calculation
			this_params->distanceMatrixSets(4, j) = 0;
			for (int elements_i = 0; elements_i < (this_params->distanceMatrixSets(1, j) - this_params->distanceMatrixSets(0, j)); elements_i++)
			{
				for (int elements_j = 0; elements_j < (this_params->distanceMatrixSets(3, j) - this_params->distanceMatrixSets(2, j) - i); elements_j++)
				{
					this_params->distanceMatrixSets(4, j) += 1;
				}
			}
			this_params->distanceMatrixCutOff.push_back(atof(XPathGetText("./@cut_off", xpathCtx)));
			this_params->distanceMatrixWriteRaw.push_back(std::strcmp(XPathGetText("./@write_matrix", xpathCtx), "true") == 0);
		};
		xmlXPathFreeObject(distanceMatrix);

#ifdef DEBUG
		std::cout << "got analysis parameters" << std::endl;
#endif // DEBUG

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
void InputParameters::PrintGenericParameters(Structs::GenericParameters &me) {
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
		std::cout << "                                    " << me.solute_cog_molecules(0, i) << " " << me.solute_cog_molecules(1, i) << " " << me.solute_cog_molecules(2, i) << " " << me.solute_cog_molecules(3, i) << " " << me.solute_cog_molecules(4, i) << " " << me.solute_cog_molecules(5, i) << std::endl;
	}
	std::cout << "  number of solutes               : " << me.solute_count << std::endl;
	std::cout << "  solutes atom numbering          : " << std::endl;
	for (int i = 0; i < me.solute_molecules.cols(); i++)
	{
		std::cout << "                                    " << me.solute_molecules(0, i) << " " << me.solute_molecules(1, i) << " " << me.solute_molecules(2, i) << " " << me.solute_molecules(3, i) << " " << me.solute_molecules(4, i) << std::endl;
	}
	std::cout << "  number of ions                  : " << me.ion_molecules.cols() << std::endl;
	std::cout << "  ions atom numbering             : " << std::endl;
	for (int i = 0; i < me.ion_molecules.cols(); i++)
	{
		std::cout << "                                    " << me.ion_molecules(0, i) << " " << me.ion_molecules(1, i) << " " << me.ion_molecules(2, i) << " " << me.ion_molecules(3, i) << " " << me.ion_molecules(4, i) << " " << me.ion_molecules(5, i) << std::endl;
	}
	std::cout << "  skip ions                       : " << me.ions_skip << std::endl;
	std::cout << "  solvent cog images              : " << me.solvent_molecules(3, 0) << std::endl;
	std::cout << "  solvent size                    : " << me.solvent_molecules(2, 0) << std::endl;
	std::cout << "  trc refernce file               : " << me.trc_reference << std::endl;
	std::cout << "  comment                         : set center of geometry in viewer to (0,0,0);" << std::endl;
	std::cout << "                                    vmd: molinfo 1 set center \"{0.0 0.0 0.0}\"" << std::endl;
}