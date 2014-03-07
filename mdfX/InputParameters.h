#pragma once
#include "gzstream.h"
#include "Structs.h"
#include <string>
#include <sstream>
#include <vector>
#include <Eigen/Dense>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

class InputParameters
{
public:
	InputParameters(std::string inputFilePath);
	~InputParameters(void);
	void ParseInputFile(Structs::GenericParameters *this_params, std::string job_id);
	//void ParseInputFile(Structs::GenericParameters *this_params, std::string job_id, std::string param_file);
	void PrintGenericParameters(Structs::GenericParameters &me);
	void PrintGenericParameters(Structs::GenericParameters &me, gz::ogzstream &outfile);
	void PrintGenericParameters(Structs::GenericParameters &me, std::stringstream &outfile);

private:
	char* XPathGetText(std::string xpath_query, xmlXPathContextPtr xpathCtx);
	std::string inputFilePath;
};

