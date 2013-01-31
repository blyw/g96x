#define DEBUG
#define POSIX

#include <signal.h>
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
    int num_thread;
    int num_frame_per_thread;
    int num_thread_real;
    int output_fragment_size;
    int output_fragment_skipframes;
    double output_fragment_skiptime;
    int solute_count;
    std::vector<std::vector<int>> solute_molecules;
    std::vector<std::vector<int>> solutes_cog_molecules;
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
};

//code checked 20130122: CORRECT!
//solute gathering i.e. cluster algorithm (nearest neighbour)
void NearestNeighbourFinder(frame *framedata, std::vector<std::string> *periodic_interactions, int FrameId, params *me){ 
    //init variables
    double distance_shortest;

    //output
    std::stringstream ss;

    //cut-off is derived from (longest bond)^2
    double cut_off = me->distance_cut_off;
    int molecule_start, molecule_end, molecule_start_previous, periodic_copies;

    //dimension expansion factor
    periodic_copies = 1;

    //get end and start of search
    for (int s = 0; s < me->solute_count; s++)
    {
        //define first atom, last atom and number of periodic copies of a solute molecule
        if (s == 0)
        {
            molecule_start = me->solute_molecules[s][0];
            periodic_copies = me->solute_molecules[s][2];
        }
        if (s == (me->solute_count - 1))
        {
            molecule_end = me->solute_molecules[s][1];
        }
        if (periodic_copies < me->solute_molecules[s][2])
        {
            periodic_copies = me->solute_molecules[s][2];
        }        
    }

    //calculated 1/2 distance matrix and pick out the interacting atoms
    for (int i = molecule_start - 1; i < molecule_end; i++)
    {
        distance_shortest = 1E20;
        for (int j = i; j < molecule_end; j++)
        {
            for (int x = 0-periodic_copies; x <= periodic_copies; x++)
            {
                for (int y = 0-periodic_copies; y <= periodic_copies; y++)
                {
                    for (int z = 0-periodic_copies; z <= periodic_copies; z++)
                    {
                        if (!(x==0 && y==0 && z==0))
                        {
                            double distance =
                                (((framedata->x[j] + x * framedata->box_length_x) - framedata->x[i]) * ((framedata->x[j] + x * framedata->box_length_x) - framedata->x[i])) +
                                (((framedata->y[j] + y * framedata->box_length_y) - framedata->y[i]) * ((framedata->y[j] + y * framedata->box_length_y) - framedata->y[i])) +
                                (((framedata->z[j] + z * framedata->box_length_z) - framedata->z[i]) * ((framedata->z[j] + z * framedata->box_length_z) - framedata->z[i]));
                            //std::cout << framedata->x[i] << "," << framedata->y[i] << "," << framedata->z[i] << " " << framedata->x[j] + x * framedata->box_length_x << "," << framedata->y[j] + y * framedata->box_length_y << "," << framedata->z[j]  + z * framedata->box_length_z << " " << sqrt(distance_shortest) << std::endl;
                            if (distance < distance_shortest)
                            {
                                distance_shortest = distance;
                            }
                        }
                    }
                }
            }

            //std::cout << " " << framedata->frame_id << " " << std::right << std::setw(10) << i  << " " << std::right << std::setw(10) << j << std::setw(50) << " [" + framedata->prefix[i] + " - " + framedata->prefix[j] + "] " << std::setw(16) << std::setprecision(9) << sqrt(distance_shortest) << std::endl;

            //check whether the atom see a periodoc copy
            if (distance_shortest <= cut_off)
            {
                ss << " " << framedata->frame_id << " " << std::right << std::setw(10) << i  << " " << std::right << std::setw(10) << j << std::setw(50) << " [" + framedata->prefix[i] + " - " + framedata->prefix[j] + "] " << std::setw(16) << std::setprecision(9) << sqrt(distance_shortest) << "\n";
            }
        }
    }

    //pass interaction list 
    if (ss.str().length() > 0)
    {        
        (*periodic_interactions)[FrameId] = ss.str();
        //std::cout << (*periodic_interactions)[FrameId] << std::endl;
        ss.clear();
    }
}

//code checked 20130122: CORRECT!
//write out the data in either CNF or PDB compatible format
void WriteOutPeriodicInteraction(frame *framedata, std::string &output, gz::ogzstream &outfile) {
    //check if it is possible to read file
    if (!outfile)
    {
        std::cerr << "cannot open otuput file" << "\n";
    }
    else {
        std::cout << "#frame_id " << framedata->frame_id << std::endl; 
        //std::cout << output << std::endl;  
        outfile << "#frame_id " << framedata->frame_id << "\n"; 
        outfile << output;
        outfile.flush();
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
        //get distance cut
        this_params->distance_cut_off = atof(XPathGetText("./variables/distance_cut_off", xpathCtx)) * atof(XPathGetText("./variables/distance_cut_off", xpathCtx));
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
    std::vector<std::string> periodic_interactions (activeFrame_count);
    std::vector<frame> activeFramesCopy (activeFrame_count);
    std::vector<std::string> periodic_interactionsCopy (activeFrame_count);

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
                                        myThreads[g] = std::thread([&activeFrames, &me, &periodic_interactions, g]() {
                                            for (int f = (g * me.num_frame_per_thread); f < ((g + 1) * me.num_frame_per_thread); f++)
                                            {
                                                NearestNeighbourFinder(&activeFrames[f], &periodic_interactions, f, &me);
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
                                    periodic_interactionsCopy = periodic_interactions;
                                    outfileThread = std::thread([&activeFramesCopy, &periodic_interactionsCopy, &me, &outfile]()
                                    {
                                        for (int g = 0; g < periodic_interactionsCopy.size(); g++)
                                        {
                                            WriteOutPeriodicInteraction(&activeFramesCopy[g], periodic_interactionsCopy[g], outfile);
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
                    workers.push_back(std::thread([&activeFrames, &me, &periodic_interactions, g]() {
                        for (int f = (g * me.num_frame_per_thread); f < ((g + 1) * me.num_frame_per_thread); f++)
                        {
                            NearestNeighbourFinder(&activeFrames[f], &periodic_interactions, f, &me);
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
                    workers.push_back(std::thread([&activeFrames, &me, &periodic_interactions, g]() {
                        for (int f = (g * me.num_frame_per_thread); f < ((g + 1) * me.num_frame_per_thread); f++)
                        {
                            NearestNeighbourFinder(&activeFrames[f], &periodic_interactions, f, &me);
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
                    WriteOutPeriodicInteraction(&activeFrames[g], periodic_interactions[g], outfile);
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