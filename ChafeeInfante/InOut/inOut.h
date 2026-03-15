#ifndef InOut_H_
#define InOut_H_
#include <iostream>
#include <fstream>
#include <string>
struct InOut
{
     std::ifstream paramsFile;
     std::ifstream sampleDynOptionsFile;
     std::ifstream initialValueFile;
     std::ifstream proofOptionsFile;
     std::ofstream sampleDynOutPutFile;

     InOut();
     void contruct();
};
#endif //_EXAMPLES_PROJECTSTARTER_OUTPUT_H_
