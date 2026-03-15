using namespace std;
#include "inOut.h"

void InOut::contruct()
{
  string direction = "ChafeeInfante/textFiles/";
 // string folderName = "textFiles/";
  this->paramsFile.open(direction +"params.txt");
  this->sampleDynOptionsFile.open(direction +"sampleDynOptions.txt") ;
  this->initialValueFile.open(direction +"initialValue.txt") ;
  this->sampleDynOutPutFile.open(direction +"sampleDynOutPut.txt");
  this->proofOptionsFile.open(direction +"proofOptions.txt");
}
InOut::InOut()
{
  contruct();
}
