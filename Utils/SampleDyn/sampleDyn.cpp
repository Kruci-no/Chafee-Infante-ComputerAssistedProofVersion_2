/*start time duration time stepNum */
# include <iostream>
# include "capd/capdlib.h"
using namespace capd;
using namespace std;
#include <vector>
/*start time duration time stepNum */

string SampleDynamics(DVector initialValue, DVector params, DMap map,double startTime, double durationTime,int stepNum){
  
  std::cout << "SampleDynamics Function:\n " ;
  std::cout << "Params " << params <<"\n";
  std::cout << "inital Time:" << startTime <<"\n";
  std::cout << "duration"<< durationTime<<"\n";
  std::cout << "dt= "<< durationTime/stepNum << '\n';
  int dim = initialValue.dimension() ;
  std::cout << "dimension "<< dim << '\n';
  std::cout << initialValue << '\n';
  
  DVector u = initialValue;
  double actualTime = startTime;
  double dt = durationTime/stepNum ;
  map.setParameters(params);
	DOdeSolver solver(map,22);
  solver.setStep(dt);
  solver.turnOffStepControl();
  string result = "{";
  for(int i = 0; i<=stepNum; i++){
    result=result+"{";
    u = solver(actualTime,u);
    for(int j=0;j<dim-1;j++){
       result=result+to_string(u[j])+",";
     }
     result=result+to_string(u[dim-1]);
     result = result +"},\n";
    //x = x +  timeMap((finalTime*1.0)/stepNum,initialValue) +",\n";
    //std::cout <<"At time:"<< finalTime * ((i * 1.0)/stepNum) <<" "
              //inOut.outPutFile << fixed <<timeMap((finalTime*1.0)/stepNum,initialValue)<< ","<< '\n';
  }
  result.pop_back();
  result.pop_back();
  result = result+"}";
  return result;
}

