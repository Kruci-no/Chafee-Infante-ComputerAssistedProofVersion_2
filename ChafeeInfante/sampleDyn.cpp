#include "InOut/inOut.h"
#include "capd/capdlib.h"
#include "GallerkinProjections/gallerkinProjections.h"
#include"../Utils/SampleDyn/sampleDyn.h"
#include<iostream>

using namespace std;
using namespace capd;
int main(){
    InOut inOut;
    DVector  params;
	inOut.paramsFile >> params;
    //cout << Params;
    DVector initialValue;
    inOut.initialValueFile >> initialValue;
    cout << initialValue;
    int dim = initialValue.dimension();
    //cout << gallerkinProjectionVecFieldSym(dim,false);
    DMap chafeeInfante(gallerkinProjectionVecFieldSym(dim,false));
    chafeeInfante.setParameters(params);
    double startTime,durationTime;
    int stepNum;
    inOut.sampleDynOptionsFile >> startTime;
    inOut.sampleDynOptionsFile >> durationTime;
    inOut.sampleDynOptionsFile >> stepNum;
    inOut.sampleDynOutPutFile <<SampleDynamics(initialValue, params, chafeeInfante,startTime,durationTime,stepNum);

}