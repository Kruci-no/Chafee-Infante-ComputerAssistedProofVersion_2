#include "InOut/inOut.h"
#include "capd/capdlib.h"
#include "ChafeeInfanteVecField/chafeeInfanteVecField.h"
#include "../DissipativePDE/Algebra/algebra.h"
#include "../DissipativePDE/Set/set.h"
#include "../DissipativePDE/Set/set.h"
#include "../DissipativePDE/SolverPDE/solverPDE.h"
#include<iostream>
using namespace std;
using namespace Algebra;
using namespace capd;
struct InfSeriesMatrix{
    interval s;
    Series rest;
    std::vector<Series> mainColumns;
    InfSeriesMatrix(){}
    InfSeriesMatrix(std::vector<Series> mainColumns,Series rest,interval s){
        this->mainColumns = mainColumns;
        this->rest = rest;
        this->s = s;
    }
    Series GetColumnSeries(int m){
        if(m <= 0)
            throw std::runtime_error("index Error\n");
        if(m <= mainColumns.size()){
            return mainColumns[m-1];
        }
        else 
        return rest * (interval(1.0)/interval(m)^s);
    }
    Series sumRest(){
        if(s <=interval(1)){
            throw std::runtime_error("s to small");
        }
        int mainColumnsNum = mainColumns.size();
        Series sumRest;
        sumRest = rest * ((interval(mainColumnsNum)^(-this->s + interval(1)))/(this->s - interval(1)));
        return sumRest;
    }
    interval sumMainAbs(){
        interval sum(0);
        for(int i=0;i<mainColumns.size();i++){
            sum += mainColumns[i].sumAbs();
        }
        return sum;
    }

    interval sumAllAbs(){
        return this->sumMainAbs()+ (this->sumRest()).sumAbs();
    }
    
};
bool cheakInclusion
    (Series u,int mainSize, int fullSize, DVector paramsDVector){
    ParamsMap params = {{"lambda",paramsDVector[0]},{"omega",2*Interval::pi()},{"A",paramsDVector[2]},{"B",paramsDVector[3]}};
    SeriesVector vec(1);
    vec[0] = u;
    //cout <<"Starting data:\n";
    //u.print();
    vec[0] = vec[0].resize(fullSize);
    IVector ivectorforRect2Set(mainSize);
    for(int i=0;i<mainSize;i++){
        ivectorforRect2Set[i] = vec[0].main[i];
    }

    InclRect2Set mainSet(ivectorforRect2Set);
    VectorField vectorFieldChaInf = getVectorFieldSym(mainSize,fullSize,params);
    Encloser encloser;
    Mover mover(vectorFieldChaInf,encloser);
    mover.setStep(1./1024);
    Set set(vec,mainSet);
    for(int i=0;i<1024;i++){
        mover.move(set,vectorFieldChaInf,true);
    }
    return set.vector[0].subset(vec[0]);

}
Series C1Computation(Series u,Series uH,int mainSize, int fullSize ,DVector paramsDVector){
    ParamsMap params = {{"lambda",paramsDVector[0]},{"omega",2*Interval::pi()},{"A",paramsDVector[2]},{"B",paramsDVector[3]}};
 
    SeriesVector vec(2);
    vec[0] = u;
    vec[1] = uH; 
    vec[0] = vec[0].resize(fullSize);
    vec[1] = vec[1].resize(fullSize);
    VectorField vectorFieldChaInfC1 = getVectorFieldC1(mainSize,fullSize,params);
    Encloser encloser;
    Mover mover(vectorFieldChaInfC1,encloser);
    mover.setStep(1./1024);
    IVector allVariables(2*mainSize);
     for(int i =0 ;i<mainSize;i++){
        allVariables[i]= vec[0].main[i];
        allVariables[i+mainSize]= vec[1].main[i];
    }
    InclRect2Set mainSet(allVariables);
    Set set(vec,mainSet);
    for(int i=0;i<1024;i++){
            mover.move(set,vectorFieldChaInfC1,true);
        }
    return set.vector[1];
}
Series TranformToSinSeries(Series sin_oddSeries){
    if(sin_oddSeries.type != SeriesType::sin_odd)
        throw std::runtime_error("Wrong Series Type");
    IVector newMain(2*sin_oddSeries.main.dimension()-1);
    for(int i=0;i<sin_oddSeries.main.dimension();i++){
        newMain[2*i] = sin_oddSeries.main[i];
    }
    sin_oddSeries.main = newMain;
    sin_oddSeries.type = SeriesType::sin;
    return Series(newMain,sin_oddSeries.C,sin_oddSeries.s,SeriesType::sin);
}
int main(){  
    InOut inOut;
	DVector  paramsDVector;
	inOut.paramsFile >> paramsDVector;
    ParamsMap params = {{"lambda",paramsDVector[0]},{"omega",2*Interval::pi()},{"A",paramsDVector[2]},{"B",paramsDVector[3]}};
    cout << "lambda = " << paramsDVector[0] <<endl ;
    cout << "A = "<< paramsDVector[2] <<endl ;
    cout << "B = "<< paramsDVector[3] <<endl;
    DVector DmainInitialData;
    inOut.initialValueFile >> DmainInitialData;
    IVector ImainInitialData(DmainInitialData);
    int mainC0Size;
    int fullC0Size;
    
    double eps;
    inOut.proofOptionsFile >> eps;
    for(int i=0;i<ImainInitialData.dimension();i++)
       ImainInitialData[i] += eps*interval(-1,1);
    double C_double;
    inOut.proofOptionsFile >> C_double ;
    interval C = interval(-1,1)*C_double;
    double s_double;
    inOut.proofOptionsFile >> s_double ;
    interval s = interval(s_double);
    Series startingSet(ImainInitialData,C,s,SeriesType::sin_odd);
    inOut.proofOptionsFile >> mainC0Size ;
    inOut.proofOptionsFile >> fullC0Size ;
    startingSet = startingSet.resize(fullC0Size);
    if(cheakInclusion(startingSet,mainC0Size, fullC0Size,paramsDVector)){
        cout << "Periodic Orbit found in Starting set\n";
        int mainC1Size ;
        int fullC1Size;
        int infiniteMatrixColumsNum;
        inOut.proofOptionsFile >> mainC1Size ;
        inOut.proofOptionsFile >> fullC1Size ;        
        inOut.proofOptionsFile >> infiniteMatrixColumsNum ;
        Series u = TranformToSinSeries(startingSet);
        Series h;
  
        InfSeriesMatrix infMatrix;
        infMatrix.mainColumns = std::vector<Series>(infiniteMatrixColumsNum);
        for(int i=0;i<infiniteMatrixColumsNum;i++){
            h= Series(IVector(infiniteMatrixColumsNum),interval(0),u.s,SeriesType::sin);
            h.main[i] = interval(1);
            infMatrix.mainColumns[i] = C1Computation(u,h,mainC1Size,fullC1Size,paramsDVector);
            cout << i << " compputed column\n";
        }
        interval sH = interval(-0);
        h= Series(IVector(infiniteMatrixColumsNum),interval(-1,1),sH,SeriesType::sin);
        infMatrix.rest = C1Computation(u,h,mainC1Size,fullC1Size,paramsDVector);
        cout << "Rest computed";
        for(int i=0;i<infiniteMatrixColumsNum;i++){
            cout <<"Column: " <<i <<"\n" ;
            infMatrix.mainColumns[i].print(); 
        }
        infMatrix.s = -sH;
        cout <<"Rest: s=" << infMatrix.s <<"\n" ;
        infMatrix.rest.print();
        interval normEstimations = 2*(infMatrix.sumMainAbs()+infMatrix.rest.sumAbs());
        cout <<"Norm estimation:\n" << normEstimations<<endl;
        if(normEstimations< interval(1)){
            cout << "Orbit is locally Attracting"<<endl;
        }
        else{
            cout << "Could not validate attracting"<<endl;
        }
    }
    else{
        cout << "could not find periodic Orbit\n";
    }


}
