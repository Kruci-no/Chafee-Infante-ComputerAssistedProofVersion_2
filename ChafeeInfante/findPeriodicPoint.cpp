#include "InOut/inOut.h"
#include "capd/capdlib.h"
#include "GallerkinProjections/gallerkinProjections.h"
#include<iostream>

using namespace std;
using namespace capd;

DVector findFixedPoint(int dim,DVector initialData, DVector params ,DMap f, double eps){
    
    f.setParameters(params);
    DOdeSolver solver(f,22);
    DTimeMap timeMap(solver);
    DVector prev =  initialData ;
    DVector next;
    int MaxIter = 1000;
    int iter = 0; 
    double err = eps+1;
    do{
        double startTime = 0; double endTime = 1;
        iter = iter +1 ;
        next = prev;
        next = timeMap(endTime,next,startTime);
		err = (next-prev).euclNorm();
        prev = next;
    }while( eps<err  &&iter < MaxIter);
    if(iter == MaxIter){
        throw std::runtime_error("Cannot Find fixed Point");
    }
    double startTime = 0; double endTime = 1;
    DMatrix initMatrix = DMatrix::Identity(dim) ;
    DMatrix monodromyMatrix(dim,dim);
    cout << "finded Fixed Point"<< next <<"\n";
    next = timeMap(endTime,next,initMatrix,monodromyMatrix,startTime);
    cout <<"Monodromi Matrix:" << monodromyMatrix << "\n";
    DMaxNorm maxNorm;
	cout <<"maxNorm of monodromyMatrix: " <<maxNorm(monodromyMatrix)<<"\n";
    cout << "finded Iterred Fixed Point"<< next <<"\n";
    DVector rV(dim), iV(dim);         // real and imaginary parts of eigenvalues
    DMatrix rVec(dim,dim), iVec(dim,dim); // real and imaginary parts of eigenvectors
    capd::alglib::computeEigenvaluesAndEigenvectors(monodromyMatrix,rV,iV,rVec,iVec);
    
    std::cout << "\n======================="
                << "\nmatrix A  = \n " << monodromyMatrix 
                << "\neigenValues = " << alglib::eigenvaluesToString(rV, iV, ", ");
    std::cout << "\neigenVectors : " << alglib::eigenvectorsToString(rVec, iVec);
    std::cout << "\neigenVectors (i-th column contains vector corresponing to i-th eigenvalue)= \n"
                << rVec << "+ i* " << iVec;
    return next;
}
DVector findFixedPointNewtonMethod(int dim,DVector initialData, DVector params ,DMap f, double eps){
    f.setParameters(params);
    DOdeSolver solver(f,22);
    DTimeMap timeMap(solver);
    DVector prev =  initialData ;
    DVector next;
    int MaxIter = 1000;
    int iter = 0; 
    double err = eps+1;
    do{
        double startTime = 0; double endTime = 1;
        iter = iter +1 ;
        next = prev;
        DMatrix initMatrix = DMatrix::Identity(dim) ;
        DMatrix monodromyMatrix(dim,dim);
        DVector f_value = timeMap(endTime,next,initMatrix,monodromyMatrix,startTime);
        DMatrix I_MinusDf = DMatrix::Identity(dim) - monodromyMatrix;
       // cout << I_MinusDf<<endl;
       // cout <<next<<endl;
        next = next + matrixAlgorithms::gauss(I_MinusDf,next-f_value);
        err = (next-prev).euclNorm();
        prev = next;
    }while( eps<err  &&iter < MaxIter);
    if(iter == MaxIter){
        throw std::runtime_error("Cannot Find fixed Point");
    }
    double startTime = 0; double endTime = 1;
    DMatrix initMatrix = DMatrix::Identity(dim) ;
    DMatrix monodromyMatrix(dim,dim);
    cout << "finded Fixed Point"<< next <<"\n";
    next = timeMap(endTime,next,initMatrix,monodromyMatrix,startTime);
    cout <<"Monodromi Matrix:" << monodromyMatrix << "\n";
    DMaxNorm maxNorm;
	cout <<"maxNorm of monodromyMatrix: " <<maxNorm(monodromyMatrix)<<"\n";
    cout << "finded Itered Fixed Point"<< next <<"\n";
    return next;

}
int main(){
    InOut inOut;
    DVector  params;
	inOut.paramsFile >> params;
    params[1] = 6.28318530718;
    //cout << Params;
    DVector initialValue;
    inOut.initialValueFile >> initialValue;
    int dim = initialValue.dimension();
    cout << "dimestion"<<  dim;
    DMap chafeeInfante(gallerkinProjectionVecField(dim,false));
    DMap chafeeInfanteSym(gallerkinProjectionVecFieldSym(dim,false));
    //cout << gallerkinProjectionVecField(dim,true) <<"\n";
   // cout << gallerkinProjectionVecFieldSym(dim,true) <<"\n";
    cout <<"Normal: ";
    findFixedPoint(dim,initialValue,params ,chafeeInfante,10e-15) ;
    cout <<"\n";
    cout <<"Symetric: " ;
    findFixedPoint(dim,initialValue,params ,chafeeInfanteSym,10e-15);
    cout<< "\n" ;
    cout <<"Newton: ";
    findFixedPointNewtonMethod(dim,initialValue,params ,chafeeInfante,10e-15) ;
    cout<<"\n";
    cout <<"NewtonSymetric: " ;
    findFixedPointNewtonMethod(dim,initialValue,params ,chafeeInfanteSym,10e-15);
    cout<< "\n" ;
}
/*
int main(){
    cout<< gallerkinProjectionVecFieldC1(3,true);

}*/