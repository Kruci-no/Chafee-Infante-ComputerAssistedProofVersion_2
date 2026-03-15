#include "chafeeInfanteVecFieldC1.h"
#include "../GallerkinProjections/gallerkinProjections.h"
#include <iostream>
using namespace std;
using namespace Algebra;
using namespace capd;


SeriesVector linear(int size,ParamsMap params){
    SeriesVector L(2);
    Series xx(1,-2,SeriesType::sin);
    Series ones(1,0,SeriesType::sin);
    xx = xx.resize(size);
    ones = ones.resize(size);
    L[0]= xx *interval(-1) + ones*params["lambda"]; 
    L[1]= xx *interval(-1) + ones*params["lambda"]; 
    return L;
}
SeriesVector linearSym(int size, ParamsMap params){
    SeriesVector L(2);
    Series xx(1,-2,SeriesType::sin_odd);
    Series ones(1,0,SeriesType::sin_odd);
    xx = xx.resize(size);
    ones = ones.resize(size);
    L[0]= xx *interval(-1) + ones*params["lambda"]; 
    return L;
}


SeriesVector nonlinearity(interval time, SeriesVector x,ParamsMap params){
    SeriesVector y(2);
    Series u_squere = squere(x[0]).resize(x[0].mainSize + 1);
    //y[0].print();
    //u_squere.print();
    Series uTo3Power = mult(u_squere,x[0]);
    y[0]= interval(-1)*(params["A"]*sin(time * params["omega"]) + params["B"]) * uTo3Power  ;
    y[1]= interval(-3)*(params["A"]*sin(time * params["omega"]) + params["B"]) * mult(u_squere,x[1])  ;
    return y;
}
IVector rest(interval time,SeriesVector x,Indexer indexer,ParamsMap params){
    auto splited_x = indexer.splitVector(x);
    SeriesVector xMain = splited_x.first ;
    SeriesVector xDiss = splited_x.second;
    SeriesVector y(2);
    Series uSquereDiss = mult(xDiss[0],x[0] + xMain[0]).resize(x[0].mainSize + 1 );
    Series uSquereMain = mult(xMain[0],xMain[0]).resize(x[0].mainSize + 1 );
    y[0] = mult(uSquereDiss,x[0]) + mult(uSquereMain,xDiss[0]);
    y[0] = interval(-1)*(params["A"]*sin(time * params["omega"]) + params["B"])*y[0] ;
    y[1] = mult(uSquereDiss,x[1]) + mult(uSquereMain,xDiss[1]);
    y[1] = interval(-3)*(params["A"]*sin(time * params["omega"]) + params["B"])*y[1] ;
    return indexer.getIVector(y);
}
Indexer makeIndexer(int size){
    std::vector<std::pair<int,int> > pairs(2*size);
    for(int i=0;i<size;i++ ){
        pairs[i].first = 0;
        pairs[i].second = i;
        pairs[i+size].first = 1;
        pairs[i+size].second = i;
    }
    Indexer indexer;
    indexer.pairs = pairs;
    return indexer;
}
VectorField getVectorFieldC1(int size,int serisesSize,ParamsMap params){
    VectorField vectorField;
    vectorField.indexer = makeIndexer(size);
    vectorField.params = params;
    vectorField.F_nonLinearity = nonlinearity;
    vectorField.F_rest = rest;
    vectorField.L = linear(serisesSize,params);
    vectorField.nMain = vectorField.indexer.size();
    vectorField.numOfNormalParams = 4;
    vectorField.f = IMap(gallerkinProjectionVecField(size,true));
    vectorField.f.setParameter("lambda",params["lambda"]);
    vectorField.f.setParameter("omega",params["omega"]);
    vectorField.f.setParameter("A",params["A"]);
    vectorField.f.setParameter("B",params["B"]);

    return vectorField;

}
/*
VectorField getVectorFieldSymC1(int size,int serisesSize,ParamsMap params){
    VectorField vectorField;
    vectorField.indexer = makeIndexer(size);
    vectorField.params = params;
    vectorField.F_nonLinearity = nonlinearity;
    vectorField.F_rest = rest;
    vectorField.L = linearSym(serisesSize,params);
    vectorField.nMain = vectorField.indexer.size();
    vectorField.numOfNormalParams = 4;
    vectorField.f = IMap(gallerkinProjectionVecFieldSym(size,true));
    vectorField.f.setParameter("lambda",params["lambda"]);
    vectorField.f.setParameter("omega",params["omega"]);
    vectorField.f.setParameter("A",params["A"]);
    vectorField.f.setParameter("B",params["B"]);
    return vectorField;
}
*/