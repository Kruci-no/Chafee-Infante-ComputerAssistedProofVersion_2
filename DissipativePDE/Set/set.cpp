#include "set.h"
using namespace std;
using namespace capd;
using namespace Algebra;
std::pair<SeriesVector, SeriesVector > Indexer::splitVector(SeriesVector& x) 
{
    SeriesVector x1 = x*interval(0);
    SeriesVector x2 = x;
    for(int i=0;i<this->pairs.size();i++){
            int l = pairs[i].first;
            int r = pairs[i].second; 
            x1[l].main[r] = x2[l].main[r];
            x2[l].main[r] = interval(0);
    }
    return std::pair<SeriesVector, SeriesVector>(x1,x2); 
}

int Indexer::size() 
{
    return pairs.size();    
}

capd::IVector Indexer::getIVector(SeriesVector& x) 
{
    IVector result(size());
        for(int i=0;i<size();i++){
            int l = pairs[i].first;
            int r = pairs[i].second; 
            result[i] = x[l].main[r];
        }
        return result;
}

void Indexer::intersectRepresetations(capd::InclRect2Set& mainModes,SeriesVector& x) 
{
    IVector main = (IVector) mainModes;
    for(int i=0;i<size();i++){
        int l = pairs[i].first;
        int r = pairs[i].second; 
        if(!intersection(main[i], x[l].main[r], x[l].main[r]) ){
            x.print();
            cout<< main<<"\n";
            throw std::runtime_error("Set intersectRepresetations - empty intersection\n");
        }
    }
}

void Indexer::makeCosistend(capd::InclRect2Set& mainModes,SeriesVector& x) 
{
    IVector main = (IVector) mainModes;
    for(int i=0;i<size();i++){
        int l = pairs[i].first;
        int r = pairs[i].second; 
        x[l].main[r] = main[i];
    }
}

Set::Set(SeriesVector vector,capd::InclRect2Set mainModes) 
:vector(vector),mainModes(mainModes)
{
}

void Set::intersectRepresetations(Indexer indexer) 
{
    indexer.intersectRepresetations(mainModes,vector);    
}

void Set::makeCosistend(Indexer indexer) 
{
    indexer.makeCosistend(mainModes,vector);
}

interval Set::getCurrentTime() 
{
    return mainModes.getCurrentTime();
}