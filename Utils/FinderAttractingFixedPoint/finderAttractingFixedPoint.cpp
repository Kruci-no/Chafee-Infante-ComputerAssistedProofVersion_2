#include "finderAttractingFixedPoint.h"
using namespace capd;
using namespace std;
DVector findFixedPoint(DVector initialData,
                        std::function<capd::DVector(capd::DVector)> f,double eps)
{
    DVector prev =  initialData ;
    DVector next;
    int MaxIter = 100000;
    int iter = 0; 
    double err = eps+1;
    do{
        iter = iter +1 ;
        next  = f(prev);
		err = (next-prev).euclNorm();

    }while(err<eps && iter < MaxIter);
    if(iter == MaxIter){
        throw std::runtime_error("Cannot Find fixed Point");
    }
    return next;

}