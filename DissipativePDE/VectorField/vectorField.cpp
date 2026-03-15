
#include "vectorField.h"

using namespace std;
using namespace capd;
using namespace Algebra;


SeriesVector VectorField::computeNonLinearity(interval time,SeriesVector x) 
{
    return F_nonLinearity(time,x,params);
}

capd::IVector VectorField::computeInclusion(interval time, SeriesVector x) 
{
    return F_rest(time,x,indexer,params);
}