#include "../DissipativePDE/Algebra/algebra.h"
#include "../DissipativePDE/Set/set.h"
#include "../DissipativePDE/VectorField/vectorField.h"
#include "../DissipativePDE/SolverPDE/solverPDE.h"
#include "../DissipativePDE/VectorFieldMaker/vectorFieldMaker.h"
#include "../LogisticModel/GallerkinProjections/gallerkinProjections.h"
#include <iostream>
//#include "vectorField.h"
//#include "solverPDE.h"
using namespace capd;
using namespace std;
using namespace Algebra;
//using namespace vectorFieldMaker;*/
int main(){
    auto var = gallerkinProjectionVecFieldSym(2,false);
    //auto var = IMap(gallerkinProjectionVecField(2,false));
    //cout << a;
    cout << var;
    IVector wx({interval(-1),interval(1),interval(2),interval(0),interval(0)});
    IVector wy({interval(1),interval(2),interval(1),interval(0)});
    Series x(wx,interval(-1,1),interval(4),SeriesType::cos_even );
    Series y(wy,interval(-1,1),interval(-2),SeriesType::sin_odd );
    x.print();
    y.print();
    //Series z = mult(x,y);
    //z.print();
    
}