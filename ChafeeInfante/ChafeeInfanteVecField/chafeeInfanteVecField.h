#include "../../DissipativePDE/Algebra/algebra.h"
#include "../../DissipativePDE/VectorField/vectorField.h"
#include "../../DissipativePDE/Set/set.h"
#include "capd/capdlib.h"
#ifndef ChafeeInfanteVecField_H
#define ChafeeInfanteVecField_H
VectorField getVectorField(int size,int serisesSize,ParamsMap params);
VectorField getVectorFieldSym(int size,int serisesSize,ParamsMap params);
VectorField getVectorFieldC1(int size,int serisesSize,ParamsMap params);

#endif