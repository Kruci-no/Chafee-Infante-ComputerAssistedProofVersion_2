#ifndef LogisticModelTools_H
#define LogisticModelTools_H
#include <iostream>
#include <vector>
#include "capd/capdlib.h"
#include "../../DissipativePDE/VectorFieldMaker/vectorFieldMaker.h"

std::string gallerkinProjectionVecField(int size,bool addInclusionTerms) ;
std::string gallerkinProjectionVecFieldSym(int size,bool addInclusionTerms);
std::string gallerkinProjectionVecFieldC1(int size,bool addInclusionTerms) ;

#endif
