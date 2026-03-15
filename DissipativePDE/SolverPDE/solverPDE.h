#ifndef __SOLVERPDE_H__
#define __SOLVERPDE_H__
#include "../Algebra/algebra.h"
#include "../Set/set.h"
#include "../VectorField/vectorField.h"
#include <vector>

struct Encloser{
    Algebra::SeriesVector enclosureExtent; 
    Algebra::SeriesVector enclosurePointWiseExtent; 
    capd::interval validatedTimeStep;
    void enclose(
                Set& x,
                VectorField& vectorField,
                capd::interval dt,
                int refineNum,
                bool constStep, 
                bool comPointWiseEnclose);
};

struct Mover{
    capd::interval step;
    Encloser encloser;
    capd::IMap perturbMap;
	capd::IMultiMap* multiMap;
	capd::CWDiffInclSolver* solver;
    capd::interval maxDecay;
    Mover(VectorField& vectorField,Encloser encloser);
    void setStep(capd::interval step);
    void move(Set& x,VectorField& vectorField,bool constStep);
    void static perturb(capd::autodiff::Node t, capd::autodiff::Node in[], int dimIn, 
                        capd::autodiff::Node out[], int dimOut, capd::autodiff::Node params[], int noParams);
};
#endif