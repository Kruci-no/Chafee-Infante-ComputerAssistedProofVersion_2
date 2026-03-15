#ifndef __VECTORFIELD_H__
#define __VECTORFIELD_H__
#include "capd/capdlib.h"
#include "../Algebra/algebra.h"
#include "../Set/set.h"
#include <utility>
#include <vector>
#include <unordered_map>
typedef std::unordered_map<std::string,capd::interval> ParamsMap;
struct VectorField{
    // liczba glownych zmiennych
    int nMain;
    //liczba normalnych parametrow
    int numOfNormalParams;
    //Funkcja do Obliczenia nielinowosci w wyrazeniu
    std::function<Algebra::SeriesVector(capd::interval time, Algebra::SeriesVector&, ParamsMap)> F_nonLinearity;
    //Funkcja do obliczania pominetych czonow z projekcji Galerkina
    std::function<capd::IVector (capd::interval time, Algebra::SeriesVector&, Indexer , ParamsMap)> F_rest;
    //funkcja do liczenia pola wektorwego dla projekcji Galerkina
	capd::IMap f;
    //Obiekt parujacy indeksy zmiennych z projekcji galerkina ze zmiennymi w SeriesVector
    Indexer indexer;
    //Parametry w słowniku
    ParamsMap params;
    //Operator liniowy L reprezentowany przez rozbiezne szeregi 
    Algebra::SeriesVector L;
    //Obliczanie nieliniowosci kozystajac z F_nonLinearity
    Algebra::SeriesVector computeNonLinearity(capd::interval time, Algebra::SeriesVector x);
    //Oblicznie inkluzji kozystajac z F_rest
    capd::IVector computeInclusion(capd::interval time, Algebra::SeriesVector x);
};
#endif 