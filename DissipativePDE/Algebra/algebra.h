#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__
#include <iostream>
#include "capd/capdlib.h"
#include <cmath>
#include <vector>


namespace Algebra{
enum class SeriesType {sin, cos, sin_odd, cos_even};
/*
Series - opisuje nieskończony szereg sinów i cosów
poczatek szeregu jest dany przez wektor main
reszta jest typu C/k^s s moze byc dowolne także ujemne
Typ szeregu opisuje jego indexowanie np szereg sin jest indexowany od 1 a szereg cos od 0
szeregi typu sin_odd i cos_even mają zerowe wartosci na odpowiednio parzystych i nieparzystych pozycjach

*/
struct Series
{
    SeriesType type; //Typ szeregu
    capd::IVector main; //Jawnie dany poczatek szeregu
    int mainSize; //Ilosc jawnie danych elementow
    int n; //Index ostatnie jawnie danego elementu szeregu
    capd::interval C; 
    capd::interval s;

    int static get_n(SeriesType type,int mainSize); 
    // Wylicza index ostatnie jawnie danego elementu 
    // na podstawie typu szeregu i jego dlugosci

    int static getSeriesIndex(SeriesType type,int mainIndex);
    // Wylicza index elementu szeregu z jego pozycji w wektorze main
    Series();
    Series(capd::IVector main,capd::interval C,capd::interval s,SeriesType type);
    Series(capd::interval C,capd::interval s,SeriesType type);
    void print() const;
    capd::interval valueAt(int i) const;
    //Zwraca wartosc i-tego elementu szeregu
    //wgledem indexu szeregu
    Series resize(int newSize);
    //zmienia liczbe jawnie danych elementow na newSize 
    //i zwraca szereg który zawiera stary
    Series upperBound();
    //Zwraca szereg górnych oszacowań na elementy Szeregu
    Series lowerBound();
    //Zwraca szereg dolnych oszacowań na elementy Szeregu
    capd::interval sumAbs();

    capd::interval tailSum();
    //Oszacowywuje sume elemetów danych niejawnie przez C/k^s
    Series refineTail(capd::interval sNew);
    //Zwraca szereg, który zawiera wyjsciowy szereg i jego nowym s jest sNew
    //By funkcja byla poprawnie zdefiniowana musimy miev ze s <= sNew
    capd::interval getNewC(capd::interval sNew) const; 
    //Wylicza jakie bedzie nowe C gdy zmienimy s na sNew 
    //w szeregu wyliczonym przez refine Tail

    bool subset(const Series& x);
    //sprawdza czy wszyskie elementy Szeregu zawieraja sie w SZeregu x
    //wymaga by oba szeregi zawierały 0
    bool subsetInterior(const Series& x);
    //nie działa bez zera
    Series xx();//działa bez zera
    Series elementWiseInverse();//działa bez zera
    Series operator-();// działa bez zera
    friend Series operator*(const capd::interval&  a,const Series& x);//działa bez zera
    friend Series operator*(const Series& x,const capd::interval&  a);//działa bez zera
    friend Series operator*(const Series& x,const Series& y);//działa bez zera ale zwraca szereg z zerem
    friend Series operator+(const Series& x,const Series& y);//działa bez zera ale zwraca szereg z zerem
    friend Series operator-(const Series& x,const Series& y);//działa bez zera ale zwraca szereg z zerem
    friend Series intersection(const Series& x,const Series& y);//działa dla szeregów z takim samym decayej
    friend Series semiIntersection(const Series& x,const Series& y);
    friend Series squere(const Series& x);//działa bez zera ale zwraca szereg z zerem
    friend Series mult(const Series& x,const Series& y);//działa bez zera ale zwraca szereg z zerem
    friend Series elementWiseMult(const Series &x,const Series &y);//działa bez zera
    friend Series hull(const Series &x,const Series &y);
    friend Series expMinusOne(const Series &x);
    friend Series exp(const Series &x,const capd::interval& sNew);
    //friend Series exp(const Series &x,capd::interval sNew);
    
};
struct SeriesVector{
    std::vector<Series> vec;
    SeriesVector(int size,capd::interval C,capd::interval s,SeriesType type);
    SeriesVector(int size);
    SeriesVector();
    void print();
    bool subset(const SeriesVector& x);
    bool subsetInterior(const SeriesVector& x);
    capd::IVector tailSum();
    Series& operator[](int i);
    SeriesVector elementWiseInverse();
    SeriesVector upperBound();
    SeriesVector lowerBound();
    friend SeriesVector operator+(const SeriesVector& x,const SeriesVector& y);
    friend SeriesVector operator*(const capd::interval&  a,const SeriesVector& x);
    friend SeriesVector operator*(const SeriesVector& x,const capd::interval&  a);
    friend SeriesVector elementWiseMult(const SeriesVector& x,const SeriesVector& y);
    friend SeriesVector expMinusOne(const SeriesVector& x);
    friend SeriesVector hull(const SeriesVector& x,const SeriesVector& y);
    friend SeriesVector semiIntersection(const SeriesVector& x,const SeriesVector& y);
    friend SeriesVector exp(const SeriesVector &x,const capd::interval& sNew);   
    };    

}



#endif // __ALGEBRA_H__