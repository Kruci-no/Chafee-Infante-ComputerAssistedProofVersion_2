#include "solverPDE.h"
using namespace std;
using namespace capd;
using namespace Algebra;
void Encloser::enclose(Set& x,VectorField& vectorField,interval dt,int refineNum = 0,bool constStep = false,bool comPointWiseEnclose = false){
    SeriesVector L = vectorField.L;
    SeriesVector eps(x.vector.vec.size(),interval(-0.01,0.01)*1,x.vector.vec[0].s,x.vector.vec[0].type);
    for(int i=0;i<eps.vec.size();i++){
        eps[i] = Series(interval(-0.01,0.01), x.vector.vec[i].s, x.vector.vec[i].type);
    }
    eps = eps ;
    enclosureExtent = eps;
    enum class Status {nonValidated,validated,end};
    Status status = Status::nonValidated;
    SeriesVector image;
    int iter = 0;
    //cout << "Szukamy Enclosure\n";
    while(status != Status::end){
        
        //cout << iter <<endl;
        if(dt < 1e-50 || iter >100){
            cout << "x";
            x.vector.print();
            cout <<"Enclosure extend";
            enclosureExtent.print();
            cout << iter <<endl ;
            throw std::runtime_error("cannot find enlosure");
        }
        interval timeDuration = x.getCurrentTime() + interval(0,1)*dt ;
        SeriesVector expLminusOne = expMinusOne(L*dt*interval(0,1));
        SeriesVector x_Extend = x.vector + enclosureExtent;
        SeriesVector nonLinearity = vectorField.computeNonLinearity(timeDuration,x_Extend);
        SeriesVector mainVectorField = dt * interval(0,1)*(elementWiseMult(L,x_Extend) + nonLinearity);
        SeriesVector rightBound = elementWiseMult(expLminusOne,x.vector.upperBound()+ elementWiseMult(nonLinearity.upperBound(),L.elementWiseInverse()));
        //rightBound.print();
        SeriesVector leftBound = elementWiseMult(expLminusOne,x.vector.lowerBound()+ elementWiseMult(nonLinearity.lowerBound(),L.elementWiseInverse()));      
        //cout << "LeftBound";
        //leftBound.lowerBound().print();   
        SeriesVector dissEstimation = hull(leftBound.lowerBound(), rightBound.upperBound() );
        //hull = semiIntersection( elementWiseMult(expLminusOne,x.vector + elementWiseMult(nonLinearity,L.elementWiseInverse())) ,hull);
        //cout << "HULL";
        //hull1.print();
        image = semiIntersection(mainVectorField,dissEstimation);
        switch(status){
            case Status::nonValidated:
                iter = iter + 1;
                if(image.subsetInterior(enclosureExtent)){
                    enclosureExtent = image;
                    if(refineNum == 0 )
                
                        status = Status::end;
                    else
                        status = Status::validated;
                    
                }
                else{
                    if(!constStep){
                        dt = dt*interval(1./2);
                    }
                    enclosureExtent = (image + eps )*interval(-0.05,1.1);// eps*interval(0,30)+x.vector*interval(-10,10));
    
                }
                break; 
            case Status::validated: 
                enclosureExtent = semiIntersection(enclosureExtent,image);
                refineNum = refineNum - 1;
                if(refineNum == 0 ){
                    status = Status::end;
                    if(comPointWiseEnclose){
                        SeriesVector dissPointWiseBounds = elementWiseMult(expLminusOne,x.vector + elementWiseMult(nonLinearity,L.elementWiseInverse()));
                        enclosurePointWiseExtent = semiIntersection(mainVectorField,dissPointWiseBounds) ; 
                        for(int i=0; i<enclosurePointWiseExtent.vec.size();i++){
                            enclosurePointWiseExtent[i] = enclosurePointWiseExtent[i].resize(x.vector[i].mainSize);
            
                        }
                    }
                }
                break; 
            case Status::end: 
                break;
        }  
        for(int i=0; i<enclosureExtent.vec.size();i++){
            enclosureExtent[i] = enclosureExtent[i].resize(x.vector[i].mainSize);         
        }
    }
    //cout <<"Enclosure Znalezione :)\n" ;
    validatedTimeStep = dt;
}

Mover::Mover(VectorField& vectorField,Encloser encloser) 
{
    this->step = interval(1./512);
    this->perturbMap = IMap(perturb, vectorField.nMain, vectorField.nMain, vectorField.nMain);
    this->multiMap = new IMultiMap(vectorField.f, this->perturbMap);
    int order = 7;
    solver = new  CWDiffInclSolver( *(this->multiMap) ,order , IMaxNorm() );
}

void Mover::setStep(interval step) 
{
    this->step = step;
}

void Mover::move(Set& x,VectorField& vectorField,bool constStep = false) 
{
    //cout << "Mover" <<endl;
    int order = 7;
    CWDiffInclSolver solver( *(this->multiMap) ,order , IMaxNorm() );
    //interval startTime = x.getCurrentTime();
    //solver = new  CWDiffInclSolver( *(this->multiMap) ,order , IMaxNorm() );
    encloser.enclose(x,vectorField,step,4,constStep);
    SeriesVector x_Extend= x.vector + encloser.enclosureExtent; 
    //encloser.enclosureExtent.print();
    
    interval dt = encloser.validatedTimeStep;
    interval timeDuration = x.getCurrentTime()+interval(0,1)*dt; 
    SeriesVector nonLinearity = vectorField.computeNonLinearity(timeDuration,x_Extend);
    SeriesVector dtL  = dt * vectorField.L;
    //dtL.print();
    SeriesVector nonLinearDuhamelTerm = elementWiseMult(expMinusOne(dtL),
                                                        elementWiseMult( vectorField.L.elementWiseInverse(),nonLinearity));

    /*SeriesVector x_Eval1 = elementWiseMult(exp(dtL,interval(0)) ,x.vector) + nonLinearDuhamelTerm;
    SeriesVector x_Eval2 = elementWiseMult(exp(dtL, -vectorField.L[0].s*(1./4)) ,x.vector) + nonLinearDuhamelTerm;
    */interval x_Eval1Sum = interval(0);
    interval x_Eval2Sum = interval(0);
    for(int i=0; i<x.vector.vec.size();i++){
        
        Series x_Eval1 = elementWiseMult(exp(dtL[i],interval(0)) ,x.vector[i]) + nonLinearDuhamelTerm[i];
        Series x_Eval2 = elementWiseMult(exp(dtL[i], -vectorField.L[i].s*(1./4)) ,x.vector[i]) + nonLinearDuhamelTerm[i];
    
        x_Eval1 = x_Eval1.resize(x.vector[i].mainSize);
        x_Eval2 = x_Eval2.resize(x.vector[i].mainSize);
        x_Eval1Sum += abs((x_Eval1.tailSum()).right()) + abs((x_Eval1.tailSum()).left());
        x_Eval2Sum += abs((x_Eval2.tailSum()).right()) + abs((x_Eval2.tailSum()).left());
        if( (x_Eval1Sum < x_Eval2Sum) || ( (abs(x.vector[i].C)).right() >= 0.1)  ){
            x.vector[i] = x_Eval1;   
            if((abs(x.vector[i].C)).right() >100){
              x.vector[i]  = x.vector[i].refineTail(x.vector[i].s-interval(0.5)) ; 
            }  
        }
        else{
            x.vector[i] = x_Eval2;
        }
    
    }
    
    IVector inclusion = vectorField.computeInclusion(timeDuration,x_Extend); 
    //cout <<"incusion\n:" <<inclusion <<"\n";
    for(int i = 0;i<vectorField.nMain; i++){
        vectorField.f.setParameter(i+vectorField.numOfNormalParams,inclusion[i].mid());
        perturbMap.setParameter(i,inclusion[i] - inclusion[i].mid());
    }
    solver.setStep(dt);
    //cout << "Przesowamy zbior" <<endl;
    x.mainModes.move( solver );
    //cout << "Zbior przesuniety :)" <<endl;
    //cout<< x.getCurrentTime();
    //cout<< x.prevTime+dt;
    //cout<< x.prevTime+dt;
    /*if(x.getCurrentTime()!=startTime+dt){
        cout<< "Time Error";
    }*/
    x.intersectRepresetations(vectorField.indexer);
}

void  Mover::perturb(capd::autodiff::Node t, capd::autodiff::Node in[], int dimIn, 
                        capd::autodiff::Node out[], int dimOut, capd::autodiff::Node params[], int noParams) 
{
    for(int k=0;k<dimIn;k++)
    {
        out[k] = params[k];
    }   
}
