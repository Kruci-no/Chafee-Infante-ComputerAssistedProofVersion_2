#include "../../DissipativePDE/VectorFieldMaker/vectorFieldMaker.h"
#include "gallerkinProjections.h"
using namespace std;
using namespace VectorFieldMaker;
//DOROBIC INKLUZJE
string gallerkinProjectionVecField(int size,bool addInclusionTerms){
	std::vector<string> u = makeStringVector(size,"u");
	auto nonLinearity = multiply("(-1)*(A*(sin(omega*t)) + B )", sinTimesCos( u ,sinSquered(u),size));
	nonLinearity.resize(u.size());
	std::vector<string> uLinear = std::vector<string>(size);
	for(int i=0;i<size;i++){
		uLinear[i] = "(("+to_string(-(i+1)*(i+1)) + ")+ lambda)" ;
	}
	
	uLinear = multiply(uLinear,u); 
	std::vector<string> ut = add(uLinear,nonLinearity);
	std::vector<string> parameters(4);
	
	if(!addInclusionTerms){
		parameters[0] = "lambda";parameters[1] = "omega"; parameters[2] = "A"; parameters[3] = "B";
	}
	else{
		parameters.resize(4+size);
		parameters[0] = "lambda";parameters[1] = "omega"; parameters[2] = "A"; parameters[3] = "B";
		std::vector<string> utI = makeStringVector(size,"utI");
		ut = add(ut,utI);
		for(int i=4;i<size + 4;i++){
			parameters[i] = "utI" + to_string(i-4);
		}
	}

	return toFormula(parameters,"t",u,ut);
}
string gallerkinProjectionVecFieldC1(int size,bool addInclusionTerms){
	std::vector<string> u = makeStringVector(size,"u");
	std::vector<string> h = makeStringVector(size,"h");
	
	
	std::vector<string> uLinear = std::vector<string>(size);
	std::vector<string> hLinear = std::vector<string>(size);
	for(int i=0;i<size;i++){
		uLinear[i] = "(("+to_string(-(i+1)*(i+1)) + ")+ lambda)" ;
		hLinear[i] = "(("+to_string(-(i+1)*(i+1)) + ")+ lambda)" ;
	}
	
	uLinear = multiply(uLinear,u); 
	hLinear = multiply(hLinear,h); 
	auto uNonLinearity = multiply("(-1)*(A*(sin(omega*t)) + B )", sinTimesCos( u ,sinSquered(u),size));
	auto hNonLinearity = multiply("(-3)*(A*(sin(omega*t)) + B )", sinTimesCos( h ,sinSquered(u),size));
	uNonLinearity.resize(u.size());
	hNonLinearity.resize(u.size());
	std::vector<string> ut = add(uLinear,uNonLinearity);
	std::vector<string> ht = add(hLinear,hNonLinearity);
	std::vector<string> parameters(4);
	
	if(!addInclusionTerms){
		parameters[0] = "lambda";parameters[1] = "omega"; parameters[2] = "A"; parameters[3] = "B";
	}
	else{
		parameters.resize(4+2*size);
		parameters[0] = "lambda";parameters[1] = "omega"; parameters[2] = "A"; parameters[3] = "B";
		std::vector<string> utI = makeStringVector(size,"utI");
		std::vector<string> htI = makeStringVector(size,"htI");
		ut = add(ut,utI);
		ht = add(ht,htI);
		for(int i=4;i<size + 4;i++){
			parameters[i] = "utI" + to_string(i-4);
			parameters[i+size] = "htI" + to_string(i-4);
		}
	}

	return toFormula(parameters,"t",merge(u,h),merge(ut,ht));
}

string gallerkinProjectionVecFieldSym(int size,bool addInclusionTerms){
	std::vector<string> u = makeStringVectorWithGap(size,2,"u");
	auto nonLinearity = multiply("(-1)*(A*sin(omega*t) + B )", sinTimesCos( u ,sinSquered(u)));
	nonLinearity.resize(u.size());
	std::vector<string> uLinear(u.size());
	for(int i=0;i<uLinear.size();i++){
		uLinear[i] = "(("+to_string(-(i+1)*(i+1)) + ")+ lambda)" ;
	}
	uLinear = multiply(u,uLinear);
	std::vector<string> ut = add(uLinear,nonLinearity);
	std::vector<string> parameters(4);
	if(!addInclusionTerms){
		parameters[0] = "lambda";parameters[1] = "omega"; parameters[2] = "A"; parameters[3] = "B";
	}
	else{
		parameters.resize(4+size);
		parameters[0] = "lambda";parameters[1] = "omega"; parameters[2] = "A"; parameters[3] = "B";
		std::vector<string> utI = makeStringVectorWithGap(size,2,"utI");
		ut = add(ut,utI);
		for(int i=4;i<size + 4;i++){
			parameters[i] = "utI" + to_string(i-4);
		}
	}
	return toFormula(parameters,"t",u,ut);
}
/*
string BrusselatorVectorField(int size) {
	std::vector<string> u = makeStringVector(size,"u");
	std::vector<string> v = makeStringVector(size,"v");
	std::vector<string> uSqeeredv = sinTimesCos(v,sinSquered(u),size);
	std::vector<string> uLinear = std::vector<string>(size);
	for(int i=0;i<size;i++){
		uLinear[i] = "(d1*("+to_string(-(i+1)*(i+1))+")" + "- B - 1)";
	}
	uLinear = multiply(uLinear,u);
	std::vector<string> vLinear = std::vector<string>(size);
	for(int i=0;i<size;i++){
		vLinear[i] = "d2*("+to_string(-(i+1)*(i+1))+")";
	}
	vLinear = add(multiply(vLinear,v), multiply("B", u));
	std::vector<string> perturbation(size);
	perturbation[0]="A";
	std::vector<string> ut = add(add(uLinear,perturbation),uSqeeredv);
	std::vector<string> vt = add(multiply("(-1)",uSqeeredv),vLinear);
	std::vector<string> parameters(4);
	parameters[0] = "d1";parameters[1] = "d2"; parameters[2] = "B"; parameters[3] = "A";
	return toFormula(parameters,shuffle(u,v),shuffle(ut,vt));
}

string BrusselatorVectorFieldSym(int size){
	std::vector<string> u = makeStringVectorWithGap(size,2,"u");
	std::vector<string> v = makeStringVectorWithGap(size,2,"v");
	std::vector<string> uSqeeredv = (sinTimesCos(v,sinSquered(u)));
	uSqeeredv.resize(u.size());
	std::vector<string> uLinear(u.size());
	for(int i=0;i<uLinear.size();i++){
		uLinear[i] = "(d1*("+to_string(-(i+1)*(i+1))+")" + "- B - 1)";
	}
	uLinear = multiply(u,uLinear);
	std::vector<string> vLinear = std::vector<string>(v.size());
	for(int i=0;i<vLinear.size();i++){
		vLinear[i] = "d2*("+to_string(-(i+1)*(i+1))+")";
	}
	vLinear = add(multiply(vLinear,v), multiply("B", u));
	std::vector<string> perturbation(u.size());
	perturbation[0]="A";
	std::vector<string> ut = add(add(uLinear,perturbation),uSqeeredv);
	std::vector<string> vt = add(multiply("(-1)",uSqeeredv),vLinear);

	std::vector<string> parameters(4);
	parameters[0] = "d1";parameters[1] = "d2"; parameters[2] = "B"; parameters[3] = "A";
	return toFormula(parameters,shuffle(u,v),shuffle(ut,vt));
}

string BrusselatorVectorFieldSym2(int size){
	std::vector<string> u = makeStringVectorWithGap(size,2,"u");
	std::vector<string> v = makeStringVectorWithGap(size,2,"v");
	std::vector<string> uSqeeredv = (sinTimesCos(v,sinSquered(u)));
	uSqeeredv.resize(u.size());
	std::vector<string> uLinear(u.size());
	for(int i=0;i<uLinear.size();i++){
		uLinear[i] = "(d1*("+to_string(-(i+1)*(i+1))+")" + "- B - 1)";
	}
	uLinear = multiply(u,uLinear);
	std::vector<string> vLinear = std::vector<string>(v.size());
	for(int i=0;i<vLinear.size();i++){
		vLinear[i] = "d2*("+to_string(-(i+1)*(i+1))+")";
	}
	vLinear = add(multiply(vLinear,v), multiply("B", u));
	std::vector<string> perturbation(u.size());
	perturbation[0]="A";
	std::vector<string> ut = add(add(uLinear,perturbation),uSqeeredv);
	std::vector<string> vt = add(multiply("(-1)",uSqeeredv),vLinear);

	int counter = 0;
	std::vector<string>  parameters(4 + 2*size);
	parameters[0] = "d1";parameters[1] = "d2"; parameters[2] = "B"; parameters[3] = "A";

	std::vector<string> inclusionut(ut.size());
	std::vector<string> inclusionvt(vt.size());
	counter = 0;
	for(int i=0;i<vt.size();i++){
		if(ut[i] != ""){
			inclusionut[i] = "utI"+to_string(counter);
			inclusionvt[i] = "vtI"+to_string(counter);
			parameters[4+2*counter] = "utI"+to_string(counter);
			parameters[4+2*counter+1] = "vtI"+to_string(counter);
			counter = counter + 1;
		}

	}
	ut = add(ut,inclusionut);
	vt = add(vt,inclusionvt);
	return toFormula(parameters,shuffle(u,v),shuffle(ut,vt));
}

string BrusselatorVectorFieldFRest(int size,int mainModesNum){
	std::vector<string> u = makeStringVectorWithGap(size,2,"u");
	std::vector<string> v = makeStringVectorWithGap(size,2,"v");
	std::vector<string> uDiss;
	std::vector<string> vDiss = v;
	std::vector<string> uMain;
	std::vector<string> vMain = v;
	int counter = mainModesNum;
	int i=0;
	split(u,uMain,uDiss,mainModesNum);
	split(v,vMain,vDiss,mainModesNum);
	std::vector<string> uSquaredDiss = add(multiply("2",sinTimesSin(uDiss,uMain)),sinSquered(uDiss));
	std::vector<string> uSquaredMain = sinSquered(uMain);
	std::vector<string> uSquaredvDiss = add(sinTimesCos(vDiss,uSquaredMain),sinTimesCos(v,uSquaredDiss));
	uSquaredvDiss.resize(2*mainModesNum - 1);

	std::vector<string> uSqeeredv = (sinTimesCos(v,sinSquered(u)));
	return toFormula(shuffle(u,v), shuffle(uSquaredvDiss,multiply("(-1)",uSquaredvDiss)) );
}

string BrusselatorVectorFieldFDiss(int size, int mainModesNum){
	std::vector<string> u = makeStringVectorWithGap(size,2,"u");
	std::vector<string> v = makeStringVectorWithGap(size,2,"v");

	std::vector<string> uSqeeredv = (sinTimesCos(v,sinSquered(u)));
	uSqeeredv.resize(u.size());
	std::vector<string> uLinear(u.size());
	for(int i=0;i<uLinear.size();i++){
		uLinear[i] = "";
	}
	uLinear = multiply(u,uLinear);
	std::vector<string> vLinear = std::vector<string>(v.size());
	for(int i=0;i<vLinear.size();i++){
		vLinear[i] = "";
	}
	vLinear = add(multiply(vLinear,v), multiply("B", u));
	std::vector<string> perturbation(u.size());
	perturbation[0]="A";
	std::vector<string> ut = add(add(uLinear,perturbation),uSqeeredv);
	std::vector<string> vt = add(multiply("(-1)",uSqeeredv),vLinear);
	std::vector<string> parameters(4);
	int counter = mainModesNum;
	for(int i=0;i<ut.size();i++){
		if(counter > 0 && ut[i]!=""){
			counter = counter - 1;
			ut[i] = "";
		}
	}
	counter = mainModesNum;
	for(int i=0;i<vt.size();i++){
		if(counter > 0 && vt[i]!=""){
			counter = counter - 1;
			vt[i] = "";
		}
	}
	parameters[0] = "d1";parameters[1] = "d2"; parameters[2] = "B"; parameters[3] = "A";
	return toFormula(parameters,shuffle(u,v),shuffle(ut,vt));
}

string BrusselatorVectorFieldFMain(int size, int mainModesNum){
	std::vector<string> u = makeStringVectorWithGap(size,2,"u");
	std::vector<string> v = makeStringVectorWithGap(size,2,"v");
	std::vector<string> uSqeeredv = (sinTimesCos(v,sinSquered(u)));
	uSqeeredv.resize(u.size());
	std::vector<string> uLinear(u.size());
	for(int i=0;i<uLinear.size();i++){
		uLinear[i] = "(d1*("+to_string(-(i+1)*(i+1))+")" + "- B - 1)";
	}
	uLinear = multiply(u,uLinear);
	std::vector<string> vLinear = std::vector<string>(v.size());
	for(int i=0;i<vLinear.size();i++){
		vLinear[i] = "d2*("+to_string(-(i+1)*(i+1))+")";
	}
	vLinear = add(multiply(vLinear,v), multiply("B", u));
	std::vector<string> perturbation(u.size());
	perturbation[0]="A";
	std::vector<string> ut = add(add(uLinear,perturbation),uSqeeredv);
	std::vector<string> vt = add(multiply("(-1)",uSqeeredv),vLinear);
	int counter=mainModesNum;
	for(int i=0;i<ut.size();i++){
		if(counter>0){
			if(ut[i]!=""){
				counter = counter - 1;
			}
		}
		else{
			ut[i]="";
			vt[i]="";
		}
	}

	std::vector<string> parameters(4);
	parameters[0] = "d1";parameters[1] = "d2"; parameters[2] = "B"; parameters[3] = "A";
	return toFormula(parameters,shuffle(u,v),shuffle(ut,vt));
}


SeriesVector linearBursselator(int size,ParamsMap params){
    SeriesVector L(2);
    Series xx(1,-2,SeriesType::sin_odd);
    Series ones(1,0,SeriesType::sin_odd);
    xx = xx.resize(size);
    ones = ones.resize(size);
    L[0]= (xx * params["d1"] + ones*(params["B"]+interval(1) ) )*interval(-1);
    L[1] = xx * params["d2"] * interval(-1); 
    return L;
}

SeriesVector BursselatorNonlinearity(SeriesVector x,ParamsMap params){
    SeriesVector y(2);
    Series sin (IVector({interval(1)}),0,x[0].s,SeriesType::sin_odd);
    Series u_squere = squere(x[0]).resize(x[1].mainSize + 1);
    Series nonLinearity = mult(u_squere,x[1]);
    y[0]= nonLinearity + params["A"]*sin ;
    y[1] = (interval(-1))*nonLinearity + params["B"]*x[0]; 
    return y;
}

Series u_squere_v(SeriesVector x,ParamsMap params){
    Series u_squere = squere(x[0]).resize(x[1].mainSize + 1);
    Series nonLinearity = mult(u_squere,x[1]);
    return nonLinearity;
}
IVector BrusselatorRest(SeriesVector x,Indexer indexer,ParamsMap params){
    auto splited_x = indexer.splitVector(x);
    SeriesVector xMain = splited_x.first ;
    SeriesVector xDiss = splited_x.second;
    SeriesVector vec(2);
    Series uSquereDiss = (2*mult(xMain[0],xDiss[0]) + squere(xDiss[0])).resize(x[1].mainSize + 1);
    Series uSquereMain = squere(xMain[0]).resize(x[1].mainSize + 1);
    Series uSquere_vDiss = mult(uSquereMain,xDiss[1]) + mult(uSquereDiss,x[1]);
    
    vec.vec[0] = uSquere_vDiss;
    vec.vec[1] = -uSquere_vDiss;
    return indexer.getIVector(vec);
}
Indexer makeBrusselatorIndexer(int u_v_size){
    std::vector<std::pair<int,int> > pairs(2 * u_v_size);
    for(int i=0;i<u_v_size;i++ ){
        pairs[2*i].first = 0;
        pairs[2*i+1].first = 1;
        pairs[2*i].second = i;
        pairs[2*i+1].second = i;
    }
    Indexer indexer;
    indexer.pairs = pairs;
    return indexer;
}
VectorField getBrusselatorVectorField(int u_v_main_size,int serisesSizes,ParamsMap params){
    VectorField vectorField;
    vectorField.indexer = makeBrusselatorIndexer(u_v_main_size);
    vectorField.params = params;
    vectorField.F_nonLinearity = BursselatorNonlinearity;
    vectorField.F_rest = BrusselatorRest;
    vectorField.L = linearBursselator(2*serisesSizes,params);
    vectorField.nMain = vectorField.indexer.size();
    vectorField.numOfNormalParams = 4;
    vectorField.f = IMap(BrusselatorVectorFieldSym2(u_v_main_size));
    vectorField.f.setParameter("A",params["A"]);
    vectorField.f.setParameter("B",params["B"]);
    vectorField.f.setParameter("d1",params["d1"]);
    vectorField.f.setParameter("d2",params["d2"]);
    return vectorField;

}
*/