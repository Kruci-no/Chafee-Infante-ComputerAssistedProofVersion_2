#include <utility>
#include "capd/capdlib.h"
capd::DVector findFixedPoint(capd::DVector initialData,
                        std::function<capd::DVector(capd::DVector)> f,double eps);