/*start time duration time stepNum */
#ifndef sampleDyn_H_
#define sampleDyn_H_
# include "capd/capdlib.h"

string SampleDynamics(capd::DVector initialValue, 
                      capd::DVector params, 
                      capd::DMap map,
                      double startTime, 
                      double durationTime,
                      int steps);

#endif //_EXAMPLES_PROJECTSTARTER_OUTPUT_H_