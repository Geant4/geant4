#include <math.h>
#include "pair.h"
#include "vector"

namespace CascadSpecialFunctions {

pair<int,double> getPositionInEnergyScale2(double e); 

pair<int,double> getPositionInEnergyScale1(double e);
 
double absorptionCrosSection(double e, int type);

double crossSection(double e, int is);

pair<int,double> getPositionInEnergyScaleEMS(double e); 
 
}
