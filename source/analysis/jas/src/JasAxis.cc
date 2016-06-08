
#ifdef G4ANALYSIS_BUILD_JAS

#include "JasAxis.h"

double JasAxis::lowerEdge() const {return 0;}
double JasAxis::upperEdge() const { return 0;}
int JasAxis::bins() const { return 0;}
double JasAxis::binLowerEdge(int) const { return 0;}
double JasAxis::binUpperEdge(int) const { return 0;}
double JasAxis::binWidth(int) const { return 0;}
double JasAxis::binCentre(int) const { return 0;}
int JasAxis::coordToIndex(double) const { return 0;}
  
#endif
