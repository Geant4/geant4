#ifndef FISSION_CONFIGURATION_H
#define FISSION_CONFIGURATION_H
#include <iostream.h>

class FissionConfiguration {

public:

FissionConfiguration() {};

FissionConfiguration(double a, double z, double ez, double ek, double ep) :
 afirst(a), zfirst(z), ezet(ez), ekin(ek), epot(ep) {};

void print() {
  cout << " new configuration " << endl
       << " a1 " << afirst << " z1 " << zfirst << " ez " << ezet <<
     " ekin " << ekin << " epot " << epot << endl;
};

double afirst;
double zfirst;
double ezet;
double ekin;
double epot;

};        

#endif // FISSION_CONFIGURATION_H 
