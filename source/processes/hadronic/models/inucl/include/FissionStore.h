#ifndef FISSION_STORE_H
#define FISSION_STORE_H

#include "FissionConfiguration.h"
#include "vector"

class FissionStore {

public:

FissionStore() {};

void addConfig(double a, double z, double ez, double ek, double ev) {
  FissionConfiguration config(a,z,ez,ek,ev);
  configurations.push_back(config);
//  config.print();
};

int size() const { return configurations.size(); };

FissionConfiguration generateConfiguration(double amax, double rand) const;

private:
 
vector<FissionConfiguration> configurations;

};

#endif // FISSION_STORE_H 
