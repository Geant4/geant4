#ifndef INUCL_PARTICLE_H
#define INUCL_PARTICLE_H

#include <iostream.h>
#include "vector"

class InuclParticle {

public:

InuclParticle() {};

InuclParticle(const vector<double>& mom) {
  setMomentum(mom);
};

void setMomentum(const vector<double>& mom) {
  momentum = mom;
};

vector<double> getMomentum() const { return momentum; };

double getMomModule() const { return sqrt(momentum[1]*momentum[1] +
  momentum[2]*momentum[2] + momentum[3]*momentum[3]); };
   
virtual void printParticle() const {
  cout << " px " << momentum[1] << " py " << momentum[2] <<
      " pz " << momentum[3] <<
      " pmod " << sqrt(momentum[1]*momentum[1] + momentum[2]*momentum[2] +
        momentum[3]*momentum[3])
      << " E " << momentum[0] << endl;
};

protected: 

vector<double> momentum;

};        

#endif // INUCL_PARTICLE_H 
