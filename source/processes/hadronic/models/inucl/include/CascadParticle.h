#ifndef CASCAD_PARTICLE_H
#define CASCAD_PARTICLE_H

#include "InuclElementaryParticle.h"

class CascadParticle {

public:

CascadParticle() {};

CascadParticle(const InuclElementaryParticle& particle, 
   const vector<double>& pos,
   int izone, double cpath) :  theParticle(particle), position(pos), 
    current_zone(izone), current_path(cpath) {
  current_path = cpath; 
  movingIn = true;
  reflectionCounter = 0;   
};

void updateParticleMomentum(const vector<double>& mom) {
  theParticle.setMomentum(mom);
};

void updatePosition(const vector<double>& pos) {
  position = pos;
};

void incrementReflectionCounter() { reflectionCounter++; reflected = true; };

void resetReflection() { reflected = false; };

void incrementCurrentPath(double npath) { current_path += npath; };

void updateZone(int izone) { current_zone = izone; };

bool movingInsideNuclei() const { return movingIn; };

double getPathToTheNextZone(double rz_in, double rz_out);

vector<double> getMomentum() const { return theParticle.getMomentum(); };

InuclElementaryParticle getParticle() const { return theParticle; };

vector<double> getPosition() const { return position; };

int getCurrentZone() const { return current_zone; };

int getNumberOfReflections() const { return reflectionCounter; };

bool young(double young_path_cut, double cpath) const { 
//
 if(current_path < 1000.) {
   return cpath < young_path_cut;
 }
  else {
   return false;
 };    
//
// return current_path + cpath < young_path_cut; 
};

bool reflectedNow() const { return reflected; };

void propagateAlongThePath(double path); 

void print() const {
  theParticle.printParticle();
  cout << " zone " << current_zone << " current_path " << current_path
    << " reflectionCounter " << reflectionCounter << endl
    << " x " <<  position[0] << " y " << position[1]
    << " z " << position[2] << endl;
};
   
private: 

InuclElementaryParticle theParticle;
vector<double> position;
int current_zone;
double current_path;
bool movingIn;
int reflectionCounter;   
bool reflected;
 
};        

#endif // CASCAD_PARTICLE_H 
