#ifndef COLLISION_OUTPUT_H
#define COLLISION_OUTPUT_H
#include <iostream.h>

#include "InuclElementaryParticle.h"
#include "InuclNuclei.h"

#include "pair.h"
#include "vector"

class CollisionOutput {

public:

CollisionOutput() {};

void reset() {
  nucleiFragments.resize(0);
  outgoingParticles.resize(0);
  
};

void addOutgoingParticle(const InuclElementaryParticle& particle) {
  outgoingParticles.push_back(particle);
};

void addOutgoingParticles(const vector<InuclElementaryParticle>& particles) {
  for(int i = 0; i < particles.size(); i++)
            outgoingParticles.push_back(particles[i]);
};

void addTargetFragment(const InuclNuclei& nuclei) {
  nucleiFragments.push_back(nuclei);
};

void addTargetFragments(const vector<InuclNuclei>& nuclea) {
  for(int i = 0; i < nuclea.size(); i++)
             nucleiFragments.push_back(nuclea[i]);
};

vector<InuclElementaryParticle> getOutgoingParticles() const {
  return outgoingParticles;
};

int numberOfNucleiFragments() const { return nucleiFragments.size(); };
 
vector<InuclNuclei> getNucleiFragments() const {
  return nucleiFragments;
};

vector<double> getTotalOutputMomentum() const {
  vector<double> tot_mom(4,0.);
  double eex_r = 0.;
  for(int i = 0; i < outgoingParticles.size(); i++) {
    vector<double> mom = outgoingParticles[i].getMomentum();
    for(int j = 0; j < 4; j++) tot_mom[j] += mom[j];
  };
  for(int i = 0; i < nucleiFragments.size(); i++) {
    vector<double> mom = nucleiFragments[i].getMomentum();
    for(int j = 0; j < 4; j++) tot_mom[j] += mom[j];
    eex_r += 0.001*nucleiFragments[i].getExitationEnergy();
  };
  tot_mom[0] += eex_r;
  return tot_mom;
};

void printCollisionOutput() const {
  cout << " Output: " << endl  
       << " Outgoing Particles: " << outgoingParticles.size() << endl;
  for(int i = 0; i < outgoingParticles.size(); i++) {
     outgoingParticles[i].printParticle(); 
  };
  cout << " Nuclei fragments: " << nucleiFragments.size() << endl;      
  for(int i = 0; i < nucleiFragments.size(); i++) {
     nucleiFragments[i].printParticle(); 
  };
};

void trivialise(InuclParticle* bullet, InuclParticle* target) {
  if(InuclNuclei* nuclei_target = dynamic_cast<InuclNuclei*>(target)) {     
    nucleiFragments.push_back(*nuclei_target);
  }
   else {
    InuclElementaryParticle* particle =
               dynamic_cast<InuclElementaryParticle*>(target);
    outgoingParticles.push_back(*particle);
  }; 
  if(InuclNuclei* nuclei_bullet = dynamic_cast<InuclNuclei*>(bullet)) {     
    nucleiFragments.push_back(*nuclei_bullet);
  }
   else {
    InuclElementaryParticle* particle =
               dynamic_cast<InuclElementaryParticle*>(bullet);
    outgoingParticles.push_back(*particle);
  }; 
};

void setOnShell(InuclParticle* bullet, InuclParticle* target);

void setRemainingExitationEnergy() { 
  eex_rest = 0.;
  for(int i = 0; i < nucleiFragments.size(); i++) 
    eex_rest += 0.001*nucleiFragments[i].getExitationEnergy();
};

double getRemainingExitationEnergy() const { return eex_rest; };

bool acceptable() const { return on_shell; };

private: 

vector<InuclElementaryParticle> outgoingParticles;

vector<InuclNuclei> nucleiFragments;

double eex_rest;

pair<pair<int,int>,int> selectPairToTune(double de) const; 

bool on_shell;

};        

#endif // COLLISION_OUTPUT_H 
