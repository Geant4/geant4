#ifndef G4COLLISION_OUTPUT_HH
#define G4COLLISION_OUTPUT_HH
#include <iostream.h>

#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"

#include "pair.h"
#include "vector"

class G4CollisionOutput {

public:

G4CollisionOutput() {};

void reset() {
  nucleiFragments.resize(0);
  outgoingParticles.resize(0);
  
};

void addOutgoingParticle(const G4InuclElementaryParticle& particle) {
  outgoingParticles.push_back(particle);
};

void addOutgoingParticles(const vector<G4InuclElementaryParticle>& particles) {
  for(G4int i = 0; i < particles.size(); i++)
            outgoingParticles.push_back(particles[i]);
};

void addTargetFragment(const G4InuclNuclei& nuclei) {
  nucleiFragments.push_back(nuclei);
};

void addTargetFragments(const vector<G4InuclNuclei>& nuclea) {
  for(G4int i = 0; i < nuclea.size(); i++)
             nucleiFragments.push_back(nuclea[i]);
};

vector<G4InuclElementaryParticle> getOutgoingParticles() const {
  return outgoingParticles;
};

G4int numberOfNucleiFragments() const { return nucleiFragments.size(); };
 
vector<G4InuclNuclei> getNucleiFragments() const {
  return nucleiFragments;
};

vector<G4double> getTotalOutputMomentum() const {
  vector<G4double> tot_mom(4,0.0);
  double eex_r = 0.0;
  for(G4int i = 0; i < outgoingParticles.size(); i++) {
    vector<G4double> mom = outgoingParticles[i].getMomentum();
    for(G4int j = 0; j < 4; j++) tot_mom[j] += mom[j];
  };
  for(G4int i = 0; i < nucleiFragments.size(); i++) {
    vector<G4double> mom = nucleiFragments[i].getMomentum();
    for(G4int j = 0; j < 4; j++) tot_mom[j] += mom[j];
    eex_r += 0.001*nucleiFragments[i].getExitationEnergy();
  };
  tot_mom[0] += eex_r;
  return tot_mom;
};

void printCollisionOutput() const {
  G4cout << " Output: " << endl  
       << " Outgoing Particles: " << outgoingParticles.size() << G4endl;
  for(G4int i = 0; i < outgoingParticles.size(); i++) {
     outgoingParticles[i].printParticle(); 
  };
  G4cout << " Nuclei fragments: " << nucleiFragments.size() << G4endl;      
  for(G4int i = 0; i < nucleiFragments.size(); i++) {
     nucleiFragments[i].printParticle(); 
  };
};

void trivialise(G4InuclParticle* bullet, G4InuclParticle* target) {
  if(G4InuclNuclei* nuclei_target = dynamic_cast<G4InuclNuclei*>(target)) {     
    nucleiFragments.push_back(*nuclei_target);
  }
   else {
    G4InuclElementaryParticle* particle =
               dynamic_cast<G4InuclElementaryParticle*>(target);
    outgoingParticles.push_back(*particle);
  }; 
  if(G4InuclNuclei* nuclei_bullet = dynamic_cast<G4InuclNuclei*>(bullet)) {     
    nucleiFragments.push_back(*nuclei_bullet);
  }
   else {
    G4InuclElementaryParticle* particle =
               dynamic_cast<G4InuclElementaryParticle*>(bullet);
    outgoingParticles.push_back(*particle);
  }; 
};

void setOnShell(G4InuclParticle* bullet, G4InuclParticle* target);

void setRemainingExitationEnergy() { 
  eex_rest = 0.0;
  for(G4int i = 0; i < nucleiFragments.size(); i++) 
    eex_rest += 0.001*nucleiFragments[i].getExitationEnergy();
};

double getRemainingExitationEnergy() const { return eex_rest; };

bool acceptable() const { return on_shell; };

private: 

vector<G4InuclElementaryParticle> outgoingParticles;

vector<G4InuclNuclei> nucleiFragments;

G4double eex_rest;

pair<pair<G4int, G4int>, G4int> selectPairToTune(G4double de) const; 

bool on_shell;

};        

#endif // G4COLLISION_OUTPUT_HH 

