//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
#ifndef G4COLLISION_OUTPUT_HH
#define G4COLLISION_OUTPUT_HH

#include "g4std/iostream"

#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"

#include "g4std/algorithm"
#include "g4std/vector"

class G4CollisionOutput {

public:

  G4CollisionOutput();

  void reset() {
    nucleiFragments.resize(0);
    outgoingParticles.resize(0);
  
  };

  void addOutgoingParticle(const G4InuclElementaryParticle& particle) {
    outgoingParticles.push_back(particle);
  };

  void addOutgoingParticles(const G4std::vector<G4InuclElementaryParticle>& particles) {
    for(G4int i = 0; i < G4int(particles.size()); i++)
      outgoingParticles.push_back(particles[i]);
  };

  void addTargetFragment(const G4InuclNuclei& nuclei) {
    nucleiFragments.push_back(nuclei);
  };

  void addTargetFragments(const G4std::vector<G4InuclNuclei>& nuclea) {
    for(G4int i = 0; i < G4int(nuclea.size()); i++)
      nucleiFragments.push_back(nuclea[i]);
  };

 G4std::vector<G4InuclElementaryParticle> getOutgoingParticles() const {
    return outgoingParticles;
  };

  G4int numberOfNucleiFragments() const { 
    return nucleiFragments.size(); 
  };
 
  G4std::vector<G4InuclNuclei> getNucleiFragments() const {
    return nucleiFragments;
  };

  G4std::vector<G4double> getTotalOutputMomentum() const {
    G4std::vector<G4double> tot_mom(4, 0.0);
    double eex_r = 0.0;
    G4int i(0);
    for(i = 0; i < G4int(outgoingParticles.size()); i++) {
      G4std::vector<G4double> mom = outgoingParticles[i].getMomentum();
      for(G4int j = 0; j < 4; j++) tot_mom[j] += mom[j];
    };
    for(i = 0; i < G4int(nucleiFragments.size()); i++) {
      G4std::vector<G4double> mom = nucleiFragments[i].getMomentum();
      for(G4int j = 0; j < 4; j++) tot_mom[j] += mom[j];
      eex_r += 0.001 * nucleiFragments[i].getExitationEnergy();
    };
    tot_mom[0] += eex_r;
    return tot_mom;
  };

  void printCollisionOutput() const {
    G4cout << " Output: " << G4endl  
	   << " Outgoing Particles: " << outgoingParticles.size() << G4endl;
    G4int i(0);
    for(i = 0; i < G4int(outgoingParticles.size()); i++) {
      outgoingParticles[i].printParticle(); 
    };
    G4cout << " Nuclei fragments: " << nucleiFragments.size() << G4endl;      
    for(i = 0; i < G4int(nucleiFragments.size()); i++) {
      nucleiFragments[i].printParticle(); 
    };
  };

  void trivialise(G4InuclParticle* bullet, 
		  G4InuclParticle* target) {
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

  void setOnShell(G4InuclParticle* bullet, 
		  G4InuclParticle* target);

  void setRemainingExitationEnergy() { 
    eex_rest = 0.0;
    for(G4int i = 0; i < G4int(nucleiFragments.size()); i++) 
      eex_rest += 0.001 * nucleiFragments[i].getExitationEnergy();
  };

  double getRemainingExitationEnergy() const { 
    return eex_rest; 
  };

  G4bool acceptable() const { 
    return on_shell; 
  };

private: 

G4int verboseLevel;
  G4std::vector<G4InuclElementaryParticle> outgoingParticles;

  G4std::vector<G4InuclNuclei> nucleiFragments;

  G4double eex_rest;

  G4std::pair<G4std::pair<G4int, G4int>, G4int> selectPairToTune(G4double de) const; 

  G4bool on_shell;

};        

#endif // G4COLLISION_OUTPUT_HH 

