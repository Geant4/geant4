//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#ifndef G4INUCL_PARTICLE_HH
#define G4INUCL_PARTICLE_HH

#ifndef GLOB
#include "globals.hh"
#endif

#include <iostream>
#include <vector>

// Notice: no cc-file for G4InuclParticle

class G4InuclParticle {

public:
  G4InuclParticle() {
    setModel(0); // default model
  };

  virtual ~G4InuclParticle() { };
 
  G4InuclParticle(const std::vector<G4double>& mom) {
    setMomentum(mom);
    setModel(0);
  };

  void setMomentum(const std::vector<G4double>& mom) {
    momentum = mom;
  };


  std::vector<G4double> getMomentum() const { 
    return momentum; 
  };

  G4double getMomModule() const { 
    return std::sqrt(momentum[1] * momentum[1] +
		     momentum[2] * momentum[2] + 
		     momentum[3] * momentum[3]); 
  };
   
  virtual void printParticle() const {
    G4cout << " px " << momentum[1] << " py " << momentum[2] <<
      " pz " << momentum[3] <<
      " pmod " << std::sqrt(momentum[1] * momentum[1] + 
			    momentum[2] * momentum[2] +
			    momentum[3] * momentum[3])
	   << " E " << momentum[0] 
           << " creator model " << modelId << G4endl;
  };

  void setModel(G4int model) {
    modelId = model;
  };

  G4int getModel() {
    return modelId;
  };

protected: 
  std::vector<G4double> momentum;

private:
  G4int modelId; // used to indicate model that created instance of G4InuclParticle

  // 0 default
  // 1 bullet
  // 2 target
  // 3 G4ElementaryParticleCollider
  // 4 G4IntraNucleiCascader
  // 5 G4NonEquilibriumEvaporator
  // 6 G4EquilibriumEvaporator
  // 7 G4Fissioner
  // 8 G4BigBanger
};        

#endif // G4INUCL_PARTICLE_HH 








