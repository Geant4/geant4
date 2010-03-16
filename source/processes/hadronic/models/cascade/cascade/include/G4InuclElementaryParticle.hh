#ifndef G4INUCL_ELEMENTARY_PARTICLE_HH
#define G4INUCL_ELEMENTARY_PARTICLE_HH
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
// $Id: G4InuclElementaryParticle.hh,v 1.20 2010-03-16 22:10:26 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly

#include "G4InuclParticle.hh"
#include "globals.hh"

class G4ParticleDefinition;

class G4InuclElementaryParticle : public G4InuclParticle {

//                     known particle types:
//      1 - proton          11 - k+         111 - quasideuteron PP
//      2 - neutron         13 - k-         112 - quasideuteron PN
//      3 - pi+             15 - k0         122 - quasideuteron NN
//      5 - pi-             17 - k0bar
//      7 - pi 0            21 - lambda 
//     10 - photon          23 - sigma+
//                          25 - sigma0
//                          27 - sigma-
//                          29 - xi0
//                          31 - xi-

public:
  G4InuclElementaryParticle() 
    : G4InuclParticle("InuclElemPart"), generation(0) {}

  explicit G4InuclElementaryParticle(G4int type) 
    : G4InuclParticle("InuclElemPart", makeDefinition(type)), generation(0) {}

  G4InuclElementaryParticle(const G4LorentzVector& mom,
			    G4int type, G4int model=0) 
    : G4InuclParticle("InuclElemPart", makeDefinition(type), mom),
      generation(0) {
    setModel(model);
  }

  G4InuclElementaryParticle(G4double ekin, G4int type) 
    : G4InuclParticle("InuclElemPart", makeDefinition(type), ekin),
      generation(0) {}

  // Copy and assignment constructors for use with std::vector<>
  G4InuclElementaryParticle(const G4InuclElementaryParticle& right)
    : G4InuclParticle(right), generation(right.generation) {}

  G4InuclElementaryParticle& operator=(const G4InuclElementaryParticle& right);

  void setType(G4int ityp);
  G4int type() const { return type(getDefinition()); }

  static G4int type(const G4ParticleDefinition* pd);

  G4bool photon() const { return (type() == 10); }

  G4bool pion() const { return (type()==3 || type()==5 || type()==7); }

  G4bool nucleon() const { return (type() <= 2); }

  G4int baryon() const { 		// Can use as a bool (!=0 ==> true)
    return getDefinition()->GetBaryonNumber();
  }

  G4bool quasi_deutron() const { return (type() > 100); }

  G4bool valid() const { return type()>0; }

  virtual void printParticle() const {
    G4InuclParticle::printParticle();
    G4cout << " Particle: type " << type() << " mass " << getMass()
	   << " ekin " << getKineticEnergy() << G4endl; 
  }

  void setGeneration(G4int gen) { generation = gen; }
  G4int getGeneration() const { return generation; }

  static G4double getStrangeness(G4int type);
  static G4double getParticleMass(G4int type);

protected:
  // Convert internal type code to standard GEANT4 pointer
  static G4ParticleDefinition* makeDefinition(G4int ityp);

private: 
  G4int generation;
};        

#endif // G4INUCL_ELEMENTARY_PARTICLE_HH 
