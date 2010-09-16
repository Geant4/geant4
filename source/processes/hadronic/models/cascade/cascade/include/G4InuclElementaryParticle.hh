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
// $Id: G4InuclElementaryParticle.hh,v 1.25 2010-09-16 05:21:00 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100409  M. Kelsey -- Drop unused string argument from ctors.
// 20100429  M. Kelsey -- Change "photon()" to "isPhoton()", use enum names
// 20100914  M. Kelsey -- Move printout to .cc file
// 20100915  M. Kelsey -- Add hyperon() identification function, ctor for
//		G4DynamicParticle

#ifndef G4INUCL_ELEMENTARY_PARTICLE_HH
#define G4INUCL_ELEMENTARY_PARTICLE_HH

#include "G4InuclParticle.hh"
#include "G4InuclParticleNames.hh"
#include "globals.hh"

class G4ParticleDefinition;

class G4InuclElementaryParticle : public G4InuclParticle {
public:
  G4InuclElementaryParticle() 
    : G4InuclParticle(), generation(0) {}

  explicit G4InuclElementaryParticle(G4int type) 
    : G4InuclParticle(makeDefinition(type)), generation(0) {}

  G4InuclElementaryParticle(const G4DynamicParticle& dynPart, G4int model=0)
    : G4InuclParticle(dynPart), generation(0) {
    setModel(model);
  }

  G4InuclElementaryParticle(const G4LorentzVector& mom,
			    G4int type, G4int model=0) 
    : G4InuclParticle(makeDefinition(type), mom), generation(0) {
    setModel(model);
  }

  G4InuclElementaryParticle(G4double ekin, G4int type) 
    : G4InuclParticle(makeDefinition(type), ekin), generation(0) {}

  // Copy and assignment constructors for use with std::vector<>
  G4InuclElementaryParticle(const G4InuclElementaryParticle& right)
    : G4InuclParticle(right), generation(right.generation) {}

  G4InuclElementaryParticle& operator=(const G4InuclElementaryParticle& right);

  void setType(G4int ityp);
  G4int type() const { return type(getDefinition()); }

  static G4int type(const G4ParticleDefinition* pd);

  G4bool isPhoton() const { return (type() == G4InuclParticleNames::photon); }

  G4bool pion() const { return (type()==G4InuclParticleNames::pionPlus ||
				type()==G4InuclParticleNames::pionMinus ||
				type()==G4InuclParticleNames::pionZero); }

  G4bool nucleon() const { return (type()==G4InuclParticleNames::proton ||
				   type()==G4InuclParticleNames::neutron); }

  G4int baryon() const { 		// Can use as a bool (!=0 ==> true)
    return getDefinition()->GetBaryonNumber();
  }

  G4bool hyperon() const {
    return (baryon() && getStrangeness() != 0.);
  }

  G4bool quasi_deutron() const { return (type() > 100); }

  G4double getStrangeness() const { return getStrangeness(type()); }

  G4bool valid() const { return type()>0; }

  virtual void printParticle() const;

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
