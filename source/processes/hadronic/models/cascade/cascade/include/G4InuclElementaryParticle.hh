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
// $Id: G4InuclElementaryParticle.hh 72074 2013-07-06 06:28:53Z mkelsey $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100409  M. Kelsey -- Drop unused string argument from ctors.
// 20100429  M. Kelsey -- Change "photon()" to "isPhoton()", use enum names
// 20100914  M. Kelsey -- Move printout to .cc file
// 20100915  M. Kelsey -- Add hyperon() identification function, ctor for
//		G4DynamicParticle
// 20110117  M. Kelsey -- Add antinucleon() and antibaryon() flag
// 20110127  M. Kelsey -- Drop generation.
// 20110214  M. Kelsey -- Replace integer "model" with enum
// 20110321  M. Kelsey -- Fix getStrangeness() to return int
// 20110721  M. Kelsey -- Add constructors to take G4ParticleDefinition as
//		input instead of type code, to allow pass-through of unusable
//		particles during rescattering.  Modify ctors to pass model to
//		base ctor.
// 20110801  M. Kelsey -- Add fill() functions to replicate ctors, allowing
//		reuse of objects as buffers; c.f. G4InuclNuclei.
// 20110922  M. Kelsey -- Add stream argument to printParticle() => print()
// 20120608  M. Kelsey -- Fix variable-name "shadowing" compiler warnings.
// 20130702  M. Kelsey -- Use static type classifiers in G4InuclParticleNames
// 20140310  M. Kelsey -- Fix constness in G4PD* passing

#ifndef G4INUCL_ELEMENTARY_PARTICLE_HH
#define G4INUCL_ELEMENTARY_PARTICLE_HH

#include "G4InuclParticle.hh"
#include "G4InuclParticleNames.hh"
#include "globals.hh"

class G4ParticleDefinition;


class G4InuclElementaryParticle : public G4InuclParticle {
public:
  G4InuclElementaryParticle() 
    : G4InuclParticle() {}

  G4InuclElementaryParticle(G4int ityp, Model model=DefaultModel) 
    : G4InuclParticle(makeDefinition(ityp), model) {}

  G4InuclElementaryParticle(const G4DynamicParticle& dynPart,
			    Model model=DefaultModel)
    : G4InuclParticle(dynPart, model) {}

  G4InuclElementaryParticle(const G4LorentzVector& mom,
			    G4int ityp, Model model=DefaultModel)
    : G4InuclParticle(makeDefinition(ityp), mom, model) {}

  G4InuclElementaryParticle(G4double ekin, G4int ityp,
			    Model model=DefaultModel) 
    : G4InuclParticle(makeDefinition(ityp), ekin, model) {}

  // WARNING:  This may create a particle without a valid type code!
  G4InuclElementaryParticle(const G4LorentzVector& mom,
			    const G4ParticleDefinition* pd,
			    Model model=DefaultModel)
    : G4InuclParticle(pd, mom, model) {}

  // Copy and assignment constructors for use with std::vector<>
  G4InuclElementaryParticle(const G4InuclElementaryParticle& right)
    : G4InuclParticle(right) {}

  G4InuclElementaryParticle& operator=(const G4InuclElementaryParticle& right);

  // Overwrite data structure (avoids creating/copying temporaries)
  void fill(G4int ityp, Model model=DefaultModel) { fill(0., ityp, model); }

  void fill(const G4LorentzVector& mom, G4int ityp, Model model=DefaultModel);

  void fill(G4double ekin, G4int ityp, Model model=DefaultModel);

  // WARNING:  This may create a particle without a valid type code!
  void fill(const G4LorentzVector& mom, const G4ParticleDefinition* pd,
	    Model model=DefaultModel);

  // Assignment and accessor functions
  void setType(G4int ityp);
  G4int type() const { return type(getDefinition()); }

  static G4int type(const G4ParticleDefinition* pd);

  // Ensure that type code refers to a known particle
  inline static G4bool valid(G4int ityp) { return ityp!=0; }
  G4bool valid() const { return valid(type()); }

  G4bool isPhoton() const { return G4InuclParticleNames::isPhoton(type()); }
  G4bool isMuon() const { return G4InuclParticleNames::isMuon(type()); }
  G4bool isElectron() const { return G4InuclParticleNames::isElectron(type()); }
  G4bool isNeutrino() const { return G4InuclParticleNames::isNeutrino(type()); }
  G4bool pion() const { return G4InuclParticleNames::pion(type()); }
  G4bool nucleon() const { return G4InuclParticleNames::nucleon(type()); }
  G4bool antinucleon() const { return G4InuclParticleNames::antinucleon(type()); }

  G4int baryon() const { 		// Can use as a bool (!=0 ==> true)
    return getDefinition()->GetBaryonNumber();
  }

  G4bool antibaryon() const { return baryon() < 0; }

  G4bool hyperon() const { return (baryon() && getStrangeness()); }

  G4bool quasi_deutron() const {
    return G4InuclParticleNames::quasi_deutron(type());
  }

  G4int getStrangeness() const { return getStrangeness(type()); }

  virtual void print(std::ostream& os) const;

  static G4int getStrangeness(G4int type);
  static G4double getParticleMass(G4int type);

protected:
  // Convert internal type code to standard GEANT4 pointer
  static const G4ParticleDefinition* makeDefinition(G4int ityp);
};        

#endif // G4INUCL_ELEMENTARY_PARTICLE_HH 
