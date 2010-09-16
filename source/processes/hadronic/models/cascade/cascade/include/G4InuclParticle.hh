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
// $Id: G4InuclParticle.hh,v 1.22 2010-09-16 05:21:00 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100112  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100409  M. Kelsey -- Drop unused string argument from ctors.
// 20100519  M. Kelsey -- Add public access to G4DynamicParticle content
// 20100715  M. Kelsey -- Add setKineticEnergy() function
// 20100915  M. Kelsey -- Add constructor to copy G4DynamicParticle input

#ifndef G4INUCL_PARTICLE_HH
#define G4INUCL_PARTICLE_HH

#include "G4DynamicParticle.hh"
#include "G4LorentzVector.hh"
#include "globals.hh"


class G4InuclParticle {
public:
  G4InuclParticle() : modelId(0) {}

  explicit G4InuclParticle(const G4DynamicParticle& dynPart)
    : pDP(dynPart), modelId(0) {}

  explicit G4InuclParticle(const G4LorentzVector& mom) : modelId(0) {
    pDP.Set4Momentum(mom*GeV/MeV);		// From Bertini to G4 units
  }

  virtual ~G4InuclParticle() {}

  // Copy and assignment constructors for use with std::vector<>
  G4InuclParticle(const G4InuclParticle& right)
    : pDP(right.pDP), modelId(right.modelId) {}

  G4InuclParticle& operator=(const G4InuclParticle& right);

  // This is no longer required, as setMomentum() handles mass adjustment
  void setEnergy() { ; }

  // These are call-throughs to G4DynamicParticle
  void setMomentum(const G4LorentzVector& mom);

  void setKineticEnergy(G4double ekin) { pDP.SetKineticEnergy(ekin*GeV/MeV); }

  void setMass(G4double mass) { pDP.SetMass(mass*GeV/MeV); }

  G4double getMass() const {
    return pDP.GetMass()*MeV/GeV;		// From G4 to Bertini units
  }

  G4double getCharge() const {
    return pDP.GetCharge();
  }

  G4double getKineticEnergy() const {
    return pDP.GetKineticEnergy()*MeV/GeV;	// From G4 to Bertini units
  }

  G4double getEnergy() const {
    return pDP.GetTotalEnergy()*MeV/GeV;	// From G4 to Bertini units
  }

  G4double getMomModule() const {
    return pDP.GetTotalMomentum()*MeV/GeV;	// From G4 to Bertini units
  }

  G4LorentzVector getMomentum() const {
    return pDP.Get4Momentum()*MeV/GeV;		// From G4 to Bertini units
  }

  virtual void printParticle() const;

  void setModel(G4int model) { modelId = model; }

  G4int getModel() const { return modelId; }

  G4ParticleDefinition* getDefinition() const {
    return pDP.GetDefinition();
  }

  const G4DynamicParticle& getDynamicParticle() const {
    return pDP;
  }

protected: 
  //  Special constructors for subclasses to set particle type correctly
  explicit G4InuclParticle(G4ParticleDefinition* pd) : modelId(0) {
    setDefinition(pd);
  }

  // FIXME: Bertini code doesn't pass valid 4-vectors, so force mass value
  //	    from supplied PartDefn, with required unit conversions
  G4InuclParticle(G4ParticleDefinition* pd, const G4LorentzVector& mom);

  // NOTE:  Momentum forced along Z direction
  G4InuclParticle(G4ParticleDefinition* pd, G4double ekin)
    : pDP(pd,G4ThreeVector(0.,0.,1.),ekin*GeV/MeV), modelId(0) {}

  void setDefinition(G4ParticleDefinition* pd) { pDP.SetDefinition(pd); }

private:
  G4DynamicParticle pDP;		// Carries all the kinematics and info

  G4int modelId;
  // used to indicate model that created instance of G4InuclParticle
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
