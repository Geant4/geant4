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
// $Id: G4InuclNuclei.hh,v 1.19 2010-07-01 22:56:43 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100112  Michael Kelsey -- Replace G4CascadeMomentum with G4LorentzVector
// 20100301  M. Kelsey -- Add function to create unphysical nuclei for use
//	     as temporary final-state fragments.
// 20100319  M. Kelsey -- Remove "using" directory and unnecessary #includes.
// 20100409  M. Kelsey -- Drop unused string argument from ctors.
// 20100630  M. Kelsey -- Add excitation energy as optional public ctor arg,
//	     remove excitation energy data member (part of G4Ions).  Add
//	     excitation energy to "getNucleiMass()" function, move print to .cc

#ifndef G4INUCL_NUCLEI_HH
#define G4INUCL_NUCLEI_HH

#include "G4InuclParticle.hh"
#include "G4LorentzVector.hh"
#include "G4ExitonConfiguration.hh"

class G4ParticleDefinition;


class G4InuclNuclei : public G4InuclParticle {
public:
  G4InuclNuclei() : G4InuclParticle() {}

  G4InuclNuclei(G4double a, G4double z, G4double exc=0.)
    : G4InuclParticle(makeDefinition(a,z,exc)) {}

  G4InuclNuclei(const G4LorentzVector& mom, G4double a, G4double z, G4double exc=0.)
    : G4InuclParticle(makeDefinition(a,z,exc), mom) {}

  G4InuclNuclei(G4double ekin, G4double a, G4double z, G4double exc) 
    : G4InuclParticle(makeDefinition(a,z,exc), ekin) {}

  virtual ~G4InuclNuclei() {}

  // Copy and assignment constructors for use with std::vector<>
  G4InuclNuclei(const G4InuclNuclei& right)
    : G4InuclParticle(right),
      theExitonConfiguration(right.theExitonConfiguration) {}

  G4InuclNuclei& operator=(const G4InuclNuclei& right);

  void setExitationEnergy(G4double e);

  void setExitonConfiguration(const G4ExitonConfiguration& config) { 
    theExitonConfiguration = config;
  }

  G4double getA() const { return getDefinition()->GetAtomicMass(); }

  G4double getZ() const { return getDefinition()->GetAtomicNumber(); }

  G4double getExitationEnergy() const {
    return getExitationEnergy(getDefinition());
  }

  G4double getExitationEnergyInGeV() const { return getExitationEnergy()/GeV; }

  const G4ExitonConfiguration& getExitonConfiguration() const {
    return theExitonConfiguration;
  }

  static G4double getNucleiMass(G4double a, G4double z, G4double exc=0.);

  static G4double getExitationEnergy(const G4ParticleDefinition* pd);

  virtual void printParticle() const;

protected:
  // Convert nuclear configuration to standard GEANT4 pointer
  static G4ParticleDefinition*
  makeDefinition(G4double a, G4double z, G4double exc);

  static G4ParticleDefinition* 
  makeNuclearFragment(G4double a, G4double z, G4double exc);

private: 
  G4ExitonConfiguration theExitonConfiguration;
};        

#endif // G4INUCL_NUCLEI_HH 






