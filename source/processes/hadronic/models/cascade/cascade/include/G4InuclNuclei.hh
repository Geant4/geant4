#ifndef G4INUCL_NUCLEI_HH
#define G4INUCL_NUCLEI_HH
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
// $Id: G4InuclNuclei.hh,v 1.16 2010-04-09 19:33:11 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100112  Michael Kelsey -- Replace G4CascadeMomentum with G4LorentzVector
// 20100301  M. Kelsey -- Add function to create unphysical nuclei for use
//	     as temporary final-state fragments.
// 20100319  M. Kelsey -- Remove "using" directory and unnecessary #includes.
// 20100409  M. Kelsey -- Drop unused string argument from ctors.

#include "G4InuclParticle.hh"
#include "G4ExitonConfiguration.hh"

class G4ParticleDefinition;

class G4InuclNuclei : public G4InuclParticle {
public:
  G4InuclNuclei() : G4InuclParticle() {}

  G4InuclNuclei(G4double a, G4double z)
    : G4InuclParticle(makeDefinition(a,z)),
      exitationEnergy(0.0) {}

  G4InuclNuclei(const G4LorentzVector& mom, G4double a, G4double z)
    : G4InuclParticle(makeDefinition(a,z), mom),
      exitationEnergy(0.0) {}

  G4InuclNuclei(G4double ekin, G4double a, G4double z) 
    : G4InuclParticle(makeDefinition(a,z), ekin),
      exitationEnergy(0.0) {}

  virtual ~G4InuclNuclei() {}

  // Copy and assignment constructors for use with std::vector<>
  G4InuclNuclei(const G4InuclNuclei& right)
    : G4InuclParticle(right), exitationEnergy(right.exitationEnergy),
      theExitonConfiguration(right.theExitonConfiguration) {}

  G4InuclNuclei& operator=(const G4InuclNuclei& right);

  void setExitationEnergy(G4double e) { exitationEnergy = e; }

  void setExitonConfiguration(const G4ExitonConfiguration& config) { 
    theExitonConfiguration = config;
  }

  G4double getA() const { return getDefinition()->GetAtomicMass(); }

  G4double getZ() const { return getDefinition()->GetAtomicNumber(); }

  G4double getExitationEnergy() const { return exitationEnergy; }

  G4double getExitationEnergyInGeV() const { 
    return 0.001 * exitationEnergy; 
  }

  const G4ExitonConfiguration& getExitonConfiguration() const {
    return theExitonConfiguration;
  }

  static G4double getNucleiMass(G4double a, G4double z);

  virtual void printParticle() const {
    G4cout << " A " << getA() << " Z " << getZ() << " mass " 
	   << getMass() << " Eex (MeV) " << exitationEnergy << G4endl;
    G4InuclParticle::printParticle();
  }

protected:
  // Convert nuclear configuration to standard GEANT4 pointer
  static G4ParticleDefinition*
  makeDefinition(G4double a, G4double z, G4double exc=0.);

  static G4ParticleDefinition* 
  makeNuclearFragment(G4double a, G4double z, G4double exc=0.);

private: 
  G4double exitationEnergy;
  G4ExitonConfiguration theExitonConfiguration;
};        

#endif // G4INUCL_NUCLEI_HH 






