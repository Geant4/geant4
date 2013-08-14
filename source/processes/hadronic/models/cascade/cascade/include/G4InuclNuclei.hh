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
// $Id: G4InuclNuclei.hh 66241 2012-12-13 18:34:42Z gunter $
//
// 20100112  Michael Kelsey -- Replace G4CascadeMomentum with G4LorentzVector
// 20100301  M. Kelsey -- Add function to create unphysical nuclei for use
//	     as temporary final-state fragments.
// 20100319  M. Kelsey -- Remove "using" directory and unnecessary #includes.
// 20100409  M. Kelsey -- Drop unused string argument from ctors.
// 20100630  M. Kelsey -- Add excitation energy as optional public ctor arg,
//	     remove excitation energy data member (part of G4Ions).  Add
//	     excitation energy to "getNucleiMass()" function, move print to .cc
// 20100711  M. Kelsey -- Add optional model ID to constructors
// 20100714  M. Kelsey -- Use G4DynamicParticle::theDynamicalMass to deal with
//	     excitation energy without instantianting "infinite" G4PartDefns.
// 20100719  M. Kelsey -- Move setExitationEnergy implementation to .cc file.
// 20100906  M. Kelsey -- Add fill() functions to rewrite contents
// 20100909  M. Kelsey -- Add function to discard exciton configuration
// 20100914  M. Kelsey -- Use integers for A and Z
// 20100915  M. Kelsey -- Add constructor to copy G4DynamicParticle input
// 20100924  M. Kelsey -- Add constructor to copy G4Fragment input, and output
//		functions to create G4Fragment.
// 20110214  M. Kelsey -- Replace integer "model" with enum
// 20110225  M. Kelsey -- Add equality operator (NOT sorting!)
// 20110721  M. Kelsey -- Follow base-class ctor change to pass model directly
// 20110722  M. Kelsey -- BUG FIX:  Deleted excitation energy in one ctor
// 20110829  M. Kelsey -- Add constructor to copy G4V3DNucleus input
// 20110919  M. Kelsey -- Add clear() to restore completely empty state
// 20110922  M. Kelsey -- Add stream argument to printParticle() => print()

#ifndef G4INUCL_NUCLEI_HH
#define G4INUCL_NUCLEI_HH

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4InuclParticle.hh"
#include "G4LorentzVector.hh"
#include "G4ExitonConfiguration.hh"

class G4Fragment;
class G4ParticleDefinition;
class G4V3DNucleus;


class G4InuclNuclei : public G4InuclParticle {
public:
  G4InuclNuclei() : G4InuclParticle() {}

  G4InuclNuclei(const G4DynamicParticle& dynPart, Model model=DefaultModel)
    : G4InuclParticle(dynPart, model) {}

  G4InuclNuclei(G4int a, G4int z, G4double exc=0., Model model=DefaultModel)
    : G4InuclParticle(makeDefinition(a,z), model) {
    setExitationEnergy(exc);
  }

  G4InuclNuclei(const G4LorentzVector& mom, G4int a, G4int z,
		G4double exc=0., Model model=DefaultModel)
    : G4InuclParticle(makeDefinition(a,z), mom, model) {
    setExitationEnergy(exc);
  }

  G4InuclNuclei(G4double ekin, G4int a, G4int z, G4double exc,
		Model model=DefaultModel)
    : G4InuclParticle(makeDefinition(a,z), ekin, model) {
    setExitationEnergy(exc);
  }

  G4InuclNuclei(const G4Fragment& aFragment, Model model=DefaultModel);

  G4InuclNuclei(G4V3DNucleus* a3DNucleus, Model model=DefaultModel);

  virtual ~G4InuclNuclei() {}

  // Copy and assignment constructors for use with std::vector<>
  G4InuclNuclei(const G4InuclNuclei& right)
    : G4InuclParticle(right),
      theExitonConfiguration(right.theExitonConfiguration) {}

  G4InuclNuclei& operator=(const G4InuclNuclei& right);

  // Equality (comparison) operator -- NOT SORTING
  bool operator==(const G4InuclNuclei& right) {
    return ( G4InuclParticle::operator==(right) && 
	     theExitonConfiguration == right.theExitonConfiguration );
  }

  // Overwrite data structure (avoids creating/copying temporaries)
  void fill(G4int a, G4int z, G4double exc=0., Model model=DefaultModel) {
    fill(0., a, z, exc, model);
  }

  void fill(const G4LorentzVector& mom, G4int a, G4int z,
	    G4double exc=0., Model model=DefaultModel);

  void fill(G4double ekin, G4int a, G4int z, G4double exc,
	    Model model=DefaultModel);

  void copy(const G4Fragment& aFragment, Model model=DefaultModel);

  void copy(G4V3DNucleus* a3DNucleus, Model model=DefaultModel);

  void clear();			// Discard all information (including A,Z)

  // Excitation energy is stored as dynamical mass of particle
  void setExitationEnergy(G4double e);

  void setExitonConfiguration(const G4ExitonConfiguration& config) { 
    theExitonConfiguration = config;
  }

  void clearExitonConfiguration() { theExitonConfiguration.clear(); }

  G4int getA() const { return getDefinition()->GetAtomicMass(); }
  G4int getZ() const { return getDefinition()->GetAtomicNumber(); }

  G4double getNucleiMass() const {
    return getDefinition()->GetPDGMass()*CLHEP::MeV/CLHEP::GeV;	// From G4 to Bertini
  }

  G4double getExitationEnergy() const {
    return (getMass()-getNucleiMass())*CLHEP::GeV/CLHEP::MeV; // Always in MeV
  }

  G4double getExitationEnergyInGeV() const { return getExitationEnergy()/CLHEP::GeV; }

  const G4ExitonConfiguration& getExitonConfiguration() const {
    return theExitonConfiguration;
  }

  static G4double getNucleiMass(G4int a, G4int z, G4double exc=0.);

  virtual void print(std::ostream& os) const;

  // Convert contents to G4Fragment for use outside package
  G4Fragment makeG4Fragment() const;
  operator G4Fragment() const;

protected:
  // Convert nuclear configuration to standard GEANT4 pointer
  static G4ParticleDefinition* makeDefinition(G4int a, G4int z);
  static G4ParticleDefinition* makeNuclearFragment(G4int a, G4int z);

private: 
  G4ExitonConfiguration theExitonConfiguration;
};        

#endif // G4INUCL_NUCLEI_HH 
