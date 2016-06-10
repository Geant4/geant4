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
// $Id: G4VPreCompoundFragment.hh 68028 2013-03-13 13:48:15Z gcosmo $
//
// J. M. Quesada (August 2008).  
// Based  on previous work by V. Lara
//
// Modif (03 September 2008) by J. M. Quesada for external choice of inverse 
// cross section option 
// JMQ (06 September 2008) Also external choice has been added for:
//                      - superimposed Coulomb barrier (if useSICB=true) 
// 20.08.2010 V.Ivanchenko added int Z and A and cleanup; added 
//                        G4ParticleDefinition to constructor, 
//                        inline method to build G4ReactionProduct; 
//                        remove string name
//                        

#ifndef G4VPreCompoundFragment_h
#define G4VPreCompoundFragment_h 1

#include "G4ios.hh"
#include <iomanip>
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Fragment.hh"
#include "G4VCoulombBarrier.hh"
#include "G4ReactionProduct.hh"
#include "G4PreCompoundParameters.hh"
#include "G4Pow.hh"

class G4VPreCompoundFragment
{
public:  

  // ============================
  // Constructors and destructor
  // ============================
    
  G4VPreCompoundFragment(const G4ParticleDefinition*,
			 G4VCoulombBarrier * aCoulombBarrier);
  
  virtual ~G4VPreCompoundFragment();
  
  // ==========
  // operators 
  // ========== 
  
  friend std::ostream& 
  operator<<(std::ostream&, const G4VPreCompoundFragment*);
  friend std::ostream& 
  operator<<(std::ostream&, const G4VPreCompoundFragment&);
  
  // =====================
  // Pure Virtual methods
  // =====================
  
  // Initialization method
  void Initialize(const G4Fragment & aFragment);
    
  // Methods for calculating the emission probability
  // ------------------------------------------------
  
  // Calculates the total (integrated over kinetic energy) emission
  // probability of a fragment
  virtual G4double CalcEmissionProbability(const G4Fragment & aFragment) = 0;
  
  virtual G4double GetKineticEnergy(const G4Fragment & aFragment) = 0;

  inline G4ReactionProduct * GetReactionProduct() const; 	

  inline G4int GetA() const;
  
  inline G4int GetZ() const;
  
  inline G4int GetRestA() const;
  
  inline G4int GetRestZ() const;

  inline G4double ResidualA13() const;
  
  inline G4double GetCoulombBarrier() const;
  
  inline G4double GetBindingEnergy() const;
  
  inline G4double GetMaximalKineticEnergy() const;
  
  inline G4double GetEnergyThreshold() const;

  inline G4double GetEmissionProbability() const;
  
  inline G4double GetNuclearMass() const;
  
  inline G4double GetRestNuclearMass() const;
  
  inline G4double GetReducedMass() const;
  
  inline const G4LorentzVector& GetMomentum() const;
  
  inline void  SetMomentum(const G4LorentzVector & value);
  
  inline const G4String GetName() const;

  //for inverse cross section choice
  inline void SetOPTxs(G4int);
  //for superimposed Coulomb Barrier for inverse cross sections
  inline void UseSICB(G4bool);

protected:

  inline G4bool IsItPossible(const G4Fragment & aFragment) const;

private:

  // default constructor
  G4VPreCompoundFragment();
  // copy constructor
  G4VPreCompoundFragment(const G4VPreCompoundFragment &right);
  const G4VPreCompoundFragment& 
  operator= (const G4VPreCompoundFragment &right);  
  G4int operator==(const G4VPreCompoundFragment &right) const;
  G4int operator!=(const G4VPreCompoundFragment &right) const;

  // =============
  // Data members
  // =============

  const G4ParticleDefinition* particle;
  G4VCoulombBarrier * theCoulombBarrierPtr;
  
  G4int theA;
  G4int theZ;
  G4int theRestNucleusA;
  G4int theRestNucleusZ;

  G4double theRestNucleusA13;
  G4double theBindingEnergy;
  G4double theMaximalKineticEnergy;
  G4double theRestNucleusMass;
  G4double theReducedMass;
  G4double theMass;

  G4LorentzVector theMomentum;
  
protected:

  G4PreCompoundParameters* theParameters;
  G4Pow* g4pow;

  G4double theEmissionProbability;
  G4double theCoulombBarrier;

  //for inverse cross section choice
  G4int OPTxs;
  //for superimposed Coulomb Barrier for inverse cross sections
  G4bool useSICB;
};

#include "G4VPreCompoundFragment.icc"

#endif
