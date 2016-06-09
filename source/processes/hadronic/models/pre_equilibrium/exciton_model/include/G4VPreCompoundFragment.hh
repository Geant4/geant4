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
//
// $Id: G4VPreCompoundFragment.hh,v 1.8.2.1 2009/03/03 13:17:04 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02-patch-04 $
//
// J. M. Quesada (August 2008).  
// Based  on previous work by V. Lara
//
// Modif (03 September 2008) by J. M. Quesada for external choice of inverse 
// cross section option 
// JMQ (06 September 2008) Also external choice has been added for:
//                      - superimposed Coulomb barrier (if useSICB=true) 

#ifndef G4VPreCompoundFragment_h
#define G4VPreCompoundFragment_h 1

#include "G4ios.hh"
#include <iomanip>
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Fragment.hh"
#include "G4VCoulombBarrier.hh"

class G4ReactionProduct;

class G4VPreCompoundFragment
{
  // ============================
  // Constructors and destructor
  // ============================
  
protected:
  // default constructor
  G4VPreCompoundFragment() {};
    
public:  
  // copy constructor
  G4VPreCompoundFragment(const G4VPreCompoundFragment &right);
    
  // constructor  
  G4VPreCompoundFragment(const G4double anA, const G4double aZ,
			 G4VCoulombBarrier * aCoulombBarrier,
			 const G4String &  aName);

  
  virtual ~G4VPreCompoundFragment();
  
  // ==========
  // operators 
  // ========== 
  
  const G4VPreCompoundFragment& 
  operator= (const G4VPreCompoundFragment &right);
  
  G4int operator==(const G4VPreCompoundFragment &right) const;
  
  G4int operator!=(const G4VPreCompoundFragment &right) const;
  
  friend std::ostream& 
  operator<<(std::ostream&, const G4VPreCompoundFragment*);
  friend std::ostream& 
  operator<<(std::ostream&, const G4VPreCompoundFragment&);
  
  // =====================
  // Pure Virtual methods
  // =====================
  virtual G4ReactionProduct * GetReactionProduct() const = 0; 	
  
  // Initialization method
  void Initialize(const G4Fragment & aFragment);
    
  // Methods for calculating the emission probability
  // ------------------------------------------------
  
  // Calculates the total (integrated over kinetic energy) emission
  // probability of a fragment
  virtual G4double CalcEmissionProbability(const G4Fragment & aFragment) = 0;
  
  virtual G4double GetKineticEnergy(const G4Fragment & aFragment) = 0;

public:
  inline G4double GetA() const;
  
  inline G4double GetZ() const;
  
  inline G4double GetRestA() const;
  
  inline G4double GetRestZ() const;
  
  inline G4double GetCoulombBarrier() const;
  
  inline G4double GetBindingEnergy() const;
  
  inline G4double GetMaximalKineticEnergy() const;
  
  inline G4double GetEnergyThreshold() const;

  inline G4double GetEmissionProbability() const;
  
  inline G4double GetNuclearMass() const;
  
  inline G4double GetRestNuclearMass() const;
  
  inline G4double GetReducedMass() const;
  
  inline const G4LorentzVector GetMomentum() const;
  
  inline void  SetMomentum(const G4LorentzVector & value);
  
  inline void  SetFragmentName(const G4String& aName);
  
  inline const G4String GetName() const;
 
  inline void ResetStage();

  inline G4int GetStage() const;

  inline void IncrementStage();

  //for inverse cross section choice
  inline void SetOPTxs(G4int);
  //for superimposed Coulomb Barrier for inverse cross sections
  inline void UseSICB(G4bool);



  // =============
  // Data members
  // =============


private:
  
  G4double theA;
  
  G4double theZ;
private:
  
  G4double theRestNucleusA;
  
  G4double theRestNucleusZ;
protected:  
  G4double theCoulombBarrier;
private:
  G4VCoulombBarrier * theCoulombBarrierPtr;
  
  G4double theBindingEnergy;

  G4double theMaximalKineticEnergy;
  
protected:
  G4double theEmissionProbability;
private:  
  G4LorentzVector theMomentum;
  
  G4String theFragmentName;

  G4int theStage; 

protected:
//for inverse cross section choice
  G4int OPTxs;
//for superimposed Coulomb Barrier for inverse cross sections
  G4bool useSICB;
};

#include "G4VPreCompoundFragment.icc"

#endif
