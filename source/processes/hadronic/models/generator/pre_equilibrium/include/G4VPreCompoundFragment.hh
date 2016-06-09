//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VPreCompoundFragment.hh,v 1.13 2002/12/12 19:17:32 gunter Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
// by V. Lara

#ifndef G4VPreCompoundFragment_h
#define G4VPreCompoundFragment_h 1

#include "G4ios.hh"
#include "g4std/iomanip"
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
  
  friend G4std::ostream& 
  operator<<(G4std::ostream&, const G4VPreCompoundFragment*);
  friend G4std::ostream& 
  operator<<(G4std::ostream&, const G4VPreCompoundFragment&);
    
  // =====================
  // Pure Virtual methods
  // =====================
  virtual G4ReactionProduct * GetReactionProduct() const = 0; 	
  
protected:
  virtual G4double 
  ProbabilityDistributionFunction(const G4double K, 
				  const G4Fragment & aFragment) = 0;
    
public:
    
  // Initialization method
  void Initialize(const G4Fragment & aFragment);
    
  // ================================================
  // Methods for calculating the emission probability
  // ================================================
  
  // Calculates the total (integrated over kinetic energy) emission
  // probability of a fragment
  G4double CalcEmissionProbability(const G4Fragment & aFragment);
  
  G4double GetKineticEnergy(const G4Fragment & aFragment);

private:	
  // This method performs integration for probability function over 
  // fragment kinetic energy
  G4double IntegrateEmissionProbability(const G4double & Low, 
					const G4double & Up, 
					const G4Fragment & aFragment);	
    
    // ============================
    // Data members access methods
    // ============================
    
public:
  inline const G4double GetA() const;
  
  inline const G4double GetZ() const;
  
  inline const G4double GetRestA() const;
  
  inline const G4double GetRestZ() const;
  
  inline const G4double GetCoulombBarrier() const;
  
  inline const G4double GetBindingEnergy() const;
  
  inline const G4double GetMaximalKineticEnergy() const;
  
  inline const G4double GetEmissionProbability() const;
  
  inline const G4double GetNuclearMass() const;
  
  inline const G4double GetRestNuclearMass() const;
  
  inline const G4double GetReducedMass() const;
  
  inline const G4LorentzVector GetMomentum() const;
  
  inline void SetMomentum(const G4LorentzVector & value);
  
  inline void SetFragmentName(const G4String& aName);
  
  inline const G4String GetName() const;
  
  
  // =============
  // Data members
  // =============
  
private:
  
  G4double theA;
  
  G4double theZ;
  
  G4double theRestNucleusA;
  
  G4double theRestNucleusZ;
  
  G4double theCoulombBarrier;
  
  G4VCoulombBarrier * theCoulombBarrierPtr;
  
  G4double theBindingEnergy;
  
  G4double theMaximalKineticEnergy;
  
  G4double theEmissionProbability;
  
  G4LorentzVector theMomentum;
  
  G4String theFragmentName;
  
};

#include "G4VPreCompoundFragment.icc"

#endif





