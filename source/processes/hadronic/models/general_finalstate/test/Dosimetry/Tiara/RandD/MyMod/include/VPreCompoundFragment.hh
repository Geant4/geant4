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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: VPreCompoundFragment.hh,v 1.1 2003-10-08 12:32:12 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by V. Lara

#ifndef VPreCompoundFragment_h
#define VPreCompoundFragment_h 1

#include "G4ios.hh"
#include <iomanip>
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Fragment.hh"
#include "G4VCoulombBarrier.hh"

//#define pctest


class G4ReactionProduct;

class VPreCompoundFragment
{
  // ============================
  // Constructors and destructor
  // ============================
  
protected:
  // default constructor
  VPreCompoundFragment() {};
    
public:  
    // copy constructor
    VPreCompoundFragment(const VPreCompoundFragment &right);
    
    // constructor  
    VPreCompoundFragment(const G4double anA, const G4double aZ,
			 G4VCoulombBarrier * aCoulombBarrier,
			 const G4String &  aName);
  
  virtual ~VPreCompoundFragment();
  
  // ==========
  // operators 
  // ========== 
  
  const VPreCompoundFragment& 
  operator= (const VPreCompoundFragment &right);
  
  G4int operator==(const VPreCompoundFragment &right) const;
  
  G4int operator!=(const VPreCompoundFragment &right) const;
  
  friend std::ostream& 
  operator<<(std::ostream&, const VPreCompoundFragment*);
  friend std::ostream& 
  operator<<(std::ostream&, const VPreCompoundFragment&);
    
  // =====================
  // Pure Virtual methods
  // =====================
  virtual G4ReactionProduct * GetReactionProduct() const = 0; 	
  
protected:
  virtual G4double 
  ProbabilityDistributionFunction(const G4double K, 
				  const G4Fragment & aFragment,G4double dLvlDensity) = 0;
    
public:
    
  // Initialization method
  void Initialize(const G4Fragment & aFragment);
    
  // ================================================
  // Methods for calculating the emission probability
  // ================================================
  
  // Calculates the total (integrated over kinetic energy) emission
  // probability of a fragment
  G4double CalcEmissionProbability(const G4Fragment & aFragment,G4double dLvlDensity);
  
  G4double GetKineticEnergy(const G4Fragment & aFragment,G4double dLvlDensity);

private:	
  // This method performs integration for probability function over 
  // fragment kinetic energy
  G4double IntegrateEmissionProbability(const G4double & Low, 
					const G4double & Up, 
					const G4Fragment & aFragment,G4double dLvlDensity);	
    
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

#include "VPreCompoundFragment.icc"

#endif





