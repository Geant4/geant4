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
// by V. Lara

#ifndef G4HETCFragment_h
#define G4HETCFragment_h 1

#include "G4VPreCompoundFragment.hh"

class G4HETCFragment : public G4VPreCompoundFragment
{
protected:
  // default constructor
  G4HETCFragment() {};
    
public:  
  // copy constructor
  G4HETCFragment(const G4HETCFragment &right);
    
  // constructor  
  G4HETCFragment(const G4double anA, const G4double aZ,
		 G4VCoulombBarrier * aCoulombBarrier,
		 const G4String &  aName);
  
  virtual ~G4HETCFragment();
  
  // ==========
  // operators 
  // ========== 
  
  const G4HETCFragment& 
  operator= (const G4HETCFragment &right);
  
  G4int operator==(const G4HETCFragment &right) const;
  
  G4int operator!=(const G4HETCFragment &right) const;
  
protected:

  virtual G4double K(const G4Fragment & aFragment) = 0;
    
  virtual G4double GetSpinFactor() = 0;
  virtual G4double GetAlpha() = 0;
  virtual G4double GetBeta() = 0;

public:
    

  G4double CalcEmissionProbability(const G4Fragment & aFragment);
  
private:	
  // This method performs integration for probability function over 
  // fragment kinetic energy
  G4double IntegrateEmissionProbability(const G4double & Low, 
					const G4double & Up, 
					const G4Fragment & aFragment);	
    
  // ============================
  // Data members access methods
  // ============================
  

  inline G4bool IsItPossible(const G4Fragment & aFragment) const;

protected:
  
  inline G4double BetaRand(const G4int N, const G4int L) const;
  
};

#include "G4HETCFragment.icc"

#endif
