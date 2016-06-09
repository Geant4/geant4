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

#ifndef G4PreCompoundIon_h
#define G4PreCompoundIon_h 1

#include "G4PreCompoundFragment.hh"

class G4PreCompoundIon : public G4PreCompoundFragment
{
protected:
  // default constructor
  G4PreCompoundIon() {}

public:

  // copy constructor
  G4PreCompoundIon(const G4PreCompoundIon &right): 
    G4PreCompoundFragment(right) {}
    
  // constructor  
  G4PreCompoundIon(const G4double anA, 
		   const G4double aZ, 
		   G4VCoulombBarrier* aCoulombBarrier,
		   const G4String & aName): 
    G4PreCompoundFragment(anA,aZ,aCoulombBarrier,aName) {}
    
  virtual ~G4PreCompoundIon() {}
    
  // operators  
  const G4PreCompoundIon & 
  operator=(const G4PreCompoundIon &right) 
  {
    if (&right != this) this->G4PreCompoundFragment::operator=(right);
    return *this;
  }
    
  G4bool operator==(const G4PreCompoundIon &right) const 
  { 
    return G4PreCompoundFragment::operator==(right);
  }
    
  G4bool operator!=(const G4PreCompoundIon &right) const 
  { 
    return G4PreCompoundFragment::operator!=(right);
  }
    
  virtual G4double ProbabilityDistributionFunction(const G4double eKin,
						   const G4Fragment& aFragment);

protected:
  G4bool IsItPossible(const G4Fragment& aFragment) 
  {
    G4int pplus = aFragment.GetNumberOfCharged();   
    G4int pneut = aFragment.GetNumberOfParticles()-pplus;
    return (pneut >= (GetA()-GetZ()) && pplus >= GetZ());
  }
  
  virtual G4double GetAlpha() = 0;
  virtual G4double GetBeta() = 0;
  virtual G4double FactorialFactor(const G4double N, const G4double P) = 0;
  virtual G4double CoalescenceFactor(const G4double A) = 0; 
    
};

#endif
