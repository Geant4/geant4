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

#ifndef G4PreCompoundNucleon_h
#define G4PreCompoundNucleon_h 1

#include "G4PreCompoundFragment.hh"

class G4PreCompoundNucleon : public G4PreCompoundFragment
{
private:
  // default constructor
  G4PreCompoundNucleon() {};

public:
  
  // copy constructor
  G4PreCompoundNucleon(const G4PreCompoundNucleon &right): 
    G4PreCompoundFragment(right) {}

  // constructor  
  G4PreCompoundNucleon(const G4double anA, 
		       const G4double aZ, 
		       G4VCoulombBarrier* aCoulombBarrier,
		       const G4String & aName): 
    G4PreCompoundFragment(anA,aZ,aCoulombBarrier,aName) {}
  
  virtual ~G4PreCompoundNucleon() {}

  // operators  
  const G4PreCompoundNucleon & 
  operator=(const G4PreCompoundNucleon &right) 
  {
    if (&right != this) this->G4PreCompoundFragment::operator=(right);
    return *this;
  }
  
  G4bool operator==(const G4PreCompoundNucleon &right) const 
  { 
    return G4PreCompoundFragment::operator==(right);
  }
    
  G4bool operator!=(const G4PreCompoundNucleon &right) const 
  { 
    return G4PreCompoundFragment::operator!=(right);
  }
    
  virtual G4double ProbabilityDistributionFunction(const G4double eKin,
						   const G4Fragment& aFragment);

    
protected:

  virtual G4double GetAlpha() = 0;
  virtual G4double GetBeta() = 0;
  virtual G4bool IsItPossible(const G4Fragment&) = 0;     
};

#endif
