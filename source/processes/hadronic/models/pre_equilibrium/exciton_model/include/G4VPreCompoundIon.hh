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
// $Id: G4VPreCompoundIon.hh,v 1.2 2005/06/04 13:48:42 jwellisc Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// by V. Lara

#ifndef G4VPreCompoundIon_h
#define G4VPreCompoundIon_h 1

#include "G4VPreCompoundFragment.hh"
#include "G4VCoulombBarrier.hh"


class G4VPreCompoundIon : public G4VPreCompoundFragment
{
protected:
  // default constructor
  G4VPreCompoundIon() {}

public:

  // copy constructor
  G4VPreCompoundIon(const G4VPreCompoundIon &right): 
    G4VPreCompoundFragment(right) {}
    
  // constructor  
  G4VPreCompoundIon(const G4double anA, 
		    const G4double aZ, 
		    G4VCoulombBarrier* aCoulombBarrier,
		    const G4String & aName): 
    G4VPreCompoundFragment(anA,aZ,aCoulombBarrier,aName) {}
    
  virtual ~G4VPreCompoundIon() {}
    
  // operators  
  const G4VPreCompoundIon & 
  operator=(const G4VPreCompoundIon &right) {
    if (&right != this) this->G4VPreCompoundFragment::operator=(right);
    return *this;
  }
    
  G4bool operator==(const G4VPreCompoundIon &right) const 
  { return G4VPreCompoundFragment::operator==(right);}
    
  G4bool operator!=(const G4VPreCompoundIon &right) const 
  { return G4VPreCompoundFragment::operator!=(right);}
    
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








