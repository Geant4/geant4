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
// $Id: VPreCompoundIon.hh,v 1.1 2003-10-08 12:32:12 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by V. Lara

#ifndef VPreCompoundIon_h
#define VPreCompoundIon_h 1

#include "VPreCompoundFragment.hh"
#include "G4VCoulombBarrier.hh"


class VPreCompoundIon : public VPreCompoundFragment
{
protected:
  // default constructor
  VPreCompoundIon() {}

public:

    // copy constructor
    VPreCompoundIon(const VPreCompoundIon &right): 
	VPreCompoundFragment(right) {}
    
    // constructor  
    VPreCompoundIon(const G4double anA, 
		      const G4double aZ, 
		      G4VCoulombBarrier* aCoulombBarrier,
		      const G4String & aName): 
	VPreCompoundFragment(anA,aZ,aCoulombBarrier,aName) {}
    
    virtual ~VPreCompoundIon() {}
    
    // operators  
    const VPreCompoundIon & 
    operator=(const VPreCompoundIon &right) {
	if (&right != this) this->VPreCompoundFragment::operator=(right);
	return *this;
    }
    
    G4bool operator==(const VPreCompoundIon &right) const 
	{ return VPreCompoundFragment::operator==(right);}
    
    G4bool operator!=(const VPreCompoundIon &right) const 
	{ return VPreCompoundFragment::operator!=(right);}
    
    virtual G4double ProbabilityDistributionFunction(const G4double eKin,
						     const G4Fragment& aFragment,G4double dLevelDensity);
    
protected:
    G4bool IsItPossible(const G4Fragment& aFragment) 
	{
	    G4int pplus = aFragment.GetNumberOfCharged();   
	    G4int pneut = aFragment.GetNumberOfParticles()-pplus;
	    return (pneut >= (GetA()-GetZ()) && pplus >= GetZ());
	}

  virtual G4double GetAlpha() = 0;
  virtual G4double GetBeta(const G4Fragment& fragment) = 0;
  virtual G4double FactorialFactor(const G4double N, const G4double P) = 0;
  virtual G4double CoalescenceFactor(const G4double A) = 0; 
    
};

#endif








