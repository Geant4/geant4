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

#ifndef G4HETCChargedFragment_h
#define G4HETCChargedFragment_h 1

#include "G4HETCFragment.hh"
#include "G4VCoulombBarrier.hh"


class G4HETCChargedFragment : public G4HETCFragment
{
protected:
  // default constructor
  G4HETCChargedFragment() {}

public:

  // copy constructor
  G4HETCChargedFragment(const G4HETCChargedFragment &right): 
    G4HETCFragment(right) {}

  // constructor  
  G4HETCChargedFragment(const G4double anA, 
			const G4double aZ, 
			G4VCoulombBarrier* aCoulombBarrier,
			const G4String & aName): 
    G4HETCFragment(anA,aZ,aCoulombBarrier,aName) {}

  virtual ~G4HETCChargedFragment() {}

  // operators  
  const G4HETCChargedFragment & 
  operator=(const G4HETCChargedFragment &right) 
  {
    if (&right != this) this->G4HETCFragment::operator=(right);
    return *this;
  }

  G4bool operator==(const G4HETCChargedFragment &right) const 
  { 
    return G4HETCFragment::operator==(right);
  }
    
  G4bool operator!=(const G4HETCChargedFragment &right) const 
  { 
    return G4HETCFragment::operator!=(right);
  }
    

  virtual G4double GetKineticEnergy(const G4Fragment & aFragment);
    
    
};

#endif
