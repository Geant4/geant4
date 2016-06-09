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
// $Id: G4StableFermiFragment.hh,v 1.1 2003/08/26 18:34:16 lara Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#ifndef G4StableFermiFragment_h
#define G4StableFermiFragment_h 1

#include "G4VFermiFragment.hh"

class G4StableFermiFragment : public G4VFermiFragment 
{
public:
  G4StableFermiFragment(const G4int anA, const G4int aZ, const G4int Pol, const G4double ExE)
    : G4VFermiFragment(anA,aZ,Pol,ExE)
    {
    }

  virtual ~G4StableFermiFragment();
  
private:
  G4StableFermiFragment();

  G4StableFermiFragment(const G4StableFermiFragment &right);
  
  const G4StableFermiFragment & operator=(const G4StableFermiFragment &right);
  G4bool operator==(const G4StableFermiFragment &right) const;
  G4bool operator!=(const G4StableFermiFragment &right) const;
  
public:

  virtual G4FragmentVector * GetFragment(const G4LorentzVector & aMomentum) const;

};


#endif


