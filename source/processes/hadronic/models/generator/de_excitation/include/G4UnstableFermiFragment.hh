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
// $Id: G4UnstableFermiFragment.hh,v 1.9 2003/06/16 17:06:11 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#ifndef G4UnstableFermiFragment_h
#define G4UnstableFermiFragment_h 1

#include <deque>

#include "G4VFermiFragment.hh"
#include "Randomize.hh"

class G4UnstableFermiFragment : public G4VFermiFragment
{
public:
  G4UnstableFermiFragment(const G4int anA, const G4int aZ, const G4int Pol, const G4double ExE):
    G4VFermiFragment(anA,aZ,Pol,ExE)
    {}; 

  ~G4UnstableFermiFragment();
  
protected:
  G4UnstableFermiFragment();
private:
  G4UnstableFermiFragment(const G4UnstableFermiFragment &right);
  
  const G4UnstableFermiFragment & operator=(const G4UnstableFermiFragment &right);
  G4bool operator==(const G4UnstableFermiFragment &right) const;
  G4bool operator!=(const G4UnstableFermiFragment &right) const;
  
public:

  std::deque<G4LorentzVector*> *
  FragmentsMomentum(G4double KinE, const G4int K, const G4double * Masses);

private:

  G4double RNKSI(const G4int K);

  G4ParticleMomentum IsotropicVector(const G4double Magnitude = 1.0);


};


#endif


