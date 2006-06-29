//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4UnstableFermiFragment.hh,v 1.3 2006-06-29 20:12:33 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#ifndef G4UnstableFermiFragment_h
#define G4UnstableFermiFragment_h 1

#include "G4VFermiFragment.hh"
#include "Randomize.hh"

class G4UnstableFermiFragment : public G4VFermiFragment
{
public:
  virtual ~G4UnstableFermiFragment();
  
protected:
  G4UnstableFermiFragment(const G4int anA, const G4int aZ, const G4int Pol, const G4double ExE)
    : G4VFermiFragment(anA,aZ,Pol,ExE)
  {
  } 

  G4UnstableFermiFragment();

private:
  G4UnstableFermiFragment(const G4UnstableFermiFragment &right);
  
  const G4UnstableFermiFragment & operator=(const G4UnstableFermiFragment &right);
  G4bool operator==(const G4UnstableFermiFragment &right) const;
  G4bool operator!=(const G4UnstableFermiFragment &right) const;
  
public:

  G4FragmentVector * GetFragment(const G4LorentzVector&) const;

protected:

  std::vector<G4double> Masses;
  std::vector<G4double> Charges;
  std::vector<G4double> AtomNum;

};


#endif


