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
// $Id: G4VFermiFragment.hh,v 1.3 2006/06/29 20:12:37 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#ifndef G4VFermiFragment_h
#define G4VFermiFragment_h 1

#include "G4FragmentVector.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

class G4VFermiFragment 
{
public:
  G4VFermiFragment(const G4int anA, const G4int aZ, const G4int Pol, const G4double ExE):
    A(anA),
    Z(aZ),
    Polarization(Pol),
    ExcitEnergy(ExE)
    {}

  virtual ~G4VFermiFragment() {};
  
protected:
  G4VFermiFragment() {};

private:

  G4VFermiFragment(const G4VFermiFragment &right);
  
  const G4VFermiFragment & operator=(const G4VFermiFragment &right);
  G4bool operator==(const G4VFermiFragment &right) const;
  G4bool operator!=(const G4VFermiFragment &right) const;
  
public:

  virtual G4FragmentVector * GetFragment(const G4LorentzVector & aMomentum) const = 0;

  G4int GetA(void) const 
  {
    return A;
  }
  
  G4int GetZ(void) const 
  {
    return Z;
  }
  
  G4int GetPolarization(void) const 
  {
    return Polarization;
  }

  G4double GetExcitationEnergy(void) const 
  {
    return ExcitEnergy;
  }

  G4double GetFragmentMass(void) const
  {
    return G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z,A);
  }

  G4double GetTotalEnergy(void) const
  {
    return this->GetFragmentMass() + this->GetExcitationEnergy();
  }

protected:

  G4int A;

  G4int Z;

  G4int Polarization;

  G4double ExcitEnergy;


};


#endif


