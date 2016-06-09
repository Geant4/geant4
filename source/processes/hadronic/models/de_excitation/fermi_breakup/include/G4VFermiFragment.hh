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
// $Id: G4VFermiFragment.hh,v 1.1 2003/08/26 18:34:20 lara Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
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


