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
// $Id: G4GEMCoulombBarrierHE.hh,v 1.2 2002/12/12 19:17:06 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)

#ifndef G4GEMCoulombBarrierHE_h
#define G4GEMCoulombBarrierHE_h 1

#include "G4VCoulombBarrier.hh"
#include "globals.hh"


class G4GEMCoulombBarrierHE : public G4VCoulombBarrier
{
public:
  G4GEMCoulombBarrierHE() : G4VCoulombBarrier(1,0) {};
  G4GEMCoulombBarrierHE(const G4int anA,const G4int aZ) :
  G4VCoulombBarrier(anA,aZ) {};
  ~G4GEMCoulombBarrierHE() {};

private:
  G4GEMCoulombBarrierHE(const G4GEMCoulombBarrierHE & right);

  const G4GEMCoulombBarrierHE & operator=(const G4GEMCoulombBarrierHE & right);
  G4bool operator==(const G4GEMCoulombBarrierHE & right) const;
  G4bool operator!=(const G4GEMCoulombBarrierHE & right) const;
  
public:
  G4double GetCoulombBarrier(const G4int ARes, const G4int ZRes, 
			     const G4double U) const;


private:

  G4double BarrierPenetrationFactor(const G4double aZ) const;

  virtual G4double CalcCompoundRadius(const G4double ARes) const;

};
#endif

