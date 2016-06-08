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
// $Id: G4Be11GEMCoulombBarrier.hh,v 1.1 2002/06/06 17:39:31 larazb Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)

#ifndef G4Be11GEMCoulombBarrier_h
#define G4Be11GEMCoulombBarrier_h 1

#include "G4GEMCoulombBarrierHE.hh"
#include "globals.hh"

class G4Be11GEMCoulombBarrier : public G4GEMCoulombBarrierHE
{
public:
  G4Be11GEMCoulombBarrier() : G4GEMCoulombBarrierHE(11,4) {};
  ~G4Be11GEMCoulombBarrier() {};

private:
  G4Be11GEMCoulombBarrier(const G4Be11GEMCoulombBarrier & right);

  const G4Be11GEMCoulombBarrier & operator=(const G4Be11GEMCoulombBarrier & right);
  G4bool operator==(const G4Be11GEMCoulombBarrier & right) const;
  G4bool operator!=(const G4Be11GEMCoulombBarrier & right) const;
  

};

#endif
