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
// $Id: G4Li9GEMCoulombBarrier.hh,v 1.1 2003/08/26 18:42:23 lara Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)

#ifndef G4Li9GEMCoulombBarrier_h
#define G4Li9GEMCoulombBarrier_h 1

#include "G4GEMCoulombBarrierHE.hh"
#include "globals.hh"

class G4Li9GEMCoulombBarrier : public G4GEMCoulombBarrierHE
{
public:
  G4Li9GEMCoulombBarrier() : G4GEMCoulombBarrierHE(9,3) {};
  ~G4Li9GEMCoulombBarrier() {};

private:
  G4Li9GEMCoulombBarrier(const G4Li9GEMCoulombBarrier & right);

  const G4Li9GEMCoulombBarrier & operator=(const G4Li9GEMCoulombBarrier & right);
  G4bool operator==(const G4Li9GEMCoulombBarrier & right) const;
  G4bool operator!=(const G4Li9GEMCoulombBarrier & right) const;
  

};

#endif
