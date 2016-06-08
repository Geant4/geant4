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
// $Id: G4AlphaCoulombBarrier.hh,v 1.4 2001/08/01 17:04:07 hpw Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)

#ifndef G4AlphaCoulombBarrier_h
#define G4AlphaCoulombBarrier_h 1

#include "G4CoulombBarrier.hh"
#include "globals.hh"

class G4AlphaCoulombBarrier : public G4CoulombBarrier
{
public:
  G4AlphaCoulombBarrier() : G4CoulombBarrier(4,2) {};
  ~G4AlphaCoulombBarrier() {};

private:
  G4AlphaCoulombBarrier(const G4AlphaCoulombBarrier & right);

  const G4AlphaCoulombBarrier & operator=(const G4AlphaCoulombBarrier & right);
  G4bool operator==(const G4AlphaCoulombBarrier & right) const;
  G4bool operator!=(const G4AlphaCoulombBarrier & right) const;
  
private:

  virtual G4double BarrierPenetrationFactor(const G4double aZ) const;


};

#endif
