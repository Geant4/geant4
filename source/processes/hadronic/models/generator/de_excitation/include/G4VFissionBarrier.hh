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
// $Id: G4VFissionBarrier.hh,v 1.3.8.1 2001/06/28 19:13:08 gunter Exp $
// GEANT4 tag $Name:  $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)

#ifndef G4VFissionBarrier_h
#define G4VFissionBarrier_h 1

#include "globals.hh"


class G4VFissionBarrier
{
public:
  G4VFissionBarrier() {};
  virtual ~G4VFissionBarrier() {};

private:
  G4VFissionBarrier(const G4VFissionBarrier & right);

  const G4VFissionBarrier & operator=(const G4VFissionBarrier & right);
  G4bool operator==(const G4VFissionBarrier & right) const;
  G4bool operator!=(const G4VFissionBarrier & right) const;
  
public:
  virtual G4double FissionBarrier(const G4int A, const G4int Z,const G4double U) = 0;
  

};

#endif
