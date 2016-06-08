// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VFissionBarrier.hh,v 1.3 2000/06/09 11:36:54 larazb Exp $
// GEANT4 tag $Name: geant4-03-00 $
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
