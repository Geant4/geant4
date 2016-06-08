// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AlphaCoulombBarrier.hh,v 1.1 2000/06/09 11:36:49 larazb Exp $
// GEANT4 tag $Name: geant4-03-01 $
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
