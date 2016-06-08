// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4He3CoulombBarrier.hh,v 1.1 2000/06/09 11:36:53 larazb Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)

#ifndef G4He3CoulombBarrier_h
#define G4He3CoulombBarrier_h 1

#include "G4CoulombBarrier.hh"
#include "globals.hh"

class G4He3CoulombBarrier : public G4CoulombBarrier
{
public:
	G4He3CoulombBarrier() : G4CoulombBarrier(3,2) {};
	~G4He3CoulombBarrier() {};

private:
	G4He3CoulombBarrier(const G4He3CoulombBarrier & right);

	const G4He3CoulombBarrier & operator=(const G4He3CoulombBarrier & right);
	G4bool operator==(const G4He3CoulombBarrier & right) const;
	G4bool operator!=(const G4He3CoulombBarrier & right) const;
  
private:

	virtual G4double BarrierPenetrationFactor(const G4double aZ) const;


};

#endif
