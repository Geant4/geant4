// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronCoulombBarrier.hh,v 1.1 2000/06/09 11:36:53 larazb Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)

#ifndef G4NeutronCoulombBarrier_h
#define G4NeutronCoulombBarrier_h 1

#include "G4CoulombBarrier.hh"
#include "globals.hh"

class G4NeutronCoulombBarrier : public G4CoulombBarrier
{
public:
	G4NeutronCoulombBarrier() : G4CoulombBarrier(1,0) {};
	~G4NeutronCoulombBarrier() {};

private:
	G4NeutronCoulombBarrier(const G4NeutronCoulombBarrier & right);

	const G4NeutronCoulombBarrier & operator=(const G4NeutronCoulombBarrier & right);
	G4bool operator==(const G4NeutronCoulombBarrier & right) const;
	G4bool operator!=(const G4NeutronCoulombBarrier & right) const;
  
private:

	virtual G4double BarrierPenetrationFactor(const G4double aZ) const
	{ return 1.0;}


};

#endif
