// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DeuteronCoulombBarrier.hh,v 1.1 2000/06/09 11:36:52 larazb Exp $
// GEANT4 tag $Name: geant4-03-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)

#ifndef G4DeuteronCoulombBarrier_h
#define G4DeuteronCoulombBarrier_h 1

#include "G4CoulombBarrier.hh"
#include "globals.hh"

class G4DeuteronCoulombBarrier : public G4CoulombBarrier
{
public:
	G4DeuteronCoulombBarrier() : G4CoulombBarrier(2,1) {};
	~G4DeuteronCoulombBarrier() {};

private:
	G4DeuteronCoulombBarrier(const G4DeuteronCoulombBarrier & right);

	const G4DeuteronCoulombBarrier & operator=(const G4DeuteronCoulombBarrier & right);
	G4bool operator==(const G4DeuteronCoulombBarrier & right) const;
	G4bool operator!=(const G4DeuteronCoulombBarrier & right) const;
  
private:

	virtual G4double BarrierPenetrationFactor(const G4double aZ) const;


};

#endif
