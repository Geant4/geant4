// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TritonCoulombBarrier.hh,v 1.1 2000/06/09 11:36:53 larazb Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)

#ifndef G4TritonCoulombBarrier_h
#define G4TritonCoulombBarrier_h 1

#include "G4CoulombBarrier.hh"
#include "globals.hh"

class G4TritonCoulombBarrier : public G4CoulombBarrier
{
public:
	G4TritonCoulombBarrier() : G4CoulombBarrier(3,1) {};
	~G4TritonCoulombBarrier() {};

private:
	G4TritonCoulombBarrier(const G4TritonCoulombBarrier & right);

	const G4TritonCoulombBarrier & operator=(const G4TritonCoulombBarrier & right);
	G4bool operator==(const G4TritonCoulombBarrier & right) const;
	G4bool operator!=(const G4TritonCoulombBarrier & right) const;
  
private:

	virtual G4double BarrierPenetrationFactor(const G4double aZ) const;


};

#endif
