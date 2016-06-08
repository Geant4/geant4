// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4CoulombBarrier.hh,v 1.1 2000/06/09 11:36:52 larazb Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)

#ifndef G4CoulombBarrier_h
#define G4CoulombBarrier_h 1

#include "G4VCoulombBarrier.hh"
#include "globals.hh"


class G4CoulombBarrier : public G4VCoulombBarrier
{
public:
	G4CoulombBarrier() : G4VCoulombBarrier(1,0) {};
	G4CoulombBarrier(const G4int anA,const G4int aZ) :
									G4VCoulombBarrier(anA,aZ) {};
	~G4CoulombBarrier() {};

private:
	G4CoulombBarrier(const G4CoulombBarrier & right);

	const G4CoulombBarrier & operator=(const G4CoulombBarrier & right);
	G4bool operator==(const G4CoulombBarrier & right) const;
	G4bool operator!=(const G4CoulombBarrier & right) const;
  
public:
	G4double GetCoulombBarrier(const G4int ARes, const G4int ZRes, 
									const G4double U) const;


private:

	virtual G4double BarrierPenetrationFactor(const G4double aZ) const {return 1.0;}


};

#endif
