// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VCoulombBarrier.hh,v 1.1 2000/06/09 11:36:54 larazb Exp $
// GEANT4 tag $Name: geant4-03-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)

#ifndef G4VCoulombBarrier_h
#define G4VCoulombBarrier_h 1

#include "globals.hh"


class G4VCoulombBarrier
{
public:
	G4VCoulombBarrier(const G4int anA, const G4int aZ);
	virtual ~G4VCoulombBarrier() {};

protected:
	G4VCoulombBarrier() : theA(1),theZ(0) {};

private:
	G4VCoulombBarrier(const G4VCoulombBarrier & right);

	const G4VCoulombBarrier & operator=(const G4VCoulombBarrier & right);
	G4bool operator==(const G4VCoulombBarrier & right) const;
	G4bool operator!=(const G4VCoulombBarrier & right) const;
  
public:
	virtual G4double GetCoulombBarrier(const G4int ARes, const G4int ZRes, 
												const G4double U) const = 0;
					
	G4int GetA(void) const {return theA;}
	G4int GetZ(void) const {return theZ;}
					
private: 
	
		G4int theA;
		G4int theZ;

};

#endif
