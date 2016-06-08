// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AlphaCoulombBarrier.cc,v 1.1 2000/06/09 11:43:34 larazb Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)


#include "G4AlphaCoulombBarrier.hh"

G4AlphaCoulombBarrier::G4AlphaCoulombBarrier(const G4AlphaCoulombBarrier & right)
{
  G4Exception("G4AlphaCoulombBarrier::copy_constructor meant to not be accessable.");
}


const G4AlphaCoulombBarrier & G4AlphaCoulombBarrier::operator=(const G4AlphaCoulombBarrier & right)
{
 G4Exception("G4AlphaCoulombBarrier::operator= meant to not be accessable.");
 return *this;
}

G4bool G4AlphaCoulombBarrier::operator==(const G4AlphaCoulombBarrier & right) const 
{
 return false;
}

G4bool G4AlphaCoulombBarrier::operator!=(const G4AlphaCoulombBarrier & right) const 
{
 return true;
}


G4double G4AlphaCoulombBarrier::BarrierPenetrationFactor(const G4double aZ) const
{
	// Data comes from 
	// Dostrovsky, Fraenkel and Friedlander
	// Physical Review, vol 116, num. 3 1959
	// 
	// const G4int size = 5;
	// const G4double Zlist[size]  = {10.0, 20.0, 30.0, 50.0, 70.0};
	// const G4double Kalpha[size] = {0.68, 0.82, 0.91, 0.97, 0.98};
	G4double K = 1.0;
	if (aZ>=70.0) {
		K = 0.98;
	} else {
		K = (((0.23684e-5*aZ) - 0.42143e-3)*aZ + 0.25222e-1)*aZ + 0.46699;
	}
	return K;
}
