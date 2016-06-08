// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4He3CoulombBarrier.cc,v 1.1 2000/06/09 11:43:36 larazb Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)


#include "G4He3CoulombBarrier.hh"

G4He3CoulombBarrier::G4He3CoulombBarrier(const G4He3CoulombBarrier & right)
{
  G4Exception("G4He3CoulombBarrier::copy_constructor meant to not be accessable.");
}


const G4He3CoulombBarrier & G4He3CoulombBarrier::operator=(const G4He3CoulombBarrier & right)
{
 G4Exception("G4He3CoulombBarrier::operator= meant to not be accessable.");
 return *this;
}

G4bool G4He3CoulombBarrier::operator==(const G4He3CoulombBarrier & right) const 
{
 return false;
}

G4bool G4He3CoulombBarrier::operator!=(const G4He3CoulombBarrier & right) const 
{
 return true;
}


G4double G4He3CoulombBarrier::BarrierPenetrationFactor(const G4double aZ) const
{
	// Data comes from 
	// Dostrovsky, Fraenkel and Friedlander
	// Physical Review, vol 116, num. 3 1959
	// 
	// const G4int size = 5;
	// const G4double Zlist[size]  = {10.0, 20.0, 30.0, 50.0, 70.0};
	// const G4double KHe3[size] = {0.68, 0.82, 0.91, 0.97, 0.98};
	//
	// K for He3 is K for alphas + 0.12
	G4double K = 1.0;
	if (aZ>=70.0) {
		K = 0.98;
	} else {
		K = (((0.23684e-5*aZ) - 0.42143e-3)*aZ + 0.25222e-1)*aZ + 0.46699;
	}
	return K+0.12;
}
