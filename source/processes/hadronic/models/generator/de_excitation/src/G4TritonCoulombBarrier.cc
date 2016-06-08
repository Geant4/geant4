// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TritonCoulombBarrier.cc,v 1.1 2000/06/09 11:43:36 larazb Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)


#include "G4TritonCoulombBarrier.hh"

G4TritonCoulombBarrier::G4TritonCoulombBarrier(const G4TritonCoulombBarrier & right)
{
  G4Exception("G4TritonCoulombBarrier::copy_constructor meant to not be accessable.");
}


const G4TritonCoulombBarrier & G4TritonCoulombBarrier::operator=(const G4TritonCoulombBarrier & right)
{
 G4Exception("G4TritonCoulombBarrier::operator= meant to not be accessable.");
 return *this;
}

G4bool G4TritonCoulombBarrier::operator==(const G4TritonCoulombBarrier & right) const 
{
 return false;
}

G4bool G4TritonCoulombBarrier::operator!=(const G4TritonCoulombBarrier & right) const 
{
 return true;
}


G4double G4TritonCoulombBarrier::BarrierPenetrationFactor(const G4double aZ) const
{
	// Data comes from 
	// Dostrovsky, Fraenkel and Friedlander
	// Physical Review, vol 116, num. 3 1959
	// 
	// const G4int size = 5;
	// const G4double Zlist[size] = {10.0, 20.0, 30.0, 50.0, 70.0};
	// const G4double Kprot[size] = {0.42, 0.58, 0.68, 0.77, 0.80};
	// 
	// K for Triton is K for protons + 0.12
	G4double K = 1.0;
	if (aZ>=70.0) {
		K = 0.80;
	} else {
		K = (((0.2357e-5*aZ) - 0.42679e-3)*aZ + 0.27035e-1)*aZ + 0.19025;
	}
	return K+0.12;
}
