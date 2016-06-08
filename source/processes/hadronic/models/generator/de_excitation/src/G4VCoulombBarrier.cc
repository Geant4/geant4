// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VCoulombBarrier.cc,v 1.1 2000/06/09 11:43:37 larazb Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)


#include "G4VCoulombBarrier.hh"
#include "g4std/strstream"

G4VCoulombBarrier::G4VCoulombBarrier(const G4int anA, const G4int aZ)
{
	if (anA >= aZ && anA > 0) {
		theA = anA;
		theZ = aZ;
	} else {
		char errMessage[1024];
		G4std::ostrstream errOs(errMessage,1024);
		errOs << "G4VCoulombBarrier::G4VCoulombBarrier: ";
		errOs << "Wrong values for ";
		errOs << "A = " << anA << " ";
		errOs << "and Z = " << aZ << G4endl;
		G4Exception(errMessage);
	}
}


G4VCoulombBarrier::G4VCoulombBarrier(const G4VCoulombBarrier & right)
{
  G4Exception("G4VCoulombBarrier::copy_constructor meant to not be accessable.");
}


const G4VCoulombBarrier & G4VCoulombBarrier::operator=(const G4VCoulombBarrier & right)
{
 G4Exception("G4VCoulombBarrier::operator= meant to not be accessable.");
 return *this;
}

G4bool G4VCoulombBarrier::operator==(const G4VCoulombBarrier & right) const 
{
 return false;
}

G4bool G4VCoulombBarrier::operator!=(const G4VCoulombBarrier & right) const 
{
 return true;
}

