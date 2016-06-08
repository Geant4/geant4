// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4CoulombBarrier.cc,v 1.1 2000/06/09 11:43:35 larazb Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)


#include "G4CoulombBarrier.hh"
#include "g4std/strstream"

G4CoulombBarrier::G4CoulombBarrier(const G4CoulombBarrier & right)
{
  G4Exception("G4CoulombBarrier::copy_constructor meant to not be accessable.");
}


const G4CoulombBarrier & G4CoulombBarrier::operator=(const G4CoulombBarrier & right)
{
 G4Exception("G4CoulombBarrier::operator= meant to not be accessable.");
 return *this;
}

G4bool G4CoulombBarrier::operator==(const G4CoulombBarrier & right) const 
{
 return false;
}

G4bool G4CoulombBarrier::operator!=(const G4CoulombBarrier & right) const 
{
 return true;
}



G4double G4CoulombBarrier::GetCoulombBarrier(const G4int ARes, const G4int ZRes, const G4double U) const 
  // Calculation of Coulomb potential energy (barrier) for outgoing fragment
{
	G4double Barrier = 0.0;
	if (ZRes > ARes || ARes < 1) {
		char errMessage[1024];
		G4std::ostrstream errOs(errMessage,1024);
		errOs << "G4CoulombBarrier::GetCoulombBarrier: ";
		errOs << "Wrong values for ";
		errOs << "residual nucleus A = " << ARes << " ";
		errOs << "and residual nucleus Z = " << ZRes << G4endl;
		G4Exception(errMessage);
	}
	if (GetA() == 1 && GetZ() == 0) {
		Barrier = 0.0;   // Neutron Coulomb Barrier is 0
	} else {
		G4double CompoundRadius = 2.173*fermi*(1.0+0.006103*G4double(GetZ())*G4double(ZRes))/
										(1.0+0.009443*G4double(GetZ())*G4double(ZRes));
		Barrier = elm_coupling/CompoundRadius * G4double(GetZ())*G4double(ZRes)/
							(pow(G4double(GetA()),1./3.) + pow(G4double(ARes),1./3.));

		// Barrier penetration coeficient
		G4double K = BarrierPenetrationFactor(ZRes);
		
		Barrier *= K;
		
		Barrier /= (1.0 + sqrt(U/(2.0*G4double(ARes))));
	}
	return Barrier;
}



