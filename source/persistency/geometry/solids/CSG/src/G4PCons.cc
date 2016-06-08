//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PCons.cc,v 1.3 2001/07/11 10:02:24 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// class G4PCons
//
// Implementation for G4PCons class
//
// History:
// 19.06.98 A.Kimura Converted G4Cons.cc

#include "G4VSolid.hh"
#include "G4PCons.hh"
#include "G4Cons.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

G4PCons::~G4PCons() 
{;}

// make a transient object
G4VSolid* G4PCons::MakeTransientObject() const {
    G4VSolid* transientObject = new G4Cons(GetName(),
					 fRmin1, fRmax1,
					 fRmin2, fRmax2,
					 fDz,
					 fSPhi, fDPhi);
    return transientObject;
}

// constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//               - note if pDPhi>2PI then reset to 2PI

G4PCons::G4PCons(const G4Cons* theCons) : G4PCSGSolid(theCons->GetName()) {

    G4double pRmin1 = theCons->GetInnerRadiusMinusZ();
    G4double pRmax1 = theCons->GetOuterRadiusMinusZ();
    G4double pRmin2 = theCons->GetInnerRadiusPlusZ();
    G4double pRmax2 = theCons->GetOuterRadiusPlusZ();
    G4double pDz    = theCons->GetZHalfLength();
    G4double pSPhi  = theCons->GetStartPhiAngle();
    G4double pDPhi  = theCons->GetDeltaPhiAngle();

// Check z-len
    if (pDz>0) {
	fDz=pDz;
    } else {
	G4Exception("Error in G4PCons::G4PCons - invalid z half-length");
    }

// Check radii
    if (pRmin1<pRmax1 && pRmin2<pRmax2 && pRmin1>=0 && pRmin2>=0) {
	fRmin1=pRmin1; fRmax1=pRmax1;
	fRmin2=pRmin2; fRmax2=pRmax2;
    } else {
	G4Exception("Error in G4PCons::G4PCons - invalid radii");
    }

// Check angles
    if (pDPhi>=2.0*M_PI) {
	fDPhi=2*M_PI;
	fSPhi=0;
    } else {
	if (pDPhi>0) {
	    fDPhi=pDPhi;
	} else {
	    G4Exception("Error in G4PCons::G4PCons - invalid pDPhi");
	}
	
// Ensure pSPhi in 0-2PI or -2PI-0 range if shape crosses 0
	if (pSPhi<0) {
	    fSPhi=2.0*M_PI-fmod(fabs(pSPhi),2.0*M_PI);
	} else {
	    fSPhi=fmod(pSPhi,2.0*M_PI);
	}
	    
	if (fSPhi+fDPhi>2.0*M_PI) fSPhi-=2.0*M_PI;
    }
}

//  ******************************* End of G4PCons.cc file **********************************
