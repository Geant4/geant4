// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PSphere.cc,v 1.3 1999/11/17 10:49:03 morita Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// 
// class G4PSphere
//
// Implementation for G4PSphere class
//
// History:
// 19.06.98 A.Kimura Converted G4Sphere.cc

#include "G4VSolid.hh"
#include "G4PSphere.hh"
#include "G4Sphere.hh"

// Destructor
G4PSphere::~G4PSphere()
{;}

// make a transient object
G4VSolid* G4PSphere::MakeTransientObject() const {
    G4VSolid* transientObject = new G4Sphere(GetName(),
					     fRmin,
					     fRmax,
					     fSPhi,
					     fDPhi,
					     fSTheta,
					     fDTheta);
    return transientObject;
}

// constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pDPhi>2PI then reset to 2PI

G4PSphere::G4PSphere(const G4Sphere* theSphere)
    : G4PCSGSolid(theSphere->GetName())
{
    G4double pRmin = fRmin;
    G4double pRmax = fRmax;
    G4double pSPhi = fSPhi;
    G4double pDPhi = fDPhi;
    G4double pSTheta = fSTheta;
    G4double pDTheta = fDTheta;

    // Check radii
    if (pRmin<pRmax&&pRmin>=0)
	{
	   fRmin=pRmin; fRmax=pRmax;
	}
    else
	{
	    G4Exception("Error in G4PSphere::G4PSphere - invalid radii");
	}

// Check angles
    if (pDPhi>=2.0*M_PI)
	{
	   fDPhi=2*M_PI;
	}
    else if (pDPhi>0)
	{
	   fDPhi=pDPhi;
	}
    else
	{
	    G4Exception("Error in G4PSphere::G4PSphere - invalid DPhi");
	}
// Convert fSPhi to 0-2PI
    if (pSPhi<0)
	{
	   fSPhi=2.0*M_PI-fmod(fabs(pSPhi),2.0*M_PI);
	}
    else
	{
	   fSPhi=fmod(pSPhi,2.0*M_PI);
	}
// Sphere is placed such that fSPhi+fDPhi>2.0*M_PI ! fSPhi could be < 0 !!? P
    if (fSPhi+fDPhi>2.0*M_PI) fSPhi-=2.0*M_PI;

// Check theta angles
    if (pSTheta<0 || pSTheta>M_PI)
	{
	    G4Exception("Error in G4PSphere::G4PSphere - stheta outside 0-PI range");
	}
    else
	{
	   fSTheta=pSTheta;
	}

    if (pDTheta+pSTheta>=M_PI)
	{
	   fDTheta=M_PI-pSTheta;
	}
    else if (pDTheta>0)
	{
	   fDTheta=pDTheta;
	}
    else
	{
	    G4Exception("Error in G4PSphere::G4PSphere - invalid pDTheta");
	}

}

// ******************************  End of G4PSphere.cc  ****************************************
