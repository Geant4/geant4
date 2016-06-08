// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PTorus.cc,v 1.3 1999/12/15 14:51:25 gunter Exp $
// GEANT4 tag $Name: geant4-03-01 $
//
// 
// class G4PTorus
//
// Implementation
//
// History:
// 19.06.98 A.Kimura Converted G4Torus.cc

#include "G4VSolid.hh"
#include "G4PTorus.hh"
#include "G4Torus.hh"

// Constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pdphi>2PI then reset to 2PI
G4PTorus::G4PTorus(const G4Torus* theTorus)
 : G4PCSGSolid(theTorus->GetName())
{
    SetAllParameters(theTorus->GetRmin(),
		     theTorus->GetRmax(),
		     theTorus->GetRtor(),
		     theTorus->GetSPhi(),
		     theTorus->GetDPhi());
}

void
G4PTorus::SetAllParameters(
	       G4double pRmin,
	       G4double pRmax,
	       G4double pRtor,
	       G4double pSPhi,
	       G4double pDPhi)
{
// Check swept radius
    if (pRtor>=pRmax)
	{
	   fRtor=pRtor;
	}
    else
	{
	    G4Exception("Error in G4PTorus::SetAllParameters - invalid swept radius");
	}

// Check radii
    if (pRmin<pRmax&&pRmin>=0)
	{
	   fRmin=pRmin; fRmax=pRmax;
	}
    else
	{
	    G4Exception("Error in G4PTorus::SetAllParameters - invalid radii");
	}

// Check angles
    if (pDPhi>=2.0*M_PI)
	{
	    fDPhi=2*M_PI;
	}
    else
	{
	    if (pDPhi>0)
		{
		    fDPhi = pDPhi;
		}
	    else
		{
		    G4Exception("Error in G4PTorus::SetAllParameters - invalid dphi");
		}
	}
	
// Ensure psphi in 0-2PI or -2PI-0 range if shape crosses 0
    fSPhi = pSPhi;

    if (fSPhi<0)
	{
	    fSPhi=2.0*M_PI-fmod(fabs(fSPhi),2.0*M_PI);
	}
    else
	{
	    fSPhi=fmod(fSPhi,2.0*M_PI);
	}

    if (fSPhi+fDPhi>2.0*M_PI)
	{
	    fSPhi-=2.0*M_PI;
	}
}

// Destructor
G4PTorus::~G4PTorus()
{;}

// --------------------------------------------------------------------------

// make a transient object
G4VSolid* G4PTorus::MakeTransientObject() const {
    G4VSolid* transientObject = new G4Torus(GetName(),
					   fRmin, fRmax,
					   fRtor,
					   fSPhi, fDPhi);
    return transientObject;
}
