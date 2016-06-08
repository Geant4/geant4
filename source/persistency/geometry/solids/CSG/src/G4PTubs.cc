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
// $Id: G4PTubs.cc,v 1.2.8.1 2001/06/28 19:11:32 gunter Exp $
// GEANT4 tag $Name:  $
//
// 
// class G4PTubs
//
// History:
// 19.06.98 A.Kimura Converted G4Tubs.cc


#include "G4VSolid.hh"
#include "G4PTubs.hh"
#include "G4Tubs.hh"

// Constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pdphi>2PI then reset to 2PI
G4PTubs::G4PTubs(const G4Tubs* theTubs)
 : G4PCSGSolid(theTubs->GetName())
{

    G4double pRMin = theTubs->GetInnerRadius();
    G4double pRMax = theTubs->GetOuterRadius();
    G4double pDz = theTubs->GetZHalfLength();
    G4double pSPhi = theTubs->GetStartPhiAngle();
    G4double pDPhi = theTubs->GetDeltaPhiAngle();

// Check z-len
    if (pDz>0)
	{
	    fDz=pDz;
	}
    else
	{
	    G4Exception("Error in G4PTubs::G4PTubs - invalid z half-length");
	}

// Check radii
    if (pRMin<pRMax&&pRMin>=0)
	{
	    fRMin=pRMin; fRMax=pRMax;
	}
    else
	{
	    G4Exception("Error in G4PTubs::G4PTubs - invalid radii");
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
		    G4Exception("Error in G4PTubs::G4PTubs - invalid dphi");
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
G4PTubs::~G4PTubs()
{;}

// make a transient object
G4VSolid* G4PTubs::MakeTransientObject() const
{
    G4VSolid* transientObject = new G4Tubs(GetName(),
					 fRMin, fRMax,
					 fDz,
					 fSPhi, fDPhi);
    return transientObject;
}
