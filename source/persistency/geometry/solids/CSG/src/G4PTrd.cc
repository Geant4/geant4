// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PTrd.cc,v 2.0 1998/07/02 16:13:29 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
//
// Implementation for G4Trd class
//
//
// History:
// 19.06.98 A.Kimura Converted G4Trd.cc


#include "G4VSolid.hh"
#include "G4PTrd.hh"
#include "G4Trd.hh"

#include "G4AffineTransform.hh"

#include <math.h>


// Constructor - check & set half widths

G4PTrd::G4PTrd(const G4Trd* theTrd) : G4PCSGSolid(theTrd->GetName())
{
    G4double pdx1 = theTrd->GetXHalfLength1();
    G4double pdx2 = theTrd->GetXHalfLength2(); 
    G4double pdy1 = theTrd->GetYHalfLength1();
    G4double pdy2 = theTrd->GetYHalfLength2();
    G4double pdz  = theTrd->GetZHalfLength();
    CheckAndSetAllParameters (pdx1, pdx2, pdy1, pdy2, pdz);
}

void
G4PTrd::CheckAndSetAllParameters (G4double pdx1,  G4double pdx2,
				 G4double pdy1,  G4double pdy2,
				 G4double pdz) {
  if (pdx1>0&&pdx2>0&&pdy1>0&&pdy2>0&&pdz>0)
    {
      fDx1=pdx1; fDx2=pdx2;
      fDy1=pdy1; fDy2=pdy2;
      fDz=pdz;
    }
  else
    {
      if (pdx1>=0 && pdx2>=0 && pdy1>=0 && pdy2>=0 && pdz>=0)
        {
          // G4double  Minimum_length= (1+per_thousand) * kCarTolerance/2.;
          //  FIX-ME : temporary solution for ZERO or very-small parameters.
          G4double  Minimum_length= kCarTolerance/2.;
          fDx1=max(pdx1,Minimum_length); 
          fDx2=max(pdx2,Minimum_length); 
          fDy1=max(pdy1,Minimum_length); 
          fDy2=max(pdy2,Minimum_length); 
          fDz=max(pdz,Minimum_length);
        }
      else
        G4Exception("Error in G4PTrd::G4PTrd - One or more parameters are < 0");
    }
}

// -------------------------------------------------------------

// Destructor

G4PTrd::~G4PTrd()
{}

// -----------------------------------------------------------------------

// make a transient object
G4VSolid* G4PTrd::MakeTransientObject() const {
    G4VSolid* transientObject = new G4Trd(GetName(),
				       fDx1, fDx2,
				       fDy1, fDy2,
				       fDz);
    return transientObject;
}

