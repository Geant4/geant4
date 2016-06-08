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
// $Id: G4PTrd.cc,v 1.4 2001/07/11 10:02:24 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $
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

#include <math.h>


// Constructor - check & set half widths

G4PTrd::G4PTrd(const G4Trd* theTrd)
 : G4PCSGSolid(theTrd->GetName())
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
				 G4double pdz)
{
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
          fDx1=G4std::max(pdx1,Minimum_length); 
          fDx2=G4std::max(pdx2,Minimum_length); 
          fDy1=G4std::max(pdy1,Minimum_length); 
          fDy2=G4std::max(pdy2,Minimum_length); 
          fDz=G4std::max(pdz,Minimum_length);
        }
      else
        G4Exception("Error in G4PTrd::G4PTrd - One or more parameters are < 0");
    }
}

// -------------------------------------------------------------

// Destructor

G4PTrd::~G4PTrd()
{;}

// -----------------------------------------------------------------------

// make a transient object
G4VSolid* G4PTrd::MakeTransientObject() const
{
    G4VSolid* transientObject = new G4Trd(GetName(),
				       fDx1, fDx2,
				       fDy1, fDy2,
				       fDz);
    return transientObject;
}

