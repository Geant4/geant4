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
// $Id: G4BREPSolidTorus.cc,v 1.5 2001-07-11 09:59:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4BREPSolidTorus.cc
//
// ----------------------------------------------------------------------

#include "G4BREPSolidTorus.hh"
#include "G4ToroidalSurface.hh"

G4BREPSolidTorus::G4BREPSolidTorus(const    G4String& name,
				   const    G4ThreeVector& origin,
				   const    G4ThreeVector& axis,
				   const    G4ThreeVector& direction,
				   G4double MinRadius,
				   G4double MaxRadius)
  : G4BREPSolid(name)
{
  SurfaceVec = new G4Surface*[1];
  SurfaceVec[0] = new G4ToroidalSurface( origin, axis, direction,
					 MinRadius, MaxRadius);
  nb_of_surfaces = 1;
  active = 1;
  Initialize();
}

G4BREPSolidTorus::~G4BREPSolidTorus()
{
}
