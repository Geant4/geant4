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
// $Id: G4BREPSolidTorus.cc,v 1.7 2003/06/16 16:52:51 gunter Exp $
// GEANT4 tag $Name: geant4-07-01 $
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
  
  // Save constructor parameters
  constructorParams.origin       = origin;
  constructorParams.axis         = axis;
  constructorParams.direction    = direction;
  constructorParams.MinRadius    = MinRadius;
  constructorParams.MaxRadius    = MaxRadius;
  
  active = 1;
  Initialize();
}

G4BREPSolidTorus::~G4BREPSolidTorus()
{
}

// Streams solid contents to output stream.
std::ostream& G4BREPSolidTorus::StreamInfo(std::ostream& os) const
{
  G4BREPSolid::StreamInfo( os )
  << "\n origin:       " << constructorParams.origin
  << "\n axis:         " << constructorParams.axis
  << "\n direction:    " << constructorParams.direction
  << "\n MinRadius:    " << constructorParams.MinRadius
  << "\n MaxRadius:    " << constructorParams.MaxRadius
  << "\n-----------------------------------------------------------\n";

  return os;
}

