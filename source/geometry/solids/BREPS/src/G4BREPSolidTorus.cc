//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4BREPSolidTorus.cc,v 1.11 2010-10-20 09:14:11 gcosmo Exp $
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
  // Save constructor parameters
  constructorParams.origin       = origin;
  constructorParams.axis         = axis;
  constructorParams.direction    = direction;
  constructorParams.MinRadius    = MinRadius;
  constructorParams.MaxRadius    = MaxRadius;
  
  active = 1;
  InitializeTorus();
}

G4BREPSolidTorus::G4BREPSolidTorus( __void__& a )
  : G4BREPSolid(a)
{
}

G4BREPSolidTorus::~G4BREPSolidTorus()
{
}

G4BREPSolidTorus::G4BREPSolidTorus(const G4BREPSolidTorus& rhs)
  : G4BREPSolid(rhs)
{
  constructorParams.origin    = rhs.constructorParams.origin;
  constructorParams.axis      = rhs.constructorParams.axis;
  constructorParams.direction = rhs.constructorParams.direction;
  constructorParams.MinRadius = rhs.constructorParams.MinRadius;
  constructorParams.MaxRadius = rhs.constructorParams.MaxRadius;
  
  InitializeTorus();
}

G4BREPSolidTorus&
G4BREPSolidTorus::operator = (const G4BREPSolidTorus& rhs) 
{
  // Check assignment to self
  //
  if (this == &rhs)  { return *this; }

  // Copy base class data
  //
  G4BREPSolid::operator=(rhs);

  // Copy data
  //
  constructorParams.origin    = rhs.constructorParams.origin;
  constructorParams.axis      = rhs.constructorParams.axis;
  constructorParams.direction = rhs.constructorParams.direction;
  constructorParams.MinRadius = rhs.constructorParams.MinRadius;
  constructorParams.MaxRadius = rhs.constructorParams.MaxRadius;
  
  InitializeTorus();

  return *this;
}  

void G4BREPSolidTorus::InitializeTorus()
{
  SurfaceVec    = new G4Surface*[1];
  SurfaceVec[0] = new G4ToroidalSurface( constructorParams.origin,
                                         constructorParams.axis,
                                         constructorParams.direction,
					 constructorParams.MinRadius,
                                         constructorParams.MaxRadius );
  nb_of_surfaces = 1;

  Initialize();
}

G4VSolid* G4BREPSolidTorus::Clone() const
{
  return new G4BREPSolidTorus(*this);
}

std::ostream& G4BREPSolidTorus::StreamInfo(std::ostream& os) const
{
  // Streams solid contents to output stream.

  G4BREPSolid::StreamInfo( os )
  << "\n origin:       " << constructorParams.origin
  << "\n axis:         " << constructorParams.axis
  << "\n direction:    " << constructorParams.direction
  << "\n MinRadius:    " << constructorParams.MinRadius
  << "\n MaxRadius:    " << constructorParams.MaxRadius
  << "\n-----------------------------------------------------------\n";

  return os;
}

