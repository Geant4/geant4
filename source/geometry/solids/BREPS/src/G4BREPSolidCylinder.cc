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
// $Id$
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4BREPSolidCylinder.cc
//
// ----------------------------------------------------------------------

#include "G4BREPSolidCylinder.hh"
#include "G4CircularCurve.hh"
#include "G4FPlane.hh"
#include "G4FCylindricalSurface.hh"

G4BREPSolidCylinder::G4BREPSolidCylinder(const G4String& name,
					 const G4ThreeVector& origin,
					 const G4ThreeVector& axis,
					 const G4ThreeVector& direction,
					       G4double radius,
					       G4double length)
  : G4BREPSolid(name)
{
  nb_of_surfaces = 3;
  active=1;
  
  // Save constructor parameters
  constructorParams.origin       = origin;
  constructorParams.axis         = axis;
  constructorParams.direction    = direction;
  constructorParams.length       = length;
  constructorParams.radius       = radius;
  
  InitializeCylinder();
}

G4BREPSolidCylinder::G4BREPSolidCylinder( __void__& a )
  : G4BREPSolid(a)
{
}

G4BREPSolidCylinder::~G4BREPSolidCylinder()
{
}

G4BREPSolidCylinder::G4BREPSolidCylinder(const G4BREPSolidCylinder& rhs)
  : G4BREPSolid(rhs)
{
  constructorParams.origin       = rhs.constructorParams.origin;
  constructorParams.axis         = rhs.constructorParams.axis;
  constructorParams.direction    = rhs.constructorParams.direction;
  constructorParams.length       = rhs.constructorParams.length;
  constructorParams.radius       = rhs.constructorParams.radius;
  
  InitializeCylinder();
}

G4BREPSolidCylinder&
G4BREPSolidCylinder::operator = (const G4BREPSolidCylinder& rhs) 
{
  // Check assignment to self
  //
  if (this == &rhs)  { return *this; }

  // Copy base class data
  //
  G4BREPSolid::operator=(rhs);

  // Copy data
  //
  constructorParams.origin       = rhs.constructorParams.origin;
  constructorParams.axis         = rhs.constructorParams.axis;
  constructorParams.direction    = rhs.constructorParams.direction;
  constructorParams.length       = rhs.constructorParams.length;
  constructorParams.radius       = rhs.constructorParams.radius;
  
  InitializeCylinder();

  return *this;
}  

void G4BREPSolidCylinder::InitializeCylinder()
{
  SurfaceVec = new G4Surface*[3];
  G4CurveVector cv;
  G4CircularCurve* tmp;

  // Creation of the cylindrical surface
  SurfaceVec[0] = new G4FCylindricalSurface(constructorParams.origin,
                                            constructorParams.axis,
                                            constructorParams.radius,
                                            constructorParams.length);
  //SurfaceVec[0]->SetBoundaries(&cv);
  //cv.clear();

  // Creation of the first circular surface, which origin is origin
  G4Point3D  ArcStart1 = G4Point3D( constructorParams.origin
                                + ( constructorParams.radius
                                *   constructorParams.direction ) );
  G4Vector3D axis1 = G4Vector3D( constructorParams.axis.cross( constructorParams.direction ) );

  tmp = new G4CircularCurve;
  tmp->Init( G4Axis2Placement3D(constructorParams.direction, axis1,
                                constructorParams.origin),
             constructorParams.radius );
  tmp->SetBounds(ArcStart1, ArcStart1);
  cv.push_back(tmp);

  SurfaceVec[1] = new G4FPlane(constructorParams.direction, axis1,
                               constructorParams.origin);
  SurfaceVec[1]->SetBoundaries(&cv);
  cv.clear();
  

  // Creation of the second circular surface
  G4Point3D  origin2   = G4Point3D( constructorParams.origin
                                + ( constructorParams.length
                                *   constructorParams.axis ) );  
  G4Point3D  ArcStart2 = origin2
                       + G4Point3D( constructorParams.radius
                                  * constructorParams.direction );
  G4Vector3D axis2     = axis1;

  tmp = new G4CircularCurve;
  tmp->Init( G4Axis2Placement3D(constructorParams.direction,
                                axis2, origin2),
             constructorParams.radius);
  tmp->SetBounds(ArcStart2, ArcStart2);
  cv.push_back(tmp);

  SurfaceVec[2] = new G4FPlane(constructorParams.direction, axis2, origin2);
  SurfaceVec[2]->SetBoundaries(&cv);
  cv.clear();

  Initialize();
}

G4VSolid* G4BREPSolidCylinder::Clone() const
{
  return new G4BREPSolidCylinder(*this);
}

std::ostream& G4BREPSolidCylinder::StreamInfo(std::ostream& os) const
{
  // Streams solid contents to output stream.

  G4BREPSolid::StreamInfo( os )
  << "\n origin:       " << constructorParams.origin
  << "\n axis:         " << constructorParams.axis
  << "\n direction:    " << constructorParams.direction
  << "\n length:       " << constructorParams.length
  << "\n radius:       " << constructorParams.radius
  << "\n-----------------------------------------------------------\n";

  return os;
}

