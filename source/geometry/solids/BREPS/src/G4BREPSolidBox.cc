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
// $Id: G4BREPSolidBox.cc,v 1.14 2010-10-20 09:14:11 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4BREPSolidBox.cc
//
// ----------------------------------------------------------------------

#include "G4BREPSolidBox.hh"
#include "G4FPlane.hh"
#include "G4Point3DVector.hh"

G4BREPSolidBox::G4BREPSolidBox(const G4String& name,
			       const G4Point3D& Pt1,
			       const G4Point3D& Pt2,
			       const G4Point3D& Pt3,
			       const G4Point3D& Pt4,
			       const G4Point3D& Pt5,
			       const G4Point3D& Pt6,
			       const G4Point3D& Pt7,
			       const G4Point3D& Pt8): G4BREPSolid(name)
{
  nb_of_surfaces=6;
  active=1; PlaneSolid=1;

  // Save the constructor parameters
  constructorParams[0] = Pt1;
  constructorParams[1] = Pt2;
  constructorParams[2] = Pt3;
  constructorParams[3] = Pt4;
  constructorParams[4] = Pt5;
  constructorParams[5] = Pt6;
  constructorParams[6] = Pt7;
  constructorParams[7] = Pt8;  

  InitializeBox();
}

G4BREPSolidBox::G4BREPSolidBox( __void__& a )
  : G4BREPSolid(a)
{
}

G4BREPSolidBox::~G4BREPSolidBox()
{
}

G4BREPSolidBox::G4BREPSolidBox(const G4BREPSolidBox& rhs)
  : G4BREPSolid(rhs), Rotation(rhs.Rotation)
{
  for (size_t i=0; i<8; ++i) { constructorParams[i]= rhs.constructorParams[i]; }
  InitializeBox();
}

G4BREPSolidBox& G4BREPSolidBox::operator = (const G4BREPSolidBox& rhs) 
{
  // Check assignment to self
  //
  if (this == &rhs)  { return *this; }

  // Copy base class data
  //
  G4BREPSolid::operator=(rhs);

  // Copy data
  //
  Rotation= rhs.Rotation;
  for (size_t i=0; i<8; ++i) { constructorParams[i]= rhs.constructorParams[i]; }
  InitializeBox();

  return *this;
}  

void G4BREPSolidBox::InitializeBox()
{
  SurfaceVec = new G4Surface*[6];
  G4Point3DVector PVec(4);
  G4int sense=0;

  PVec[0] = constructorParams[0];
  PVec[1] = constructorParams[1];
  PVec[2] = constructorParams[2];
  PVec[3] = constructorParams[3];  
  SurfaceVec[0] = new G4FPlane(&PVec);

  PVec[2] = constructorParams[5];
  PVec[3] = constructorParams[4];  
  SurfaceVec[1] = new G4FPlane(&PVec,0,sense);

  PVec[0] = constructorParams[1];
  PVec[1] = constructorParams[5];
  PVec[2] = constructorParams[6];
  PVec[3] = constructorParams[2];  
  SurfaceVec[2] = new G4FPlane(&PVec);

  PVec[0] = constructorParams[2];
  PVec[1] = constructorParams[6];
  PVec[2] = constructorParams[7];
  PVec[3] = constructorParams[3];  
  SurfaceVec[3] = new G4FPlane(&PVec);

  PVec[0] = constructorParams[0];
  PVec[1] = constructorParams[4];
  PVec[2] = constructorParams[7];
  PVec[3] = constructorParams[3];  
  SurfaceVec[4] = new G4FPlane(&PVec,0,sense);

  PVec[0] = constructorParams[4];
  PVec[1] = constructorParams[5];
  PVec[2] = constructorParams[6];
  PVec[3] = constructorParams[7];  
  SurfaceVec[5] = new G4FPlane(&PVec,0,sense);

  Initialize();
}

EInside G4BREPSolidBox::Inside(register const G4ThreeVector& Pt) const
{
  G4Point3D Point(Pt);

  // Get the bounding box extent
  G4Point3D min = bbox->GetBoxMin();
  min += G4Point3D(-0.5*kCarTolerance,-0.5*kCarTolerance,-0.5*kCarTolerance);

  G4Point3D max = bbox->GetBoxMax();
  max += G4Point3D(0.5*kCarTolerance,0.5*kCarTolerance,0.5*kCarTolerance);

  if( (Point.x() < min.x() || Point.x() > max.x()) ||
      (Point.y() < min.y() || Point.y() > max.y()) ||
      (Point.z() < min.z() || Point.z() > max.z())    )
    return kOutside;

  if( (Point.x() > min.x() && Point.x() < max.x())&&
      (Point.y() > min.y() && Point.y() < max.y())&&
      (Point.z() > min.z() && Point.z() < max.z())   )
    return kInside;
  
  return kSurface;
}

G4VSolid* G4BREPSolidBox::Clone() const
{
  return new G4BREPSolidBox(*this);
}

std::ostream& G4BREPSolidBox::StreamInfo(std::ostream& os) const
{
     // Streams solid contents to output stream.

     G4BREPSolid::StreamInfo( os )
     << "\n"
     << "   Pt1: " << constructorParams[0]
     << "   Pt2: " << constructorParams[1]
     << "   Pt3: " << constructorParams[2]
     << "   Pt4: " << constructorParams[3]
     << "\n   Pt5: " << constructorParams[4]
     << "   Pt6: " << constructorParams[5]
     << "   Pt7: " << constructorParams[6]
     << "   Pt8: " << constructorParams[7]
     << "\n-----------------------------------------------------------\n";

  return os;
}

