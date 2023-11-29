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
//
// 
// John Allison  15th July 2012
// Model that knows how to draw an arrow.

#include "G4ArrowModel.hh"

#include "G4PhysicalConstants.hh"
#include "G4VGraphicsScene.hh"
#include "G4VisAttributes.hh"
#include "G4Tubs.hh"
#include "G4Tet.hh"
#include "G4Polyhedron.hh"
#include "G4Vector3D.hh"
#include "G4Point3D.hh"
#include "G4Transform3D.hh"
#include "G4GeometryTolerance.hh"

#include <cmath>

G4ArrowModel::~G4ArrowModel ()
{
  delete fpHeadPolyhedron;
  delete fpShaftPolyhedron;
}

G4ArrowModel::G4ArrowModel
(G4double x1, G4double y1, G4double z1,
 G4double x2, G4double y2, G4double z2,
 G4double width, const G4Colour& colour,
 const G4String& description,
 G4int lineSegmentsPerCircle,
 const G4Transform3D& transform)
: fpShaftPolyhedron(nullptr)
, fpHeadPolyhedron(nullptr)
, fTransform(transform)
{
  fType = "G4ArrowModel";
  fGlobalTag = fType;
  fGlobalDescription = fType + ": " + description;
  fExtent = G4VisExtent
    (std::min(x1,x2),
     std::max(x1,x2),
     std::min(y1,y2),
     std::max(y1,y2),
     std::min(z1,z2),
     std::max(z1,z2));

  // Force number of line segments per circle (aka number of rotation steps)
  G4int tempN = G4Polyhedron::GetNumberOfRotationSteps();
  G4Polyhedron::SetNumberOfRotationSteps(lineSegmentsPerCircle);

  // Make a cylinder slightly shorter than the arrow length so that it
  // doesn't stick out of the head.
  const G4double tolerance = G4GeometryTolerance::GetInstance()->GetRadialTolerance();
  
  G4double totalLength = std::hypot(x2-x1, y2-y1, z2-z1);
  if (totalLength < tolerance)
    {totalLength = tolerance;}
  
  G4double shaftRadius = width/6.;
  if (shaftRadius < tolerance)
    {shaftRadius = tolerance;}

  // case 1 - arrow length >> width -> arrow head is width and 1.5x width tall
  // case 2 - arrow length <  width -> arrow head is made to be 0.5x length
  G4double arrowLength = std::min(1.5*width, 0.5*totalLength);

  G4double shaftLength = totalLength - arrowLength;
  if (shaftLength < 2*tolerance)
    {shaftLength = 2*tolerance;}

  const G4Tubs shaft("shaft",0.,shaftRadius,0.5*shaftLength,0.,twopi);
  fpShaftPolyhedron = shaft.CreatePolyhedron();
  // translate the polyhedron down w.r.t. the centre of the whole arrow
  if (fpShaftPolyhedron)
    {fpShaftPolyhedron->Transform(G4Translate3D(0,0,-0.5*arrowLength));}

  // Locate the head at +halfShaftLength.
  const G4double zHi  = 0.5*totalLength;
  const G4double zLow = zHi - arrowLength;
  const G4double rExt = 0.5*width;
  const G4double xExt = std::sqrt(3.)*rExt/2.;
  const G4Tet head("head",
                   G4ThreeVector(0.,0.,zHi),
                   G4ThreeVector(0.,rExt,zLow),
                   G4ThreeVector(xExt,-rExt/2.,zLow),
                   G4ThreeVector(-xExt,-rExt/2.,zLow));
  fpHeadPolyhedron = head.CreatePolyhedron();

  // Transform to position
  const G4Vector3D arrowDirection = G4Vector3D(x2-x1,y2-y1,z2-z1).unit();
  const G4double theta = arrowDirection.theta();
  const G4double phi = arrowDirection.phi();
  const G4Point3D arrowCentre(0.5*(x1+x2),0.5*(y1+y2),0.5*(z1+z2));
  const G4Transform3D tr =
    G4Translate3D(arrowCentre) * G4RotateZ3D(phi) * G4RotateY3D(theta);
  if (fpShaftPolyhedron) fpShaftPolyhedron->Transform(tr);
  if (fpHeadPolyhedron) fpHeadPolyhedron->Transform(tr);

  G4VisAttributes va;
  va.SetColour(colour);
  va.SetForceSolid(true);
  if (fpShaftPolyhedron) fpShaftPolyhedron->SetVisAttributes(va);
  if (fpHeadPolyhedron) fpHeadPolyhedron->SetVisAttributes(va);

  // Restore number of line segments per circle
  G4Polyhedron::SetNumberOfRotationSteps(tempN);
}

void G4ArrowModel::DescribeYourselfTo (G4VGraphicsScene& sceneHandler)
{
  if (fpShaftPolyhedron && fpHeadPolyhedron) {
    sceneHandler.BeginPrimitives(fTransform);
    sceneHandler.AddPrimitive(*fpShaftPolyhedron);
    sceneHandler.AddPrimitive(*fpHeadPolyhedron);
    sceneHandler.EndPrimitives();
  }
}
