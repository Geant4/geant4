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
  G4double shaftLength = std::sqrt
    (std::pow(x2-x1,2)+std::pow(y2-y1,2)+std::pow(z2-z1,2));
  if (shaftLength < tolerance) shaftLength = tolerance;
  G4double shaftRadius = width/2.;
  if (shaftRadius > shaftLength/100.) shaftRadius = shaftLength/100.;
  if (shaftRadius < tolerance) shaftRadius = tolerance;
  const G4double halfShaftLength = shaftLength/2.;
  const G4double halfReduction = 4.*shaftRadius;
  G4double halfLength = halfShaftLength - halfReduction;
  if (halfLength < tolerance) halfLength = tolerance;
  const G4Tubs shaft("shaft",0.,shaftRadius,halfLength,0.,twopi);
  fpShaftPolyhedron = shaft.CreatePolyhedron();
  // Move it a little so that the tail is at z = -halfShaftLength.
  if (fpShaftPolyhedron)
    fpShaftPolyhedron->Transform(G4Translate3D(0,0,-halfReduction));

  // Locate the head at +halfShaftLength.
  const G4double zHi  = halfShaftLength;
  const G4double zLow = halfShaftLength - 12.*shaftRadius;
  const G4double rExt = 8. * shaftRadius;
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
