// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AxesModel.cc,v 1.2 2001-05-03 11:16:42 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  3rd April 2001
// Model which knows how to draw axes.

#include "G4AxesModel.hh"

#include "G4ModelingParameters.hh"
#include "G4VGraphicsScene.hh"
#include "G4VisAttributes.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"

G4AxesModel::~G4AxesModel () {}

G4AxesModel::G4AxesModel
(G4double x0, G4double y0, G4double z0, G4double length):
  fX0(x0), fY0(y0), fZ0(z0), fLength(length) {
  fGlobalTag = "G4AxesModel";
  fGlobalDescription = fGlobalTag;
}

void G4AxesModel::DescribeYourselfTo (G4VGraphicsScene& scene) {

  G4Polyline x_axis, y_axis, z_axis;

  G4Colour cx(1.0, 0.0, 0.0); // color for x-axis
  G4Colour cy(0.0, 1.0, 0.0); // color for y-axis
  G4Colour cz(0.0, 0.0, 1.0); // color for z-axis

  G4VisAttributes ax(cx);     // VA for x-axis
  G4VisAttributes ay(cy);     // VA for y-axis
  G4VisAttributes az(cz);     // VA for z-axis

  //----- Draw x-axis
  x_axis.SetVisAttributes(&ax);
  x_axis.push_back(G4Point3D(fX0,fY0,fZ0));
  x_axis.push_back(G4Point3D((fX0 + fLength),fY0,fZ0));
  scene.AddPrimitive(x_axis);

  //----- Draw y-axis
  y_axis.SetVisAttributes(&ay);
  y_axis.push_back(G4Point3D(fX0,fY0,fZ0));
  y_axis.push_back(G4Point3D(fX0,(fY0 + fLength),fZ0));
  scene.AddPrimitive(y_axis);

  //----- Draw z-axis
  z_axis.SetVisAttributes(&az);
  z_axis.push_back(G4Point3D(fX0,fY0,fZ0));
  z_axis.push_back(G4Point3D(fX0,fY0,(fZ0 + fLength)));
  scene.AddPrimitive(z_axis);
}
