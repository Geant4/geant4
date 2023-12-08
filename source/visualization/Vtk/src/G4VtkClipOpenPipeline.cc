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

#include "G4VtkClipOpenPipeline.hh"

#include "G4VtkViewer.hh"
#include "G4VtkVisContext.hh"

#include "vtkActor.h"
#include "vtkClipPolyData.h"
#include "vtkPlane.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkProperty.h"
#include "vtkSmartPointer.h"
#include <vtkClipClosedSurface.h>

G4VtkClipOpenPipeline::G4VtkClipOpenPipeline(G4String nameIn, const G4VtkVisContext& vcIn,
                                             vtkSmartPointer<vtkPolyDataAlgorithm> filter,
                                             G4bool useVcColour)
  : G4VVtkPipeline(nameIn, G4String("G4VtkClipOpenPipeline"), vcIn, true, vcIn.fViewer->renderer)
{
  // create implicit function for clipping
  plane = vtkSmartPointer<vtkPlane>::New();
  plane->SetOrigin(0, 0, 0);
  plane->SetNormal(-1, 0, 0);

  // clipper
  clipper = vtkSmartPointer<vtkClipPolyData>::New();
  clipper->SetClipFunction(plane);
  clipper->SetInputConnection(filter->GetOutputPort());
  // clipper->SetScalarModeToColors();
  // clipper->PassPointDataOn();

  // calculate normals with a feature angle of 45 degrees
  auto filterNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
  filterNormals->SetFeatureAngle(45);
  filterNormals->SetInputConnection(clipper->GetOutputPort());

  // mapper
  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(filterNormals->GetOutputPort());
  mapper->SetColorModeToDirectScalars();
  mapper->ScalarVisibilityOn();

  // add to actor
  actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->SetVisibility(1);

  // colour parameters
  if (useVcColour) {
    actor->GetProperty()->SetOpacity(vc.alpha);
    actor->GetProperty()->SetColor(vc.red, vc.green, vc.blue);
  }

  // set actor properties from vis context
  if (vc.fDrawingStyle == G4ViewParameters::hsr) {
    actor->GetProperty()->SetRepresentationToSurface();
  }
  else if (vc.fDrawingStyle == G4ViewParameters::hlr) {
    actor->GetProperty()->SetRepresentationToSurface();
  }
  else if (vc.fDrawingStyle == G4ViewParameters::wireframe) {
    actor->GetProperty()->SetRepresentationToWireframe();
  }

  // shading parameters
  actor->GetProperty()->SetAmbient(0.2);
  actor->GetProperty()->SetDiffuse(0.7);
  actor->GetProperty()->SetSpecular(0.1);
  actor->GetProperty()->SetSpecularPower(1);

  // add to renderer
  vc.fViewer->renderer->AddActor(actor);
}

void G4VtkClipOpenPipeline::SetPlane(const G4Plane3D& planeIn)
{
  auto normal = planeIn.normal();
  auto point = planeIn.point();
  this->SetPlane(point.x(), point.y(), point.z(), normal.x(), normal.y(), normal.z());
}

void G4VtkClipOpenPipeline::SetPlane(G4double x, G4double y, G4double z, G4double nx, G4double ny,
                                     G4double nz)
{
  plane->SetOrigin(x, y, z);
  plane->SetNormal(nx, ny, nz);
  Modified();
}

void G4VtkClipOpenPipeline::TransformPlane(G4double dx, G4double dy, G4double dz, G4double r00,
                                           G4double r01, G4double r02, G4double r10, G4double r11,
                                           G4double r12, G4double r20, G4double r21, G4double r22)
{
  auto o = plane->GetOrigin();
  auto n = plane->GetNormal();

  SetPlane(r00 * o[0] + r01 * o[1] + r02 * o[2] + dx, r10 * o[0] + r11 * o[1] + r12 * o[2] + dy,
           r20 * o[0] + r21 * o[1] + r22 * o[2] + dz, r00 * n[0] + r01 * n[1] + r02 * n[2],
           r10 * n[0] + r11 * n[1] + r12 * n[2], r20 * n[0] + r21 * n[1] + r22 * n[2]);
}

void G4VtkClipOpenPipeline::Enable()
{
  actor->SetVisibility(1);
}
void G4VtkClipOpenPipeline::Disable()
{
  actor->SetVisibility(0);
}