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

#include "G4VtkClipperClosedPipeline.hh"

#include "vtkActor.h"
#include "vtkAppendPolyData.h"
#include "vtkClipPolyData.h"
#include "vtkFeatureEdges.h"
#include "vtkPlane.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkStripper.h"
#include "vtkTriangleFilter.h"

G4VtkClipperClosedPipeline::G4VtkClipperClosedPipeline(
  vtkSmartPointer<vtkPolyDataAlgorithm> filterIn)
{
  finalFilter = filterIn;

  // create implicit function for clipping
  plane = vtkSmartPointer<vtkPlane>::New();
  plane->SetOrigin(0, 0, 0);
  plane->SetNormal(1, 0, 0);

  // clipper
  clipper = vtkSmartPointer<vtkClipPolyData>::New();
  clipper->SetClipFunction(plane);
  clipper->SetInputConnection(finalFilter->GetOutputPort());

  // calculate normals with a feature angle of 45 degrees
  auto filterNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
  filterNormals->SetFeatureAngle(45);
  filterNormals->SetInputConnection(clipper->GetOutputPort());

  // boundary edges
  auto boundaryEdges = vtkSmartPointer<vtkFeatureEdges>::New();
  boundaryEdges->SetInputConnection(clipper->GetOutputPort());
  boundaryEdges->BoundaryEdgesOn();
  boundaryEdges->FeatureEdgesOff();
  boundaryEdges->NonManifoldEdgesOff();
  boundaryEdges->ManifoldEdgesOff();

  // boundary strips
  auto boundaryStrips = vtkSmartPointer<vtkStripper>::New();
  boundaryStrips->SetInputConnection(boundaryEdges->GetOutputPort());

  // boundary polygon
  auto boundaryPoly = vtkSmartPointer<vtkPolyData>::New();
  boundaryPoly->SetPoints(boundaryStrips->GetOutput()->GetPoints());
  boundaryPoly->SetPolys(boundaryStrips->GetOutput()->GetLines());

  // boundary triangle filter
  auto boundaryTriangle = vtkSmartPointer<vtkTriangleFilter>::New();
  boundaryTriangle->AddInputData(boundaryPoly);

  // calculate normals with a feature angle of 45 degrees
  auto boundaryNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
  boundaryNormals->SetFeatureAngle(45);
  boundaryNormals->SetInputConnection(boundaryTriangle->GetOutputPort());

  // append filter
  appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
  appendFilter->AddInputConnection(boundaryNormals->GetOutputPort());
  appendFilter->AddInputConnection(filterNormals->GetOutputPort());

  // mapper
  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(appendFilter->GetOutputPort());
  mapper->SetColorModeToDirectScalars();
  mapper->ScalarVisibilityOn();

  // add to actor
  actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->SetVisibility(1);
}

void G4VtkClipperClosedPipeline::Modified()
{
  appendFilter->Update();
};
