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

#include "G4VtkPolydataPipeline.hh"

#include "G4Normal3D.hh"
#include "G4Point3D.hh"
#include "G4Polyhedron.hh"
#include "G4Polyline.hh"
#include "G4ViewParameters.hh"
#include "G4VtkViewer.hh"
#include "G4VtkVisContext.hh"

#include "vtkActor.h"
#include "vtkCleanPolyData.h"
#include "vtkLine.h"
#include "vtkMatrix4x4.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkProperty.h"
#include "vtkTriangleFilter.h"

std::size_t G4VtkPolydataPipeline::MakeHash(const G4Polyhedron& polyhedron,
                                            const G4VtkVisContext& vc)
{
  std::size_t hash = std::hash<G4Polyhedron>{}(polyhedron);

  std::hash_combine(hash, std::hash<G4double>{}(vc.fTransform.dx()));
  std::hash_combine(hash, std::hash<G4double>{}(vc.fTransform.dy()));
  std::hash_combine(hash, std::hash<G4double>{}(vc.fTransform.dz()));

  std::hash_combine(hash, std::hash<G4double>{}(vc.fTransform.xx()));
  std::hash_combine(hash, std::hash<G4double>{}(vc.fTransform.xy()));
  std::hash_combine(hash, std::hash<G4double>{}(vc.fTransform.xz()));

  std::hash_combine(hash, std::hash<G4double>{}(vc.fTransform.yx()));
  std::hash_combine(hash, std::hash<G4double>{}(vc.fTransform.yy()));
  std::hash_combine(hash, std::hash<G4double>{}(vc.fTransform.yz()));

  std::hash_combine(hash, std::hash<G4double>{}(vc.fTransform.zx()));
  std::hash_combine(hash, std::hash<G4double>{}(vc.fTransform.zy()));
  std::hash_combine(hash, std::hash<G4double>{}(vc.fTransform.zz()));

  return hash;
}

G4VtkPolydataPipeline::G4VtkPolydataPipeline(G4String nameIn, const G4VtkVisContext& vc)
  : G4VVtkPipeline(nameIn, "G4VtkPolydataPipeline", vc, false, vc.fViewer->renderer)
{
  // Set pipeline type
  SetTypeName(G4String("G4VtkPolydataPipeline"));

  polydataPoints = vtkSmartPointer<vtkPoints>::New();
  polydataCells = vtkSmartPointer<vtkCellArray>::New();
  polydata = vtkSmartPointer<vtkPolyData>::New();

  polydata->SetPoints(polydataPoints);
  polydata->SetPolys(polydataCells);

  // clean input polydata
  auto filterClean = vtkSmartPointer<vtkCleanPolyData>::New();
  filterClean->PointMergingOn();
  filterClean->AddInputData(polydata);
  AddFilter(filterClean);

  // ensure triangular mesh
  auto filterTriangle = vtkSmartPointer<vtkTriangleFilter>::New();
  filterTriangle->SetInputConnection(filterClean->GetOutputPort());
  AddFilter(filterTriangle);

  // calculate normals with a feature angle of 45 degrees
  auto filterNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
  filterNormals->SetFeatureAngle(45);
  filterNormals->SetInputConnection(filterTriangle->GetOutputPort());
  AddFilter(filterNormals);

  // mapper
  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(GetFinalFilter()->GetOutputPort());
  mapper->SetColorModeToDirectScalars();

  // add to actor
  actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->SetVisibility(1);

  // set actor properties from vis context
  if (vc.fDrawingStyle == G4ViewParameters::hsr) {
  }
  else if (vc.fDrawingStyle == G4ViewParameters::hlr) {
  }
  else if (vc.fDrawingStyle == G4ViewParameters::wireframe) {
    actor->GetProperty()->SetRepresentationToWireframe();
  }

  // add to renderer
  vc.fViewer->renderer->AddActor(GetActor());
}

void G4VtkPolydataPipeline::Enable()
{
  actor->SetVisibility(1);
}

void G4VtkPolydataPipeline::Disable()
{
  actor->SetVisibility(0);
}

void G4VtkPolydataPipeline::Print()
{
  G4cout << "G4VtkPolydataPipeline filters (";
  for (const auto& f : filters)
    G4cout << f->GetInformation() << ",";
  G4cout << ")" << G4endl;

  G4VVtkPipeline::Print();
}

void G4VtkPolydataPipeline::Modified()
{
  actor->Modified();
  polydata->Modified();
  mapper->Update();

  G4VVtkPipeline::Modified();
}

void G4VtkPolydataPipeline::Clear()
{
  renderer->RemoveActor(actor);
  G4VVtkPipeline::Clear();
}

void G4VtkPolydataPipeline::SetPolydata(const G4Polyhedron& polyhedron)
{
  G4bool notLastFace;
  int iVert = 0;
  do {
    G4Point3D vertex[4];
    G4int edgeFlag[4];
    G4Normal3D normals[4];
    G4int nEdges;
    notLastFace = polyhedron.GetNextFacet(nEdges, vertex, edgeFlag, normals);

    vtkSmartPointer<vtkIdList> poly = vtkSmartPointer<vtkIdList>::New();
    // loop over vertices
    for (int i = 0; i < nEdges; i++) {
      polydataPoints->InsertNextPoint(vertex[i].x(), vertex[i].y(), vertex[i].z());
      poly->InsertNextId(iVert);
      iVert++;
    }
    polydataCells->InsertNextCell(poly);

  } while (notLastFace);
}

void G4VtkPolydataPipeline::SetPolydata(const G4Polyline& polyline)
{
  // Data data
  const size_t nLines = polyline.size();

  for (size_t i = 0; i < nLines; ++i) {
    auto id = polydataPoints->InsertNextPoint(polyline[i].x(), polyline[i].y(), polyline[i].z());

    if (i < nLines - 1) {
      vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
      line->GetPointIds()->SetId(0, id);
      line->GetPointIds()->SetId(1, id + 1);
      polydataCells->InsertNextCell(line);
    }
  }
}

void G4VtkPolydataPipeline::SetPolydataData(const G4Point3D& p)
{
  SetPolydataData(p.x(), p.y(), p.z());
}

void G4VtkPolydataPipeline::SetPolydataData(double x, double y, double z)
{
  polydataPoints->InsertNextPoint(x, y, z);
}

void G4VtkPolydataPipeline::SetActorTransform(G4double dx, G4double dy, G4double dz, G4double r00,
                                              G4double r01, G4double r02, G4double r10,
                                              G4double r11, G4double r12, G4double r20,
                                              G4double r21, G4double r22)
{
  // create transform
  auto transform = vtkSmartPointer<vtkMatrix4x4>::New();
  double transformArray[16] = {r00, r01, r02, dx, r10, r11, r12, dy,
                               r20, r21, r22, dz, 0.,  0.,  0.,  1.};
  transform->DeepCopy(transformArray);
  actor->SetUserMatrix(transform);
}

void G4VtkPolydataPipeline::SetActorColour(G4double r, G4double g, G4double b, G4double a)
{
  actor->GetProperty()->SetColor(r, g, b);
  actor->GetProperty()->SetOpacity(a);
}
