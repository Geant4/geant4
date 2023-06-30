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

#include "G4VtkPolydataInstanceBakePipeline.hh"

#include "G4VtkViewer.hh"
#include "G4VtkVisContext.hh"

#include "vtkActor.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"

std::size_t G4VtkPolydataInstanceBakePipeline::MakeHash(const G4Polyhedron& polyhedron,
                                                        const G4VtkVisContext& vc)
{
  // Get view parameters that the user can force through the vis attributes, thereby over-riding the
  // current view parameter.
  const G4VisAttributes* pVA =
    vc.fViewer->GetApplicableVisAttributes(polyhedron.GetVisAttributes());
  G4Colour colour = pVA->GetColour();

  // Hash the vis attributes
  std::size_t hash = std::hash<G4double>{}(colour.GetAlpha());
  std::size_t shash = std::hash<G4double>{}(vc.fDrawingStyle);

  std::hash_combine(hash, shash);

  return hash;
}

G4VtkPolydataInstanceBakePipeline::G4VtkPolydataInstanceBakePipeline(G4String nameIn,
                                                                     const G4VtkVisContext& vc)
  : G4VtkPolydataInstancePipeline(nameIn, vc)
{
  iVert = 0;
  iFace = 0;
  polydataPointData = vtkSmartPointer<vtkDoubleArray>::New();
  polydataPointData->SetNumberOfComponents(3);
  polydataPointData->SetName("Colors");
  polydata->GetPointData()->SetScalars(polydataPointData);

  mapper->SetColorModeToDirectScalars();

  actor->GetProperty()->SetOpacity(vc.alpha);
}

void G4VtkPolydataInstanceBakePipeline::SetPolydata(const G4Polyhedron& polyhedron)
{
  polyhedronPrototype = &polyhedron;
}

void G4VtkPolydataInstanceBakePipeline::Print()
{
  G4cout << "G4VtkPolydataInstanceBakePipeline " << GetName() << G4endl;
  G4VtkPolydataPipeline::Print();
}

void G4VtkPolydataInstanceBakePipeline::addInstance(G4double dx, G4double dy, G4double dz,
                                                    G4double r00, G4double r01, G4double r02,
                                                    G4double r10, G4double r11, G4double r12,
                                                    G4double r20, G4double r21, G4double r22,
                                                    G4double r, G4double g, G4double b, G4double a,
                                                    const G4String& name)
{
  G4VtkPolydataInstancePipeline::addInstance(dx, dy, dz, r00, r01, r02, r10, r11, r12, r20, r21,
                                             r22, r, g, b, a, name);

  vtkIdType vStart;
  vtkIdType vEnd;
  vtkIdType fStart;
  vtkIdType fEnd;

  vStart = iVert;
  fStart = iFace;

  G4bool notLastFace;
  do {
    G4Point3D vertex[4];
    G4int edgeFlag[4];
    G4Normal3D normals[4];
    G4int nEdges;
    notLastFace = polyhedronPrototype->GetNextFacet(nEdges, vertex, edgeFlag, normals);

    vtkSmartPointer<vtkIdList> poly = vtkSmartPointer<vtkIdList>::New();
    // loop over vertices
    for (int i = 0; i < nEdges; i++) {
      // note : G4Transform3D does not have a matrix element constructor
      G4Point3D bakedVertex =
        G4Point3D(vertex[i].x() * r00 + vertex[i].y() * r01 + vertex[i].z() * r02 + dx,
                  vertex[i].x() * r10 + vertex[i].y() * r11 + vertex[i].z() * r12 + dy,
                  vertex[i].x() * r20 + vertex[i].y() * r21 + vertex[i].z() * r22 + dz);
      polydataPoints->InsertNextPoint(bakedVertex.x(), bakedVertex.y(), bakedVertex.z());
      poly->InsertNextId(iVert);
      iVert++;
      double cols[] = {r, g, b, 0.5};
      polydataPointData->InsertNextTuple(cols);
    }

    polydataCells->InsertNextCell(poly);
    iFace++;

  } while (notLastFace);

  vEnd = iVert;
  fEnd = iFace;

  // Add vertex and face to maps to allow for removal later
  instanceVertexMap[name] = std::pair<vtkIdType, vtkIdType>(vStart, vEnd);
  instanceFaceMap[name] = std::pair<vtkIdType, vtkIdType>(fStart, fEnd);
}

void G4VtkPolydataInstanceBakePipeline::removeInstance(const G4String& /*name*/) {}
