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

#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkNew.h"
#include "vtkTetra.h"
#include "vtkVoxel.h"
#include "vtkStaticCleanUnstructuredGrid.h"

#include "G4VtkUnstructuredGridPipeline.hh"
#include "G4VtkVisContext.hh"
#include "G4VtkViewer.hh"

#include "G4Mesh.hh"
#include "G4LogicalVolume.hh"
#include "G4VNestedParameterisation.hh"
#include "G4ModelingParameters.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4PseudoScene.hh"
#include "G4Transform3D.hh"
#include "G4VSolid.hh"
#include "G4Tet.hh"
#include "G4Material.hh"

G4VtkUnstructuredGridPipeline::G4VtkUnstructuredGridPipeline(G4String nameIn, const G4VtkVisContext& vc) :
  G4VVtkPipeline(nameIn, "G4VtkUnstructuredPipeline", vc, false, vc.fViewer->renderer) {

  points      = vtkSmartPointer<vtkPoints>::New();
  pointValues = vtkSmartPointer<vtkDoubleArray>::New();
  cellValues  = vtkSmartPointer<vtkDoubleArray>::New();

  cellValues->SetNumberOfComponents(4);
  cellValues->SetName("Colors");

  unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  unstructuredGrid->SetPoints(points);
  unstructuredGrid->GetCellData()->SetScalars(cellValues);

  mapper = vtkSmartPointer<vtkDataSetMapper>::New();
  mapper->SetScalarModeToUseCellData();
  mapper->SetColorModeToDirectScalars();
  mapper->SetInputData(unstructuredGrid);

  // create actor
  actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->SetVisibility(1);

  // add to renderer
  vc.fViewer->renderer->AddActor(actor);
}

void G4VtkUnstructuredGridPipeline::PseudoSceneForTetCells::AddSolid(const G4VSolid& solid) {
  if (fpPVModel->GetCurrentDepth() == fDepth) {  // Leaf-level cells only
    // Need to know it's a tet !!!! or implement G4VSceneHandler::AddSolid (const G4Tet&) !!!!
    try {
      const G4Tet& tet = dynamic_cast<const G4Tet&>(solid);
      const auto* pVisAtts = fpPVModel->GetCurrentLV()->GetVisAttributes();
      const auto& colour = pVisAtts->GetColour();

      G4ThreeVector p0, p1, p2, p3;
      tet.GetVertices(p0, p1,p2, p3);

      fpPoints->InsertNextPoint(p0.x(),p0.y(),p0.z());
      fpPoints->InsertNextPoint(p1.x(),p1.y(),p1.z());
      fpPoints->InsertNextPoint(p2.x(),p2.y(),p2.z());
      fpPoints->InsertNextPoint(p3.x(),p3.y(),p3.z());

      vtkNew<vtkTetra> tetra;
      tetra->GetPointIds()->SetId(0, 4*iCell);
      tetra->GetPointIds()->SetId(1, 4*iCell+1);
      tetra->GetPointIds()->SetId(2, 4*iCell+2);
      tetra->GetPointIds()->SetId(3, 4*iCell+3);

      fpGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());

      double cols[] = {colour.GetRed(),
                       colour.GetGreen(),
                       colour.GetBlue(),
                       colour.GetAlpha()};
      fpCellValues->InsertNextTuple(cols);

      iCell++;

    }
    catch (const std::bad_cast&) {
      G4ExceptionDescription ed;
      ed << "Called for a mesh that is not a tetrahedron mesh: " << solid.GetName();
      G4Exception("PseudoSceneForTetVertices","visman0108",JustWarning,ed);
    }
  }
}

void G4VtkUnstructuredGridPipeline::PseudoSceneForCubicalCells::AddSolid(const G4Box& box) {
  if (fpPVModel->GetCurrentDepth() == fDepth) {  // Leaf-level cells only
    try {
      const auto* pVisAtts = fpPVModel->GetCurrentLV()->GetVisAttributes();
      const auto& colour = pVisAtts->GetColour();
      const G4ThreeVector& position = fpCurrentObjectTransformation->getTranslation();

      fpPoints->InsertNextPoint(position.x()-box.GetXHalfLength(),position.y()-box.GetYHalfLength(),position.z()-box.GetZHalfLength());
      fpPoints->InsertNextPoint(position.x()-box.GetXHalfLength(),position.y()+box.GetYHalfLength(),position.z()-box.GetZHalfLength());
      fpPoints->InsertNextPoint(position.x()+box.GetXHalfLength(),position.y()+box.GetYHalfLength(),position.z()-box.GetZHalfLength());
      fpPoints->InsertNextPoint(position.x()+box.GetXHalfLength(),position.y()-box.GetYHalfLength(),position.z()-box.GetZHalfLength());

      fpPoints->InsertNextPoint(position.x()-box.GetXHalfLength(),position.y()-box.GetYHalfLength(),position.z()+box.GetZHalfLength());
      fpPoints->InsertNextPoint(position.x()-box.GetXHalfLength(),position.y()+box.GetYHalfLength(),position.z()+box.GetZHalfLength());
      fpPoints->InsertNextPoint(position.x()+box.GetXHalfLength(),position.y()+box.GetYHalfLength(),position.z()+box.GetZHalfLength());
      fpPoints->InsertNextPoint(position.x()+box.GetXHalfLength(),position.y()-box.GetYHalfLength(),position.z()+box.GetZHalfLength());

      vtkNew<vtkVoxel> voxel;
      voxel->GetPointIds()->SetId(0, 8*iCell);
      voxel->GetPointIds()->SetId(1, 8*iCell+1);
      voxel->GetPointIds()->SetId(2, 8*iCell+2);
      voxel->GetPointIds()->SetId(3, 8*iCell+3);
      voxel->GetPointIds()->SetId(4, 8*iCell+4);
      voxel->GetPointIds()->SetId(5, 8*iCell+5);
      voxel->GetPointIds()->SetId(6, 8*iCell+6);
      voxel->GetPointIds()->SetId(7, 8*iCell+7);

      fpGrid->InsertNextCell(voxel->GetCellType(), voxel->GetPointIds());

      double cols[] = {colour.GetRed(),
                       colour.GetGreen(),
                       colour.GetBlue(),
                       colour.GetAlpha()};
      fpCellValues->InsertNextTuple(cols);

      iCell++;

    }
    catch (const std::bad_cast&) {
      G4ExceptionDescription ed;
      ed << "Called for a mesh that is not a cubical mesh: " << box.GetName();
      G4Exception("PseudoSceneForCubicalVertices", "visman0108", JustWarning, ed);
    }
  }
}

void G4VtkUnstructuredGridPipeline::SetUnstructuredGridData(const G4Mesh &mesh) {

  // Modelling parameters
  G4ModelingParameters tmpMP;
  tmpMP.SetCulling(true);             // This avoids drawing transparent...
  tmpMP.SetCullingInvisible(true);    // ... or invisble volumes.

  // Physical volume model
  const auto& container = mesh.GetContainerVolume();
  const G4bool useFullExtent = true;  // To avoid calculating the extent
  G4PhysicalVolumeModel tmpPVModel(container,
                                   G4PhysicalVolumeModel::UNLIMITED,
                                   G4Transform3D(),  // so that positions are in local coordinates
                                   &tmpMP,
                                   useFullExtent);

  // Graphics scene
  if(mesh.GetMeshType() == G4Mesh::tetrahedron) {
    PseudoSceneForTetCells pseudoScene(&tmpPVModel, mesh.GetMeshDepth(), points, cellValues,
                                       unstructuredGrid);
    tmpPVModel.DescribeYourselfTo(pseudoScene);

  }
  else if(mesh.GetMeshType() == G4Mesh::nested3DRectangular) {
    PseudoSceneForCubicalCells pseudoScene(&tmpPVModel, mesh.GetMeshDepth(), points, cellValues,
                                           unstructuredGrid);
    tmpPVModel.DescribeYourselfTo(pseudoScene);
  }
  else if(mesh.GetMeshType() == G4Mesh::rectangle) {
    PseudoSceneForCubicalCells pseudoScene(&tmpPVModel, mesh.GetMeshDepth(), points, cellValues,
                                           unstructuredGrid);
    tmpPVModel.DescribeYourselfTo(pseudoScene);
  }

}