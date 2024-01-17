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
#include "vtkCharArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetMapper.h"
#include "vtkUnstructuredGridVolumeMapper.h"
#include "vtkUnstructuredGridVolumeRayCastMapper.h"
#include "vtkUnstructuredGridVolumeZSweepMapper.h"
#include "vtkOpenGLProjectedTetrahedraMapper.h"
#include "vtkActor.h"
#include "vtkColorTransferFunction.h"
#include "vtkPiecewiseFunction.h"
#include "vtkLookupTable.h"
#include "vtkVolumeProperty.h"
#include "vtkVolume.h"
#include "vtkNew.h"
#include "vtkTetra.h"
#include "vtkVoxel.h"
// #include "vtkCleanUnstructuredGrid.h"
#include "vtkDataSetTriangleFilter.h"
#include "vtkClipDataSet.h"

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

G4VtkUnstructuredGridPipeline::G4VtkUnstructuredGridPipeline(G4String nameIn, const G4VtkVisContext& vcIn) :
  G4VVtkPipeline(nameIn, "G4VtkUnstructuredPipeline", vcIn, false, vcIn.fViewer->renderer) {

  points      = vtkSmartPointer<vtkPoints>::New();

  pointColourValues = vtkSmartPointer<vtkDoubleArray>::New();
  cellColourValues  = vtkSmartPointer<vtkDoubleArray>::New();
  pointColourValues->SetNumberOfComponents(4);
  cellColourValues->SetNumberOfComponents(4);

  pointColourIndices = vtkSmartPointer<vtkDoubleArray>::New();
  cellColourIndices  = vtkSmartPointer<vtkDoubleArray>::New();
  pointColourIndices->SetNumberOfComponents(1);
  cellColourIndices->SetNumberOfComponents(1);

  colourLUT = vtkSmartPointer<vtkDiscretizableColorTransferFunction>::New();
  colourLUT->DiscretizeOn();

  unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  unstructuredGrid->SetPoints(points);
  unstructuredGrid->GetPointData()->SetScalars(pointColourValues);
  unstructuredGrid->GetCellData()->SetScalars(cellColourValues);

  // clean filter
#if 0
  clean = vtkSmartPointer<vtkStaticCleanUnstructuredGrid>::New();
  clean->SetInputData(unstructuredGrid);
  clean->ToleranceIsAbsoluteOff();
  clean->SetTolerance(1e-6);
#endif

  // Clip filter
  clip = vtkSmartPointer<vtkClipDataSet>::New();
  vtkNew<vtkPlane> plane;
  clip->SetClipFunction(plane);
  clip->SetInputData(unstructuredGrid);

  // Triangle filter
  auto tri = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
  tri->SetInputData(unstructuredGrid);

  // create dataset mapper
  mapper = vtkSmartPointer<vtkDataSetMapper>::New();
  mapper->SetScalarModeToUseCellData();
  mapper->SetColorModeToDirectScalars();
  mapper->SetInputData(unstructuredGrid);
  //mapper->SetInputConnection(clip->GetOutputPort());
  //mapper->SetInputConnection(tri->GetOutputPort());

  // create volume mapper
  volumeMapper = vtkSmartPointer<vtkUnstructuredGridVolumeRayCastMapper>::New();
  volumeMapper->SetScalarModeToUseCellData();
  volumeMapper->SetInputConnection(tri->GetOutputPort());

  // create volume properties
  auto volumeProp = vtkSmartPointer<vtkVolumeProperty>::New();
  volumeProp->SetColor(colourLUT);

  // create actor
  actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->SetVisibility(1);

  // create volume
  volume = vtkSmartPointer<vtkVolume>::New();
  volume->SetMapper(volumeMapper);
  //volume->SetProperty(volumeProp);
  volume->SetVisibility(1);

  // add to renderer
  vc.fViewer->renderer->AddActor(actor);
  //vc.fViewer->renderer->AddVolume(volume);
}

void G4VtkUnstructuredGridPipeline::PseudoSceneForTetCells::AddSolid(const G4VSolid& solid) {
  if (fpPVModel->GetCurrentDepth() == fDepth) {  // Leaf-level cells only
    // Need to know it's a tet !!!! or implement G4VSceneHandler::AddSolid (const G4Tet&) !!!!
    try {
      const G4Tet& tet = dynamic_cast<const G4Tet&>(solid);
      const auto* pVisAtts = fpPVModel->GetCurrentLV()->GetVisAttributes();
      const auto& colour = pVisAtts->GetColour();

      G4ThreeVector p0, p1, p2, p3;
      tet.GetVertices(p0, p1, p2, p3);

      G4int idx[4] = {-1, -1, -1, -1};

      if (iCell > 0 && iCell < 100000000) {
#if 0
        G4ThreeVector pts[4] = {p0, p1, p2, p3};
        for (G4int i = 0; i < 4; i++) {
          auto h = MakeHash(pts[i]);
          if (fpPointMap.find(h) == fpPointMap.end()) {
            fpPointMap.insert(std::make_pair(h, iPoint));
            fpPoints->InsertNextPoint(pts[i].x(), pts[i].y(), pts[i].z());
            fpPointVector.push_back(pts[i]);
            idx[i] = iPoint;
            iPoint++;
          }
          else {
            idx[i] = fpPointMap[h];
          }
        }
#endif

#if 1
        fpPoints->InsertNextPoint(p0.x(), p0.y(), p0.z());
        fpPoints->InsertNextPoint(p1.x(), p1.y(), p1.z());
        fpPoints->InsertNextPoint(p2.x(), p2.y(), p2.z());
        fpPoints->InsertNextPoint(p3.x(), p3.y(), p3.z());
        iPoint += 4;

        idx[0] = 4 * iCellAdd;
        idx[1] = 4 * iCellAdd + 1;
        idx[2] = 4 * iCellAdd + 2;
        idx[3] = 4 * iCellAdd + 3;
#endif

        vtkNew<vtkTetra> tetra;
        tetra->GetPointIds()->SetId(0, idx[0]);
        tetra->GetPointIds()->SetId(1, idx[1]);
        tetra->GetPointIds()->SetId(2, idx[2]);
        tetra->GetPointIds()->SetId(3, idx[3]);

        fpGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());

        auto hash = MakeHash(colour);
        unsigned short int iColour = 0;
        if (fpColourMap.find(hash) == fpColourMap.end()) {
          iColour = fpColourMap.size();
          fpColourMap.insert(std::make_pair(hash,iColour));
          double c[4] = {colour.GetRed(), colour.GetGreen(), colour.GetBlue(), colour.GetAlpha()};
          fpColourLUT->SetIndexedColorRGBA(iColour,c);
        }
        else {
          iColour = fpColourMap[hash];
        }

        fpCellColourValues->InsertNextTuple4(colour.GetRed(), colour.GetGreen(), colour.GetBlue(), colour.GetAlpha());
        fpCellColourIndices->InsertNextTuple1(iColour);

        iCellAdd++;
      }
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
      fpPoints->InsertNextPoint(position.x()+box.GetXHalfLength(),position.y()-box.GetYHalfLength(),position.z()-box.GetZHalfLength());
      fpPoints->InsertNextPoint(position.x()+box.GetXHalfLength(),position.y()+box.GetYHalfLength(),position.z()-box.GetZHalfLength());

      fpPoints->InsertNextPoint(position.x()-box.GetXHalfLength(),position.y()-box.GetYHalfLength(),position.z()+box.GetZHalfLength());
      fpPoints->InsertNextPoint(position.x()-box.GetXHalfLength(),position.y()+box.GetYHalfLength(),position.z()+box.GetZHalfLength());
      fpPoints->InsertNextPoint(position.x()+box.GetXHalfLength(),position.y()-box.GetYHalfLength(),position.z()+box.GetZHalfLength());
      fpPoints->InsertNextPoint(position.x()+box.GetXHalfLength(),position.y()+box.GetYHalfLength(),position.z()+box.GetZHalfLength());

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
      fpCellColourValues->InsertNextTuple(cols);
      //fpCellColourIndices->InsertNextTuple1(iColour);


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
    PseudoSceneForTetCells pseudoScene(&tmpPVModel, mesh.GetMeshDepth(), points,
                                       pointColourValues, cellColourValues,
                                       pointColourIndices, cellColourIndices,
                                       colourLUT, unstructuredGrid);
    tmpPVModel.DescribeYourselfTo(pseudoScene);
    //G4cout << pseudoScene.GetNumberOfPoints() << " " << pseudoScene.GetNumberOfCells() << " " << pseudoScene.GetNumberOfAddedCells() << G4endl;
    //pseudoScene.DumpColourMap();
  }
  else if(mesh.GetMeshType() == G4Mesh::nested3DRectangular) {
    PseudoSceneForCubicalCells pseudoScene(&tmpPVModel, mesh.GetMeshDepth(), points,
                                           pointColourValues, cellColourValues,
                                           pointColourIndices, cellColourIndices,
                                           colourLUT, unstructuredGrid);
    tmpPVModel.DescribeYourselfTo(pseudoScene);
  }
  else if(mesh.GetMeshType() == G4Mesh::rectangle) {
    PseudoSceneForCubicalCells pseudoScene(&tmpPVModel, mesh.GetMeshDepth(), points,
                                           pointColourValues, cellColourValues,
                                           pointColourIndices, cellColourIndices,
                                           colourLUT, unstructuredGrid);
    tmpPVModel.DescribeYourselfTo(pseudoScene);
  }

}