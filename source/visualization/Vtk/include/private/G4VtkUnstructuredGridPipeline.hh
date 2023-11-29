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

#ifndef G4VTKUNSTRUCTUREDGRIDPIPELINE_HH
#define G4VTKUNSTRUCTUREDGRIDPIPELINE_HH

#include "G4VVtkPipeline.hh"
#include "G4PseudoScene.hh"

#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"
#include "vtkStaticCleanUnstructuredGrid.h"

class vtkDoubleArray;
class vtkDataSetMapper;
class vtkActor;

class G4Mesh;

class G4VtkUnstructuredGridPipeline : public G4VVtkPipeline
{
  public:
    G4VtkUnstructuredGridPipeline(G4String name, const G4VtkVisContext& vc);

    ~G4VtkUnstructuredGridPipeline() = default;

    void Modified() {};
    void Clear() {};
    void Print() {};

    void Enable() {};
    void Disable() {};

    void SetUnstructuredGridData(const G4Mesh &mesh);

  protected:
    class PseudoSceneVtkBase: public G4PseudoScene {
      public:
        PseudoSceneVtkBase(G4PhysicalVolumeModel* pvModel,  // input
                           G4int depth,
                           vtkPoints *points,
                           vtkDoubleArray *cellValues,
                           vtkUnstructuredGrid *unstructuredGrid)                     // input...the following are outputs by reference
          : fpPVModel(pvModel), fDepth(depth), fpPoints(points),
            fpGrid(unstructuredGrid), fpCellValues(cellValues)
        {}
      protected:
        using G4PseudoScene::AddSolid;  // except for...
        void AddSolid(const G4VSolid& /*solid*/) override {};
        void AddSolid(const G4Box& /*box*/) override{};
        void ProcessVolume(const G4VSolid&) override {
          // Do nothing if uninteresting solids found, e.g., the container if not marked invisible.
        }
        G4PhysicalVolumeModel* fpPVModel;
        G4int fDepth;
        vtkSmartPointer<vtkPoints> fpPoints;
        vtkSmartPointer<vtkUnstructuredGrid> fpGrid;
        vtkSmartPointer<vtkDoubleArray> fpCellValues;
        G4int iCell = 0;
    };
    class PseudoSceneForTetCells: public PseudoSceneVtkBase {
      public:
        PseudoSceneForTetCells(G4PhysicalVolumeModel* pvModel,  // input
                               G4int depth,
                               vtkPoints *points,
                               vtkDoubleArray *cellValues,
                               vtkUnstructuredGrid *unstructuredGrid)                     // input...the following are outputs by reference
          : PseudoSceneVtkBase(pvModel,depth,points,cellValues,unstructuredGrid)
        {}
      private:
        using G4PseudoScene::AddSolid;  // except for...
        void AddSolid(const G4VSolid& solid) override;
        void ProcessVolume(const G4VSolid&) override {}
    };
    class PseudoSceneForCubicalCells: public PseudoSceneVtkBase {
      public:
        PseudoSceneForCubicalCells(G4PhysicalVolumeModel* pvModel,  // input
                                    G4int depth,
                                    vtkPoints *points,
                                    vtkDoubleArray *cellValues,
                                    vtkUnstructuredGrid *unstructuredGrid)                     // input...the following are outputs by reference
          : PseudoSceneVtkBase(pvModel,depth,points,cellValues,unstructuredGrid)
        {}
      private:
        using G4PseudoScene::AddSolid;  // except for...
        void AddSolid(const G4Box& box) override;
        void ProcessVolume(const G4VSolid&) override {}
    };

  private:
    vtkSmartPointer<vtkPoints> points;
    vtkSmartPointer<vtkDoubleArray> pointValues;
    vtkSmartPointer<vtkDoubleArray> cellValues;
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
    vtkSmartPointer<vtkStaticCleanUnstructuredGrid> cleanUnstructuredGrid;
    vtkSmartPointer<vtkDataSetMapper> mapper;
    vtkSmartPointer<vtkActor> actor;
};

#endif  // G4VTKUNSTRUCTUREDGRIDPIPELINE_HH
