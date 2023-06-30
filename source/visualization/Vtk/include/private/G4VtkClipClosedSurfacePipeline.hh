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

#ifndef G4VTKCLIPPERPIPELINE_HH
#define G4VTKCLIPPERPIPELINE_HH

#include "G4VVtkPipeline.hh"

#include "vtkSmartPointer.h"

class G4VtkVisContext;

class vtkPlane;
class vtkPolyDataMapper;
class vtkClipClosedSurface;
class vtkPolyDataAlgorithm;
class vtkActor;

class G4VtkClipClosedSurfacePipeline : public G4VVtkPipeline
{
  public:
    G4VtkClipClosedSurfacePipeline(G4String name, const G4VtkVisContext& vc,
                                   vtkSmartPointer<vtkPolyDataAlgorithm> filter,
                                   G4bool useVcColour = false);
    ~G4VtkClipClosedSurfacePipeline() override = default;

    void SetPlane(const G4Plane3D& plane);
    void SetPlane(G4double x, G4double y, G4double z, G4double nx, G4double ny, G4double nz);
    void TransformPlane(G4double dx, G4double dy, G4double dz, G4double r00, G4double r01,
                        G4double r02, G4double r10, G4double r11, G4double r12, G4double r20,
                        G4double r21, G4double r22);

    void Enable() override;
    void Disable() override;

    void Print() override { G4VVtkPipeline::Print(); };
    void Modified() override { G4VVtkPipeline::Modified(); };
    void Clear() override
    {
      if (renderer != nullptr) {
        renderer->RemoveActor(actor);
      }
      G4VVtkPipeline::Clear();
    };

    vtkSmartPointer<vtkActor> GetActor() { return actor; }

  private:
    vtkSmartPointer<vtkPlane> plane;
    vtkSmartPointer<vtkClipClosedSurface> clipper;
    vtkSmartPointer<vtkPolyDataMapper> mapper;
    vtkSmartPointer<vtkActor> actor;
};

#endif  // G4VTKCLIPPERPIPELINE_HH
