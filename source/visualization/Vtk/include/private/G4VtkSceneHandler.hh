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
// John Allison  5th April 2001
// A template for a simplest possible graphics driver.
//?? Lines or sections marked like this require specialisation for your driver.

#ifndef G4VTKSCENEHANDLER_HH
#define G4VTKSCENEHANDLER_HH

#include "G4VSceneHandler.hh"
#include "G4VtkStore.hh"
#include "G4VtkViewer.hh"

#include <map>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra-semi"
#include "vtkActor.h"
#include "vtkBillboardTextActor3D.h"
#include "vtkCleanPolyData.h"
#include "vtkDoubleArray.h"
#include "vtkFeatureEdges.h"
#include "vtkFollower.h"
#include "vtkGlyph2D.h"
#include "vtkGlyph3D.h"
#include "vtkImageData.h"
#include "vtkLine.h"
#include "vtkMatrix4x4.h"
#include "vtkNamedColors.h"
#include "vtkOpenGLGPUVolumeRayCastMapper.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkProperty.h"
#include "vtkRegularPolygonSource.h"
#include "vtkScalarsToColors.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkTensorGlyph.h"
#include "vtkTensorGlyphColor.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkTriangleFilter.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkVolume.h"
#pragma GCC diagnostic pop

class G4VtkSceneHandler : public G4VSceneHandler
{
  public:
    G4VtkSceneHandler(G4VGraphicsSystem& system, const G4String& name);
    ~G4VtkSceneHandler() override = default;

    ////////////////////////////////////////////////////////////////
    // Required implementation of pure virtual functions...

    using G4VSceneHandler::AddPrimitive;
    void AddPrimitive(const G4Polyline&) override;
    void AddPrimitive(const G4Text&) override;
    void AddPrimitive(const G4Circle&) override;
    void AddPrimitive(const G4Square&) override;
    void AddPrimitive(const G4Polyhedron&) override;

    // Further optional AddPrimitive methods.  Explicitly invoke base
    // class methods if not otherwise defined to avoid warnings about
    // hiding of base class methods.
    void AddPrimitive(const G4Polymarker& polymarker) override
    {
      G4VSceneHandler::AddPrimitive(polymarker);
    }

    using G4VSceneHandler::AddSolid;
    void AddSolid(const G4Box& box) override;

    using G4VSceneHandler::AddCompound;
    void AddCompound(const G4Mesh& mesh) override;

    void Modified();
    void ClearStore() override;
    void ClearTransientStore() override;
    G4VtkVisContext MakeDefaultVisContext();
    G4VtkStore& GetStore() { return store; }
    G4VtkStore& GetTransientStore() { return transientStore; }

    virtual void Print();

    void SetPolyhedronPipeline(const G4String& str);

  protected:
    static G4int fSceneIdCount;  // Counter for Vtk scene handlers.

    G4VtkStore store = G4VtkStore("perm");
    G4VtkStore transientStore = G4VtkStore("trans");

    G4String polyhedronPipelineType;

  private:
    friend class G4VtkViewer;
};

#endif
