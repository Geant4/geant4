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
// Created by Stewart Boogert on 13/11/2021.
//

#ifndef G4VTKSTORE_HH
#define G4VTKSTORE_HH

#include "G4Transform3D.hh"
#include "G4ViewParameters.hh"

#include "vtkActor.h"
#include "vtkActor2D.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkClipClosedSurface.h"
#include "vtkClipPolyData.h"
#include "vtkDataSetMapper.h"
#include "vtkDoubleArray.h"
#include "vtkFeatureEdges.h"
#include "vtkNamedColors.h"
#include "vtkPlane.h"
#include "vtkPlaneCollection.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataMapper2D.h"
#include "vtkPolyDataNormals.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkStripper.h"
#include "vtkStructuredGrid.h"
#include "vtkTensorGlyphColor.h"
#include "vtkTriangleFilter.h"
#include "vtkVertexGlyphFilter.h"
#include <vtkSmartPointer.h>

#include <map>
#include <vector>

class G4VViewer;
class G4VtkViewer;
class G4Polyline;
class G4Text;
class G4Circle;
class G4Square;
class G4Polyhedron;
class G4VisAttributes;
class G4VtkVisContext;
class G4Mesh;

class G4VtkPolydataPipeline;
class G4VtkPolydataPolylinePipeline;
class G4VtkPolydataPolyline2DPipeline;
class G4VtkPolydataSpherePipeline;
class G4VtkPolydataPolygonPipeline;
class G4VtkPolydataInstanceTensorPipeline;
class G4VtkPolydataInstanceAppendPipeline;
class G4VtkPolydataInstanceBakePipeline;
class G4VtkClipCloseSurfacePipeline;
class G4VtkImagePipeline;
class G4VtkTextPipeline;
class G4VtkText2DPipeline;
class G4VtkUnstructuredGridPipeline;

template<typename>
class vtkSmartPointer;

class G4VtkTensorGlyphPolydataPipeline;

class G4VtkStore
{
  public:
    G4VtkStore(G4String name);
    ~G4VtkStore();

    void Modified();
    void Clear();
    void ClearNonG4();
    void Print();

    void AddPrimitive(const G4Polyline& polyline, const G4VtkVisContext& vc);
    void AddPrimitive(const G4Text& text, const G4VtkVisContext& vc);
    void AddPrimitive(const G4Circle& circle, const G4VtkVisContext& vc);
    void AddPrimitive(const G4Square& square, const G4VtkVisContext& vc);
    void AddPrimitiveSeparate(const G4Polyhedron& polyhedron, const G4VtkVisContext& vc);
    void AddPrimitiveTensorGlyph(const G4Polyhedron& polyhedron, const G4VtkVisContext& vc);
    void AddPrimitiveAppend(const G4Polyhedron& polyhedron, const G4VtkVisContext& vc);
    void AddPrimitiveTransformBake(const G4Polyhedron& polyhedron, const G4VtkVisContext& vc);
    void AddCompound(const G4Mesh& mesh,  const G4VtkVisContext& vc);

    void UpdatePlanePipelines(G4String name, G4String type, const G4Plane3D);

    void AddClipper(G4String name, const G4Plane3D& plane);
    void UpdateClipper(G4String name, const G4Plane3D& plane);
    void RemoveClipper(G4String name);

    void AddCutter(G4String name, const G4Plane3D& plane);
    void UpdateCutter(G4String name, const G4Plane3D& plane);
    void RemoveCutter(G4String name);

    void AddNonG4ObjectImage(const G4String& fileName, const G4VtkVisContext& vc);
    void AddNonG4ObjectPolydata(const G4String fileName, const G4VtkVisContext& vc);

    void GetBounds(G4double maxBound[6]);

    void AddToRenderer(vtkRenderer* renderer);

    std::map<std::size_t, std::shared_ptr<G4VtkPolydataPolylinePipeline>>& GetPolylinePipeMap()
    {
      return polylinePipeMap;
    }
    std::map<std::size_t, std::shared_ptr<G4VtkPolydataPolyline2DPipeline>>& GetPolyline2DPipeMap()
    {
      return polyline2DPipeMap;
    }
    std::map<std::size_t, std::shared_ptr<G4VtkPolydataSpherePipeline>>& GetPolydataSpherePipeMap()
    {
      return circlePipeMap;
    }
    std::map<std::size_t, std::shared_ptr<G4VtkPolydataPolygonPipeline>>&
    GetPolydataPolygonPipeMap()
    {
      return squarePipeMap;
    }
    std::map<std::size_t, std::shared_ptr<G4VtkPolydataPipeline>>& GetSeparatePipeMap()
    {
      return separatePipeMap;
    }
    std::map<std::size_t, std::shared_ptr<G4VtkPolydataInstanceTensorPipeline>>& GetTensorPipeMap()
    {
      return tensorGlyphPipeMap;
    }
    std::map<std::size_t, std::shared_ptr<G4VtkPolydataInstanceAppendPipeline>>& GetAppendPipeMap()
    {
      return appendPipeMap;
    }
    std::map<std::size_t, std::shared_ptr<G4VtkPolydataInstanceBakePipeline>>& GetBakePipeMap()
    {
      return bakePipeMap;
    }

  private:
    G4String name;
    std::map<std::size_t, std::shared_ptr<G4VtkPolydataPolylinePipeline>> polylinePipeMap;
    std::map<std::size_t, std::shared_ptr<G4VtkPolydataPolyline2DPipeline>> polyline2DPipeMap;
    std::map<std::size_t, std::shared_ptr<G4VtkPolydataSpherePipeline>> circlePipeMap;
    std::map<std::size_t, std::shared_ptr<G4VtkPolydataPolygonPipeline>> squarePipeMap;
    std::map<std::size_t, std::shared_ptr<G4VtkTextPipeline>> textPipeMap;
    std::map<std::size_t, std::shared_ptr<G4VtkText2DPipeline>> text2DPipeMap;

    std::map<std::size_t, std::shared_ptr<G4VtkPolydataPipeline>> separatePipeMap;
    std::map<std::size_t, std::shared_ptr<G4VtkPolydataInstanceTensorPipeline>> tensorGlyphPipeMap;
    std::map<std::size_t, std::shared_ptr<G4VtkPolydataInstanceAppendPipeline>> appendPipeMap;
    std::map<std::size_t, std::shared_ptr<G4VtkPolydataInstanceBakePipeline>> bakePipeMap;

    std::map<G4String, std::shared_ptr<G4VtkUnstructuredGridPipeline>> ugridPipeMap;

    std::map<G4String, std::shared_ptr<G4VtkImagePipeline>> imagePipeMap;
    std::map<std::size_t, std::shared_ptr<G4VtkPolydataPipeline>> sideloadPipeMap;
};

#endif  // G4VTKSTORE_HH
