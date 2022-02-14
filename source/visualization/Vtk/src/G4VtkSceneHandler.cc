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

#include "G4VtkSceneHandler.hh"

#include "G4PhysicalVolumeModel.hh"
#include "G4LogicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Polyline.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Polyhedron.hh"
#include "G4Mesh.hh"
#include "G4PseudoScene.hh"
#include "G4UnitsTable.hh"

#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4Box.hh"

#include <stdlib.h>

namespace std
{
  inline void hash_combine(std::size_t) {}

  template <typename T, typename... Rest>
  inline void hash_combine(std::size_t &seed, const T &v, Rest... rest) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    std::hash_combine(seed, rest...);
  }

  template<> struct hash<G4VisAttributes> {
    std::size_t operator()(const G4VisAttributes &va) const {
      using std::size_t;
      using std::hash;

      std::size_t h = 0;

      std::hash_combine(h,va.IsVisible());
      std::hash_combine(h,va.IsDaughtersInvisible());
      std::hash_combine(h,va.GetColour().GetRed());
      std::hash_combine(h,va.GetColour().GetGreen());
      std::hash_combine(h,va.GetColour().GetBlue());
      std::hash_combine(h,va.GetColour().GetAlpha());
      std::hash_combine(h,static_cast<int>(va.GetLineStyle()));

      return h;
    }
  };

  template<> struct hash<G4Polyhedron> {
    std::size_t operator()(const G4Polyhedron &ph) const {
      using std::size_t;
      using std::hash;

      G4bool notLastFace;
      G4Point3D vertex[4];
      G4int edgeFlag[4];
      G4Normal3D normals[4];
      G4int nEdges;

      std::size_t h = 0;

      do {
        notLastFace = ph.GetNextFacet(nEdges, vertex, edgeFlag, normals);

        for (int i = 0; i < nEdges; i++) {
          std::size_t hx = std::hash<double>()(vertex[i].x());
          std::size_t hy = std::hash<double>()(vertex[i].y());
          std::size_t hz = std::hash<double>()(vertex[i].z());
          std::hash_combine(h,hx);
          std::hash_combine(h,hy);
          std::hash_combine(h,hz);
        }
      } while (notLastFace);

      return h;
    }
  };
}

G4int G4VtkSceneHandler::fSceneIdCount = 0;
// Counter for XXX scene handlers.

G4VtkSceneHandler::G4VtkSceneHandler(G4VGraphicsSystem& system,
                                     const G4String& name) :
                                     G4VSceneHandler(system, fSceneIdCount++, name)
{}

#ifdef G4VTKDEBUG
void G4VtkSceneHandler::PrintThings() {
  G4cout << "  with transformation " << fObjectTransformation.xx() << G4endl;
  if (fpModel) {
    G4cout << " from " << fpModel->GetCurrentDescription()
           << " (tag " << fpModel->GetCurrentTag()
           << ')';
  }
  else {
    G4cout << "(not from a model)";
  }
  G4PhysicalVolumeModel* pPVModel = dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
  if (pPVModel) {
    G4cout << "\n  current physical volume:        " << pPVModel->GetCurrentPV()->GetName()
           << "\n  current logical volume :        " << pPVModel->GetCurrentLV()->GetName() // There might be a problem with the LV pointer if this is a G4LogicalVolumeModel
           << "\n  current depth of geometry tree: " << pPVModel->GetCurrentDepth();
  }
  G4cout << G4endl;
}
#endif

void G4VtkSceneHandler::AddPrimitive(const G4Polyline& polyline) {

  G4VSceneHandler::MarkerSizeType sizeType;
  // GetMarkerSize(text,sizeType);
  if(fProcessing2D) {sizeType = screen;}
  else              {sizeType = world;}

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA = fpViewer -> GetApplicableVisAttributes(polyline.GetVisAttributes());
  G4Color  colour            = pVA->GetColour();
  G4double opacity           = colour.GetAlpha();
  G4double lineWidth         = pVA->GetLineWidth();

#ifdef G4VTKDEBUG
  G4cout << "=================================" << G4endl;
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Polyline& polyline) called>             vis hash: " << sizeType << " " << std::hash<G4VisAttributes>{}(*pVA) << G4endl;
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Polyline& polyline) called>             sizeType: " << sizeType << G4endl;
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Polyline& polyline) called>            isVisible: " << pVA->IsVisible() << G4endl;
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Polyline& polyline) called> isDaughtersInvisible: " << pVA->IsDaughtersInvisible() << G4endl;
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Polyline& polyline) called>               colour: " << colour.GetRed() << " " << colour.GetGreen() << " " << colour.GetBlue()  << G4endl;
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Polyline& polyline) called>                alpha: " << colour.GetAlpha() << G4endl;
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Polyline& polyline) called>            lineWidth: " << lineWidth << G4endl;
#endif

  if(sizeType == world) {
    std::size_t hash = std::hash<G4VisAttributes>{}(*pVA);
    if (polylineVisAttributesMap.find(hash) == polylineVisAttributesMap.end()) {
      polylineVisAttributesMap.insert(std::pair<std::size_t, const G4VisAttributes *>(hash, pVA));

      vtkSmartPointer <vtkPoints> data = vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer <vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
      vtkSmartPointer <vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer <vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      vtkSmartPointer <vtkActor> actor = vtkSmartPointer<vtkActor>::New();

      polyData->SetPoints(data);
      polyData->SetLines(lines);
      mapper->SetInputData(polyData);
      actor->SetMapper(mapper);

      // Setup actor and mapper
      actor->GetProperty()->SetLineWidth(lineWidth);
      actor->GetProperty()->SetColor(colour.GetRed(), colour.GetGreen(),
                                     colour.GetBlue());
      actor->GetProperty()->SetOpacity(opacity);
      actor->SetVisibility(1);
      // actor->GetProperty()->BackfaceCullingOn();
      // actor->GetProperty()->FrontfaceCullingOn();

      auto pVtkViewer = dynamic_cast<G4VtkViewer*>(fpViewer);
      pVtkViewer->renderer->AddActor(actor);

      polylineDataMap.insert(
        std::pair<std::size_t, vtkSmartPointer<vtkPoints>>(hash, data));
      polylineLineMap.insert(
        std::pair<std::size_t, vtkSmartPointer<vtkCellArray>>(hash, lines));
      polylinePolyDataMap.insert(
        std::pair<std::size_t, vtkSmartPointer<vtkPolyData>>(hash, polyData));
      polylinePolyDataMapperMap.insert(
        std::pair<std::size_t, vtkSmartPointer<vtkPolyDataMapper>>(hash,
                                                                   mapper));
      polylinePolyDataActorMap.insert(
        std::pair<std::size_t, vtkSmartPointer<vtkActor>>(hash, actor));
    }

    // Data data
    const size_t nLines = polyline.size();

    for (size_t i = 0; i < nLines; ++i) {
      auto id = polylineDataMap[hash]->InsertNextPoint(polyline[i].x(),
                                                       polyline[i].y(),
                                                       polyline[i].z());
      if (i < nLines - 1) {
        vtkSmartPointer <vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, id);
        line->GetPointIds()->SetId(1, id + 1);
        polylineLineMap[hash]->InsertNextCell(line);
      }
    }
  }
  else if (sizeType == screen ) {

  }
}

void G4VtkSceneHandler::AddPrimitive(const G4Text& text) {

  G4VSceneHandler::MarkerSizeType sizeType;
  if(fProcessing2D) {sizeType = screen;}
  else              {sizeType = world;}

  const G4VisAttributes* pVA = fpViewer -> GetApplicableVisAttributes(text.GetVisAttributes());
  G4Color colour             = pVA->GetColour();
  G4double opacity           = colour.GetAlpha();
  // G4Text::Layout layout      = text.GetLayout();
  // G4double xOffset           = text.GetXOffset();
  // G4double yOffset           = text.GetYOffset();

  double x = text.GetPosition().x();
  double y = text.GetPosition().y();
  double z = text.GetPosition().z();

#ifdef G4VTKDEBUG
  G4cout << "=================================" << G4endl;
  G4cout << "G4VtkSeneHandler::AddPrimitive(const G4Text& text) called>     text: " << text.GetText() << " sizeType:" << sizeType << " " << fProcessing2D << G4endl;
  G4cout << "G4VtkSeneHandler::AddPrimitive(const G4Text& text) called>   colour: " << colour.GetRed() << " " << colour.GetBlue() << " " << colour.GetGreen()  << G4endl;
  G4cout << "G4VtkSeneHandler::AddPrimitive(const G4Text& text) called>    alpha: " << colour.GetAlpha() << G4endl;
  G4cout << "G4VtkSeneHandler::AddPrimitive(const G4Text& text) called> position: " << x << " " << y << " " << z << G4endl;
#endif

  switch (sizeType) {
    default:
    case (screen): {
      vtkSmartPointer <vtkTextActor> actor = vtkSmartPointer<vtkTextActor>::New();
      actor->SetInput(text.GetText().c_str());
      actor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
      // actor->SetTextScaleModeToViewport();
      actor->SetPosition((x+1.)/2.0, (y+1.)/2.);
      actor->GetTextProperty()->SetFontSize(text.GetScreenSize());
      actor->GetTextProperty()->SetColor(colour.GetRed(), colour.GetBlue(), colour.GetGreen());
      actor->GetTextProperty()->SetOpacity(opacity);

      auto *pVtkViewer = dynamic_cast<G4VtkViewer *>(fpViewer);
      pVtkViewer->renderer->AddActor(actor);
      break;
    }
    case world: {
      vtkSmartPointer <vtkBillboardTextActor3D> actor = vtkSmartPointer<vtkBillboardTextActor3D>::New();
      actor->SetInput(text.GetText().c_str());
      actor->SetPosition(x, y, z);
      actor->GetTextProperty()->SetFontSize(text.GetScreenSize());
      actor->GetTextProperty()->SetColor(colour.GetRed(), colour.GetBlue(), colour.GetGreen());
      actor->GetTextProperty()->SetOpacity(opacity);

      auto *pVtkViewer = dynamic_cast<G4VtkViewer*>(fpViewer);
      pVtkViewer->renderer->AddActor(actor);
      break;
    }
  }
}

void G4VtkSceneHandler::AddPrimitive(const G4Circle& circle) {

  MarkerSizeType sizeType;
  G4double size= GetMarkerSize(circle, sizeType);
  if(fProcessing2D) {sizeType = screen;}
  else              {sizeType = world;}

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes *pVA = fpViewer->GetApplicableVisAttributes(circle.GetVisAttributes());
  G4Color colour = pVA->GetColour();
  G4double opacity = colour.GetAlpha();
  // G4bool isVisible = pVA->IsVisible();

#ifdef G4VTKDEBUG
  G4cout << "=================================" << G4endl;
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Circle& circle) called>           " << " radius:" << size << " sizeType:" << sizeType << G4endl;
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Circle& circle) called>   colour: " << colour.GetRed() << " " << colour.GetBlue() << " " << colour.GetGreen()  << G4endl;
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Circle& circle) called>    alpha: " << colour.GetAlpha() << G4endl;
#endif

  if (sizeType == world) {
    std::size_t hash = std::hash<G4VisAttributes>{}(*pVA);
    if (circleVisAttributesMap.find(hash) == circleVisAttributesMap.end()) {
      circleVisAttributesMap.insert(std::pair<std::size_t, const G4VisAttributes *>(hash, pVA));

      vtkSmartPointer<vtkPoints> data = vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkVertexGlyphFilter> filter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
      vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();

      polyData->SetPoints(data);
      filter->SetInputData(polyData);
      mapper->SetInputConnection(filter->GetOutputPort());
      actor->SetMapper(mapper);

      // Setup actor and mapper
      actor->GetProperty()->SetColor(colour.GetRed(), colour.GetGreen(), colour.GetBlue());
      actor->GetProperty()->SetOpacity(opacity);
      actor->SetVisibility(1);
      actor->GetProperty()->SetRenderPointsAsSpheres(true);
      actor->GetProperty()->SetPointSize(size*5);

      auto *pVtkViewer = dynamic_cast<G4VtkViewer *>(fpViewer);
      pVtkViewer->renderer->AddActor(actor);

      circleDataMap.insert(std::pair<std::size_t,vtkSmartPointer<vtkPoints>>(hash, data));
      circlePolyDataMap.insert(std::pair<std::size_t, vtkSmartPointer < vtkPolyData>>(hash, polyData));
      circleFilterMap.insert(std::pair<std::size_t, vtkSmartPointer<vtkVertexGlyphFilter>>(hash, filter));
      circlePolyDataMapperMap.insert(std::pair<std::size_t, vtkSmartPointer < vtkPolyDataMapper>>(hash, mapper));
      circlePolyDataActorMap.insert(std::pair<std::size_t, vtkSmartPointer < vtkActor>>(hash, actor));
    }

    // Data data point
    const CLHEP::HepRotation rot = fObjectTransformation.getRotation();
    G4Point3D posPrime = rot*circle.GetPosition();
    circleDataMap[hash]->InsertNextPoint(fObjectTransformation.dx() + posPrime.x(),
                                         fObjectTransformation.dy() + posPrime.y(),
                                         fObjectTransformation.dz() + posPrime.z());
  }
  else if (sizeType == screen) {

  }
}

void G4VtkSceneHandler::AddPrimitive(const G4Square& square) {

  MarkerSizeType sizeType;
  G4double size = GetMarkerSize (square, sizeType);

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA = fpViewer -> GetApplicableVisAttributes(square.GetVisAttributes());
  G4Color colour    = pVA->GetColour();
  G4double opacity  = colour.GetAlpha();
  // G4bool  isVisible = pVA->IsVisible();

  // Draw in world coordinates.
  vtkSmartPointer<vtkRegularPolygonSource> polygonSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
  polygonSource->SetNumberOfSides(4);
  polygonSource->SetRadius(size);


#ifdef G4VTKDEBUG
  G4cout << "=================================" << G4endl;
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Square& square) called" << G4endl;
  G4cout << square.GetPosition().x() << " " << square.GetPosition().y() << " " << square.GetPosition().z() << G4endl;
  //PrintThings();
#endif

  if (sizeType == world) {
    std::size_t hash = std::hash<G4VisAttributes>{}(*pVA);
    if (squareVisAttributesMap.find(hash) == squareVisAttributesMap.end()) {
      squareVisAttributesMap.insert(std::pair<std::size_t, const G4VisAttributes *>(hash, pVA));

      vtkSmartPointer<vtkPoints>              data = vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkPolyData>        polyData = vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkVertexGlyphFilter> filter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
      vtkSmartPointer<vtkPolyDataMapper>    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      vtkSmartPointer<vtkActor>              actor = vtkSmartPointer<vtkActor>::New();

      polyData->SetPoints(data);
      filter->SetInputData(polyData);
      mapper->SetInputConnection(filter->GetOutputPort());
      actor->SetMapper(mapper);

      // Setup actor and mapper
      actor->GetProperty()->SetColor(colour.GetRed(), colour.GetGreen(), colour.GetBlue());
      actor->GetProperty()->SetOpacity(opacity);
      actor->SetVisibility(1);
      actor->GetProperty()->SetPointSize(size*5);

      auto *pVtkViewer = dynamic_cast<G4VtkViewer *>(fpViewer);
      pVtkViewer->renderer->AddActor(actor);

      squareDataMap.insert(std::pair<std::size_t,vtkSmartPointer<vtkPoints>>(hash, data));
      squarePolyDataMap.insert(std::pair<std::size_t, vtkSmartPointer < vtkPolyData>>(hash, polyData));
      squareFilterMap.insert(std::pair<std::size_t, vtkSmartPointer<vtkVertexGlyphFilter>>(hash, filter));
      squarePolyDataMapperMap.insert(std::pair<std::size_t, vtkSmartPointer < vtkPolyDataMapper>>(hash, mapper));
      squarePolyDataActorMap.insert(std::pair<std::size_t, vtkSmartPointer < vtkActor>>(hash, actor));
    }

    // Data data point
    const CLHEP::HepRotation rot = fObjectTransformation.getRotation();
    G4Point3D posPrime = rot*square.GetPosition();
    squareDataMap[hash]->InsertNextPoint(fObjectTransformation.dx() + posPrime.x(),
                                         fObjectTransformation.dy() + posPrime.y(),
                                         fObjectTransformation.dz() + posPrime.z());
  }
  else if (sizeType == screen) {

  }
}
void G4VtkSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) {
  AddPrimitiveTensorGlyph(polyhedron);
}

void G4VtkSceneHandler::AddPrimitiveTensorGlyph(const G4Polyhedron& polyhedron) {
  // only have a single polydata object for each LV but bake in the PV
  // transformation so clippers and cutters can be implemented

  //Get colour, etc..
  if (polyhedron.GetNoFacets() == 0) {return;}

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes *pVA = fpViewer->GetApplicableVisAttributes(polyhedron.GetVisAttributes());
  G4Color                       colour = pVA->GetColour();
  // G4bool                     isVisible = pVA->IsVisible();
  // G4double                   lineWidth = pVA->GetLineWidth();
  // G4VisAttributes::LineStyle lineStyle = pVA->GetLineStyle();

  // G4double lineWidthScale = drawing_style.GetGlobalLineWidthScale();

  // Get view parameters that the user can force through the vis attributes, thereby over-riding the current view parameter.
  G4ViewParameters::DrawingStyle drawing_style = GetDrawingStyle(pVA);
  //G4bool isAuxEdgeVisible = GetAuxEdgeVisible (pVA);

#ifdef G4VTKDEBUG
  G4cout << "=================================" << G4endl;
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) called> " << G4endl;
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) called>     colour:" << colour.GetRed() << " " << colour.GetBlue() << " " << colour.GetGreen() << G4endl;
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) called>      alpha:" << colour.GetAlpha() << G4endl;
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) called>  lineWidth:" << lineWidth << G4endl;
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) called>  lineStyle:" << lineStyle << G4endl;
#endif

  //std::size_t vhash = std::hash<G4VisAttributes>{}(*pVA);

  std::size_t vhash = 0;
  std::hash_combine(vhash,static_cast<int>(drawing_style));

  std::size_t phash = std::hash<G4Polyhedron>{}(polyhedron);

  std::size_t hash = 0;
  std::hash_combine(hash, phash);
  std::hash_combine(hash, vhash);

  if (polyhedronPolyDataMap.find(hash) == polyhedronPolyDataMap.end()) {

    vtkSmartPointer<vtkPoints>     points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray>   polys = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();

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
        points->InsertNextPoint(vertex[i].x(), vertex[i].y(), vertex[i].z());
        poly->InsertNextId(iVert);
        iVert++;
      }
      polys->InsertNextCell(poly);

    } while (notLastFace);

    polydata->SetPoints(points);
    polydata->SetPolys(polys);

    polyhedronDataMap.insert(std::pair<std::size_t, vtkSmartPointer<vtkPoints>>(hash, points));
    polyhedronPolyMap.insert(std::pair<std::size_t, vtkSmartPointer<vtkCellArray>>(hash, polys));
    polyhedronPolyDataMap.insert(std::pair<std::size_t, vtkSmartPointer<vtkPolyData>>(hash, polydata));
    polyhedronPolyDataCountMap.insert(std::pair<std::size_t, std::size_t>(hash, 0));

    vtkSmartPointer<vtkPoints>          instancePosition = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkDoubleArray>     instanceRotation = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray>       instanceColors = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkPolyDataMapper>    instanceMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor>              instanceActor = vtkSmartPointer<vtkActor>::New();
    instanceColors->SetName("colors");

    instanceColors->SetNumberOfComponents(4);
    instanceRotation->SetNumberOfComponents(9);

    vtkSmartPointer<vtkPolyData> instancePolyData = vtkSmartPointer<vtkPolyData>::New();
    instancePolyData->SetPoints(instancePosition);
    instancePolyData->GetPointData()->SetTensors(instanceRotation);
    instancePolyData->GetPointData()->SetVectors(instanceColors);
    instancePolyData->GetPointData()->SetScalars(instanceColors);

    vtkSmartPointer<vtkCleanPolyData> filterClean = vtkSmartPointer<vtkCleanPolyData>::New();
    filterClean->PointMergingOn();
    filterClean->AddInputData(polydata);

    vtkSmartPointer<vtkTriangleFilter> filterTriangle = vtkSmartPointer<vtkTriangleFilter>::New();
    filterTriangle->SetInputConnection(filterClean->GetOutputPort());

    vtkSmartPointer<vtkPolyDataNormals> filterNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
    filterNormals->SetFeatureAngle(45);
    filterNormals->SetInputConnection(filterTriangle->GetOutputPort());

    vtkSmartPointer<vtkFeatureEdges> filterEdge = vtkSmartPointer<vtkFeatureEdges>::New();
    filterEdge->SetFeatureEdges(1);
    filterEdge->SetManifoldEdges(0);
    filterEdge->SetBoundaryEdges(0);
    filterEdge->SetFeatureAngle(45);  // TODO need to have a function and command to set  this
    filterEdge->SetInputConnection(filterTriangle->GetOutputPort());

    vtkSmartPointer<vtkTensorGlyphColor> tensorGlyph = vtkSmartPointer<vtkTensorGlyphColor>::New();
    tensorGlyph->SetInputData(instancePolyData);
    tensorGlyph->SetSourceConnection(filterNormals->GetOutputPort());
    tensorGlyph->ColorGlyphsOn();
    tensorGlyph->ScalingOff();
    tensorGlyph->ThreeGlyphsOff();
    tensorGlyph->ExtractEigenvaluesOff();
    tensorGlyph->SetColorModeToScalars();
    tensorGlyph->Update();

    instanceMapper->SetInputData(tensorGlyph->GetOutput());
    instanceMapper->SetColorModeToDirectScalars();

    instanceActor->SetMapper(instanceMapper);
    // instanceActor->GetProperty()->SetLineWidth(10);
    instanceActor->SetVisibility(1);

    if(drawing_style == G4ViewParameters::hsr) {

    }

    if(drawing_style == G4ViewParameters::hlr) {
    }

    if(drawing_style == G4ViewParameters::wireframe) {
      instanceActor->GetProperty()->SetRepresentationToWireframe();
    }

    auto *pVtkViewer = dynamic_cast<G4VtkViewer *>(fpViewer);
    pVtkViewer->renderer->AddActor(instanceActor);

    instancePositionMap.insert(std::pair<std::size_t, vtkSmartPointer<vtkPoints>>(hash, instancePosition));
    instanceRotationMap.insert(std::pair<std::size_t, vtkSmartPointer<vtkDoubleArray>>(hash, instanceRotation));
    instanceColoursMap.insert(std::pair<std::size_t, vtkSmartPointer<vtkDoubleArray>>(hash, instanceColors));
    instancePolyDataMap.insert(std::pair<std::size_t, vtkSmartPointer<vtkPolyData>>(hash,instancePolyData));
    instanceActorMap.insert(std::pair<std::size_t, vtkSmartPointer<vtkActor>>(hash, instanceActor));
    instanceTensorGlyphMap.insert(std::pair<std::size_t, vtkSmartPointer<vtkTensorGlyphColor>>(hash, tensorGlyph));
  }

  polyhedronPolyDataCountMap[hash]++;

  double red   = colour.GetRed();
  double green = colour.GetGreen();
  double blue  = colour.GetBlue();
  double alpha = colour.GetAlpha();

  instanceColoursMap[hash]->InsertNextTuple4(red, green, blue, alpha);

  instancePositionMap[hash]->InsertNextPoint(fObjectTransformation.dx(),
                                              fObjectTransformation.dy(),
                                              fObjectTransformation.dz());

  G4Transform3D fInvObjTrans = fObjectTransformation.inverse();

  instanceRotationMap[hash]->InsertNextTuple9(fInvObjTrans.xx(), fInvObjTrans.xy(),fInvObjTrans.xz(),
                                               fInvObjTrans.yx(), fInvObjTrans.yy(),fInvObjTrans.yz(),
                                               fInvObjTrans.zx(), fInvObjTrans.zy(),fInvObjTrans.zz());

}

void G4VtkSceneHandler::AddPrimitiveBakedTransform(const G4Polyhedron& polyhedron)
{
  // only have a single polydata object for each LV but bake in the PV
  // transformation so clippers and cutters can be implemented
  AddPrimitiveTensorGlyph(polyhedron);
}
void G4VtkSceneHandler::Modified() {
  for (auto it = polylineDataMap.begin(); it != polylineDataMap.end(); it++)
  {
    it->second->Modified();
  }

  for (auto it = polylineLineMap.begin(); it != polylineLineMap.end(); it++)
  {
    it->second->Modified();
  }
  for (auto it = circleDataMap.begin(); it != circleDataMap.end(); it++)
  {
    it->second->Modified();
  }

  for (auto it = squareDataMap.begin(); it != squareDataMap.end(); it++)
  {
    it->second->Modified();
  }

  for(auto it = instancePositionMap.begin(); it != instancePositionMap.end(); it++)
  {
    it->second->Modified();
  }

  for(auto it = instanceRotationMap.begin(); it != instanceRotationMap.end(); it++)
  {
    it->second->Modified();
  }

  for(auto it = instanceColoursMap.begin(); it != instanceColoursMap.end(); it++)
  {
    it->second->Modified();
  }

  for(auto it = instancePolyDataMap.begin(); it != instancePolyDataMap.end(); it++)
  {
    it->second->Modified();
  }

  for(auto it = instanceActorMap.begin(); it != instanceActorMap.end(); it++)
  {
    it->second->Modified();
  }

  for(auto it = instanceTensorGlyphMap.begin(); it != instanceTensorGlyphMap.end(); it++)
  {
    it->second->Update();
  }

#ifdef G4VTKDEBUG
  G4cout << "G4VtkSceneHandler::Modified()    polyline styles: " << polylineVisAttributesMap.size() << G4endl;
  for (auto it = polylineLineMap.begin(); it != polylineLineMap.end(); it++)
  {
    G4cout << "G4VtkSceneHandler::Modified()   polyline segments: "
           << it->second->GetNumberOfCells() << G4endl;
  }
  G4cout << "G4VtkSceneHandler::Modified()       circle styles: " << circleVisAttributesMap.size() << G4endl;
  for (auto it = circleDataMap.begin(); it != circleDataMap.end(); it++)
  {
    G4cout << "G4VtkSceneHandler::Modified()             circles: "
           << it->second->GetNumberOfPoints() << G4endl;
  }

  G4cout << "G4VtkSceneHandler::Modified()      square styles: " << squareVisAttributesMap.size() << G4endl;
  for (auto it = squareDataMap.begin(); it != squareDataMap.end(); it++)
  {
    G4cout << "G4VtkSceneHanler::Modified()           squares: "
           << it->second->GetNumberOfPoints() << G4endl;
  }

  G4cout << "G4VtkSceneHandler::Modified()   unique polyhedra: " << polyhedronDataMap.size() << G4endl;

  int nPlacements = 0;
  int nCells = 0;
  for (auto it = polyhedronPolyDataMap.begin(); it != polyhedronPolyDataMap.end(); it++) {
    G4cout << "G4VtkSceneHandler::Modified()  polyhedronPolyData: " << it->second->GetPoints()->GetNumberOfPoints() << " " << it->second->GetPolys()->GetNumberOfCells() << " " << polyhedronPolyDataCountMap[it->first] <<G4endl;
    nCells      += it->second->GetPolys()->GetNumberOfCells()*polyhedronPolyDataCountMap[it->first];
    nPlacements += polyhedronPolyDataCountMap[it->first];
  }
  G4cout << "G4VtkSceneHandler::Modified()  polyhedronPolyData: " << nPlacements << " " << nCells << G4endl;
#endif
}

void G4VtkSceneHandler::Clear() {

  polylineVisAttributesMap.clear();
  polylineDataMap.clear();
  polylineLineMap.clear();
  polylinePolyDataMap.clear();
  polylinePolyDataMapperMap.clear();
  polylinePolyDataActorMap.clear();

  for (auto& v : polylineDataMap)
    {v.second->Reset();}
  for (auto& v : polylineLineMap)
    {v.second->Reset();}

  circleVisAttributesMap.clear();
  circleDataMap.clear();
  circlePolyDataMap.clear();
  circleFilterMap.clear();
  circlePolyDataMapperMap.clear();
  circlePolyDataActorMap.clear();

  squareVisAttributesMap.clear();
  squareDataMap.clear();
  squarePolyDataMap.clear();
  squareFilterMap.clear();
  squarePolyDataMapperMap.clear();
  squarePolyDataActorMap.clear();

  for (auto& v : squareDataMap)
    {v.second->Reset();}

  polyhedronVisAttributesMap.clear();
  polyhedronDataMap.clear();
  polyhedronPolyMap.clear();
  polyhedronPolyDataMap.clear();
  polyhedronPolyDataCountMap.clear();

  instancePositionMap.clear();
  instanceRotationMap.clear();
  instanceColoursMap.clear();
  instancePolyDataMap.clear();
  instanceTensorGlyphMap.clear();
  instanceActorMap.clear();
}

void G4VtkSceneHandler::AddSolid (const G4Box& box) {

  G4VSceneHandler::AddSolid(box);

#if 0
  return;

  const G4VModel* pv_model = GetModel();
  if (!pv_model) { return ; }

  G4PhysicalVolumeModel* pPVModel = dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
  if (!pPVModel) { return ; }

  //-- debug information
  if(1) {
    G4VPhysicalVolume *pv = pPVModel->GetCurrentPV();
    G4LogicalVolume   *lv = pv->GetLogicalVolume();
    G4cout << "name=" << box.GetName() << " volumeType=" << pv->VolumeType() << " pvName=" << pv->GetName() << " lvName=" << lv->GetName() << " multiplicity=" << pv->GetMultiplicity() << " isparametrised=" << pv->IsParameterised() << " isreplicated=" << pv->IsReplicated() << " parametrisation=" << pv->GetParameterisation() << G4endl;
  }

  if(0) {
    G4Material *mat = pPVModel->GetCurrentMaterial();
    G4String name   = mat->GetName();
    G4double dens   = mat->GetDensity()/(g/cm3);
    G4int copyNo    = pPVModel->GetCurrentPV()->GetCopyNo();
    G4int depth     = pPVModel->GetCurrentDepth();
    G4cout << "    name    : " << box.GetName() << G4endl;
    G4cout << "    copy no.: " << copyNo << G4endl;
    G4cout << "    depth   : " << depth << G4endl;
    G4cout << "    density : " << dens << " [g/cm3]" << G4endl;
    G4cout << "    location: " << pPVModel->GetCurrentPV()->GetObjectTranslation() << G4endl;
    G4cout << "    Multiplicity        : " << pPVModel->GetCurrentPV()->GetMultiplicity() << G4endl;
    G4cout << "    Is replicated?      : " << pPVModel->GetCurrentPV()->IsReplicated() << G4endl;
    G4cout << "    Is parameterised?   : " << pPVModel->GetCurrentPV()->IsParameterised() << G4endl;
    G4cout << "    top phys. vol. name : " << pPVModel->GetTopPhysicalVolume()->GetName() << G4endl;
  }
#endif
}

void G4VtkSceneHandler::AddCompound(const G4Mesh& mesh)
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkSceneHandler::AddCompound" << G4endl;
#endif 

  if(mesh.GetMeshType() != G4Mesh::rectangle || mesh.GetMeshDepth() != 3)
  {
    G4VSceneHandler::AddCompound(mesh);
  }

  auto container = mesh.GetContainerVolume();
  auto rep1 = container->GetLogicalVolume()->GetDaughter(0);
  auto rep2 = rep1->GetLogicalVolume()->GetDaughter(0);
  auto rep3 = rep2->GetLogicalVolume()->GetDaughter(0);

  EAxis    rep1_axis, rep2_axis, rep3_axis;
  G4int    rep1_nReplicas, rep2_nReplicas, rep3_nReplicas;
  G4double rep1_width, rep2_width, rep3_width;
  G4double rep1_offset, rep2_offset, rep3_offset;
  G4bool   rep1_consuming, rep2_consuming, rep3_consuming;

  rep1->GetReplicationData(rep1_axis,rep1_nReplicas, rep1_width, rep1_offset, rep1_consuming);
  rep2->GetReplicationData(rep2_axis,rep2_nReplicas, rep2_width, rep2_offset, rep2_consuming);
  rep3->GetReplicationData(rep3_axis,rep3_nReplicas, rep3_width, rep3_offset, rep3_consuming);

  // Instantiate a temporary G4PhysicalVolumeModel
  G4ModelingParameters tmpMP;
  tmpMP.SetCulling(false);            // This avoids drawing transparent...
  tmpMP.SetCullingInvisible(false);   // ... or invisble volumes.
  const G4bool useFullExtent = false;  // To avoid calculating the extent

  G4PhysicalVolumeModel tmpPVModel(container, G4PhysicalVolumeModel::UNLIMITED,
                                   G4Transform3D(), &tmpMP, useFullExtent);

  // Instantiate a pseudo scene so that we can make a "private" descent and fill vtkImageData
  vtkSmartPointer<vtkImageData> imagedata = vtkSmartPointer<vtkImageData>::New();
  imagedata->SetDimensions(rep1_nReplicas+1,rep2_nReplicas+1,rep3_nReplicas+1);
  imagedata->SetSpacing(rep1_width,rep2_width,rep2_width);
  imagedata->SetOrigin(0,0,0);
  imagedata->AllocateScalars(VTK_DOUBLE, 1);

  G4double halfX = 0., halfY = 0., halfZ = 0.;
  struct PseudoScene : public G4PseudoScene
  {
    PseudoScene(
      G4PhysicalVolumeModel* pvModel,  // input...the following are output
      vtkImageData *vtkId,
      G4int nx, G4int ny, G4int nz,
      G4double& halfX, G4double& halfY, G4double& halfZ) : fpPVModel(pvModel)
      , fVtkId(vtkId)
      , fNx(nx)
      , fNy(ny)
      , fNz(nz)
      , fHalfX(halfX)
      , fHalfY(halfY)
      , fHalfZ(halfZ)
    {}
    using G4PseudoScene::AddSolid;  // except for...
    void AddSolid(const G4Box& box)
    {
      const G4Colour& colour = fpPVModel->GetCurrentLV()->GetVisAttributes()->GetColour();
      // G4double density = fpPVModel->GetCurrentLV()->GetMaterial()->GetDensity();
      const G4ThreeVector& position = fpCurrentObjectTransformation->getTranslation();
      fHalfX = box.GetXHalfLength();
      fHalfY = box.GetYHalfLength();
      fHalfZ = box.GetZHalfLength();
      fVtkId->SetScalarComponentFromDouble(int(position.x()/(2*fHalfX))+fNx,
                                           int(position.y()/(2*fHalfY))+fNy,
                                           int(position.z()/(2*fHalfZ))+fNz,
                                           0,
                                           (colour.GetRed() + colour.GetBlue() + colour.GetGreen())/3.0);
    }
    G4PhysicalVolumeModel* fpPVModel;
    vtkImageData *fVtkId;
    G4int fNx, fNy, fNz;
    G4double &fHalfX, &fHalfY, &fHalfZ;
  }
  // construct the pseudoScene
  pseudoScene(&tmpPVModel, imagedata, rep1_nReplicas, rep2_nReplicas, rep3_nReplicas, halfX, halfY, halfZ);

  // Make private descent into the nested parameterisation
  tmpPVModel.DescribeYourselfTo(pseudoScene);

  vtkSmartPointer<vtkOpenGLGPUVolumeRayCastMapper> volumeMapper = vtkSmartPointer<vtkOpenGLGPUVolumeRayCastMapper>::New();
  volumeMapper->SetInputData(imagedata);

  vtkNew<vtkVolume> volume;
  volume->SetMapper(volumeMapper);

  vtkSmartPointer<vtkMatrix4x4> vtkMatrix = vtkSmartPointer<vtkMatrix4x4>::New();

  vtkMatrix->SetElement(0,0,fObjectTransformation.xx());
  vtkMatrix->SetElement(0,1,fObjectTransformation.xy());
  vtkMatrix->SetElement(0,2,fObjectTransformation.xz());

  vtkMatrix->SetElement(1,0,fObjectTransformation.yx());
  vtkMatrix->SetElement(1,1,fObjectTransformation.yy());
  vtkMatrix->SetElement(1,2,fObjectTransformation.yz());

  vtkMatrix->SetElement(2,0,fObjectTransformation.zx());
  vtkMatrix->SetElement(2,1,fObjectTransformation.zy());
  vtkMatrix->SetElement(2,2,fObjectTransformation.zz());

  vtkMatrix->SetElement(3,0, fObjectTransformation.dx());
  vtkMatrix->SetElement(3,1, fObjectTransformation.dy());
  vtkMatrix->SetElement(3,2, fObjectTransformation.dz());
  vtkMatrix->SetElement(3,3, 1);

  volume->SetUserMatrix(vtkMatrix);

  auto *pVtkViewer = dynamic_cast<G4VtkViewer *>(fpViewer);
  pVtkViewer->renderer->AddVolume(volume);
}
