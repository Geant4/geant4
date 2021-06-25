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
// John Allison  6th October 2020

#ifdef G4VIS_BUILD_TOOLSSG_DRIVER

#include "G4ToolsSGSceneHandler.hh"

#include "G4TransportationManager.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Polyhedron.hh"
#include "G4Text.hh"

//#define G4TOOLSSG_DEBUG

#include <tools/sg/separator>
#include <tools/sg/matrix>
#include <tools/sg/rgba>
#include <tools/sg/draw_style>
#include <tools/sg/atb_vertices>
#include <tools/sg/markers>
#ifdef TOOLS_USE_FREETYPE
#include <tools/sg/text_freetype_marker>
#include <tools/sg/strings>
#else
#include <tools/sg/text_hershey_marker>
#endif

#include "G4ToolsSGNode.hh"

#include <utility>

G4int G4ToolsSGSceneHandler::fSceneIdCount = 0;

G4ToolsSGSceneHandler::G4ToolsSGSceneHandler
(G4VGraphicsSystem& system, const G4String& name)
:parent(system, fSceneIdCount++, name)
{
#ifdef G4TOOLSSG_DEBUG
  G4cout << "G4ToolsSGSceneHandler::G4ToolsSGSceneHandler called" << G4endl;
#endif
  fpToolsSGScene = new tools::sg::separator;
  EstablishBaseNodes();
}

G4ToolsSGSceneHandler::~G4ToolsSGSceneHandler()
{
  delete fpToolsSGScene;
}

void G4ToolsSGSceneHandler::EstablishBaseNodes()
{
  fpTransientObjects  = new tools::sg::separator;
  fpToolsSGScene->add(fpTransientObjects);

  fpPersistentObjects = new tools::sg::separator;
  fpToolsSGScene->add(fpPersistentObjects);

  // Physical volume objects for each world hang from POs
  G4TransportationManager* transportationManager = G4TransportationManager::GetTransportationManager ();
  size_t nWorlds = transportationManager->GetNoWorlds();
  std::vector<G4VPhysicalVolume*>::iterator iterWorld = transportationManager->GetWorldsIterator();
  fpPhysicalVolumeObjects.resize(nWorlds);
  for (size_t i = 0; i < nWorlds; ++i, ++iterWorld) {
    G4VPhysicalVolume* _world = (*iterWorld);
    auto entity = new G4ToolsSGNode;
    fpPersistentObjects->add(entity);
    entity->SetPVNodeID(G4PhysicalVolumeModel::G4PhysicalVolumeNodeID(_world));
    fpPhysicalVolumeObjects[i] = entity;
  }
}

tools::sg::separator* G4ToolsSGSceneHandler::GetOrCreateNode()
{ // Retrieve or create a G4ToolsSGNode node suitable for next solid or primitive

  // For time being, avoid errors in MT mode - see G4ToolsSGViewer::SwitchToMasterThread
#ifdef G4MULTITHREADED
  if (!G4Threading::IsMasterThread()) return nullptr;
#endif

  if (fReadyForTransients) {  // All transients hang from this node
    tools::sg::separator* sep = new tools::sg::separator;
    fpTransientObjects->add(sep);
    return sep;
  }

  G4PhysicalVolumeModel* pPVModel = dynamic_cast<G4PhysicalVolumeModel*>(fpModel);

  if (!pPVModel) {  // Persistent objects (e.g., axes)
    tools::sg::separator* sep = new tools::sg::separator;
    fpPersistentObjects->add(sep);
    return sep;
  }

  // So this is a G4PhysicalVolumeModel
  typedef G4PhysicalVolumeModel::G4PhysicalVolumeNodeID PVNodeID;
  typedef std::vector<PVNodeID> PVPath;
  //const PVPath& drawnPVPath = pPVModel->GetDrawnPVPath();
  const PVPath& fullPVPath  = pPVModel->GetFullPVPath();
  //G4int currentDepth = pPVModel->GetCurrentDepth();
  //G4VPhysicalVolume* pCurrentPV = pPVModel->GetCurrentPV();
  //G4LogicalVolume* pCurrentLV = pPVModel->GetCurrentLV();
  //G4Material* pCurrentMaterial = pPVModel->GetCurrentMaterial();
  // Note: pCurrentMaterial may be zero (parallel world).

  // Find appropriate root
  const size_t nWorlds = fpPhysicalVolumeObjects.size();
  size_t iWorld = 0;
  for (; iWorld < nWorlds; ++iWorld) {
    if (fullPVPath[0].GetPhysicalVolume() ==
        fpPhysicalVolumeObjects[iWorld]->GetPVNodeID().GetPhysicalVolume()) break;
  }
  if (iWorld == nWorlds) {
    G4Exception("G4ToolsSGSceneHandler::GetOrCreateNode", "ToolsSG-0000", FatalException,
                "World mis-match - not possible(!?)");
  }

  // (Re-)establish pv path of root entity
  G4ToolsSGNode* _world = fpPhysicalVolumeObjects[iWorld];
  _world->SetPVNodeID(fullPVPath[0]);

  // Provide nodes as required - may be a new node or a pre-existing node
  G4ToolsSGNode* node = _world;  // Working variable - default to world
  const size_t depth = fullPVPath.size();
  size_t iDepth = 1;
  while (iDepth < depth) {
    const auto& children = node->children();
    const G4int nChildren = children.size();
    G4int iChild = 0;
    G4ToolsSGNode* child = nullptr;
    for (; iChild < nChildren; ++iChild) {
      child = static_cast<G4ToolsSGNode*>(children[iChild]);
      if (child->GetPVNodeID() == fullPVPath[iDepth]) break;
    }
    if (iChild != nChildren) {  // Existing node found
      node = child;  // Must be the ancestor of new node (subsequent iteration)
    } else {
      // Add a new node as child of node
      G4ToolsSGNode* newNode = new G4ToolsSGNode;
      node->add(newNode);
      newNode->SetPVNodeID(fullPVPath[iDepth]);
      node = newNode;
    }
    ++iDepth;
  }
  return node;
}

void G4ToolsSGSceneHandler::PreAddSolid
(const G4Transform3D& objectTransformation,
 const G4VisAttributes& visAttribs)
{  
  parent::PreAddSolid(objectTransformation, visAttribs);
}

void G4ToolsSGSceneHandler::PostAddSolid()
{
  parent::PostAddSolid();
}

void G4ToolsSGSceneHandler::BeginPrimitives2D(const G4Transform3D& objectTransformation)
{
  // The x,y coordinates of the primitives passed to AddPrimitive are
  // intrepreted as screen coordinates, -1 < x,y < 1.  The
  // z-coordinate is ignored.
  static G4bool first = true;
  if (first) {
    first = false;
    G4Exception("G4ToolsSGSceneHandler::BeginPrimitives2D", "ToolsSG-0001",
                JustWarning,
                "2D drawing not yet implemented");
  }
  parent::BeginPrimitives2D (objectTransformation);
}

void G4ToolsSGSceneHandler::EndPrimitives2D ()
{
  parent::EndPrimitives2D ();
}

void G4ToolsSGSceneHandler::BeginPrimitives(const G4Transform3D& objectTransformation)
{
  parent::BeginPrimitives(objectTransformation);
}

void G4ToolsSGSceneHandler::EndPrimitives ()
{
  parent::EndPrimitives ();
}

void G4ToolsSGSceneHandler::AddPrimitive(const G4Polyline& a_polyline)
{
#ifdef G4TOOLSSG_DEBUG
  G4cout << "debug : G4ToolsSGSceneHandler::AddPrimitive(const G4Polyline&) : \n" << a_polyline << G4endl;
#endif
  if (a_polyline.size() == 0) return;

  auto currentNode = GetOrCreateNode();
  if (!currentNode) return;  // Node not available

  fpVisAttribs = fpViewer->GetApplicableVisAttributes(a_polyline.GetVisAttributes());

  // Transformation
 {tools::sg::matrix* mtx = new tools::sg::matrix;
  G4Transform3D& elem = fObjectTransformation;
  mtx->mtx.value().set_matrix(elem(0,0),elem(0,1),elem(0,2),elem(0,3),
                              elem(1,0),elem(1,1),elem(1,2),elem(1,3),
                              elem(2,0),elem(2,1),elem(2,2),elem(2,3),
                                      0,        0,        0,        1);
  currentNode->add(mtx);}

  // Material (colour, etc.)
 {const auto& colour = GetColour();
  tools::sg::rgba* mat = new tools::sg::rgba();
  mat->color =
    tools::colorf(float(colour.GetRed()),
                  float(colour.GetGreen()),
                  float(colour.GetBlue()),
                  float(colour.GetAlpha()));
  currentNode->add(mat);}

 {tools::sg::draw_style* ds = new tools::sg::draw_style;
  ds->style = tools::sg::draw_lines;
  ds->line_width = 1;
  currentNode->add(ds);}

  // Geometry :
  tools::sg::vertices* vtxs = new tools::sg::vertices;
  vtxs->mode = tools::gl::line_strip();  //polyline
  currentNode->add(vtxs);
  
 {for (size_t i = 0; i < a_polyline.size(); ++i) {
    vtxs->add(float(a_polyline[i].x()),float(a_polyline[i].y()),float(a_polyline[i].z()));
  }}
  
}

void G4ToolsSGSceneHandler::AddPrimitive (const G4Scale& a_scale) {parent::AddPrimitive(a_scale);}

void G4ToolsSGSceneHandler::AddPrimitive (const G4Polymarker& a_polymarker)
{
#ifdef G4TOOLSSG_DEBUG
  ::printf("debug G4ToolsSGSceneHandler::AddPrimitive(const G4Polymarker&) : %lu, type %d\n",
	   a_polymarker.size(),a_polymarker.GetMarkerType());
#endif
  if (a_polymarker.size() == 0) return;
  auto currentNode = GetOrCreateNode();
  if (!currentNode) return;  // Node not available

  fpVisAttribs = fpViewer->GetApplicableVisAttributes(a_polymarker.GetVisAttributes());

  // Transformation
 {tools::sg::matrix* mtx = new tools::sg::matrix;
  G4Transform3D& elem = fObjectTransformation;
  mtx->mtx.value().set_matrix(elem(0,0),elem(0,1),elem(0,2),elem(0,3),
                              elem(1,0),elem(1,1),elem(1,2),elem(1,3),
                              elem(2,0),elem(2,1),elem(2,2),elem(2,3),
                                      0,        0,        0,        1);
  currentNode->add(mtx);}

  // Material (colour, etc.)
 {const auto& colour = GetColour();
  tools::sg::rgba* mat = new tools::sg::rgba();
  mat->color =
    tools::colorf(float(colour.GetRed()),
                  float(colour.GetGreen()),
                  float(colour.GetBlue()),
                  float(colour.GetAlpha()));
  currentNode->add(mat);}

  switch (a_polymarker.GetMarkerType()) {
    default:
    case G4Polymarker::dots:{
      //::printf("debug : GB : Add Markers : +++++++++++++++++++++++++++++++++++++++++++ : dots\n");
      tools::sg::draw_style* ds = new tools::sg::draw_style;
      ds->style = tools::sg::draw_points;
      ds->point_size = 10;
      currentNode->add(ds);

      tools::sg::vertices* vtxs = new tools::sg::vertices;
      vtxs->mode = tools::gl::points();
     {for (size_t i = 0; i < a_polymarker.size(); ++i) {
        vtxs->add(float(a_polymarker[i].x()),float(a_polymarker[i].y()),float(a_polymarker[i].z()));
      }}
      currentNode->add(vtxs);
    }break;
    case G4Polymarker::circles:{
      //::printf("debug : GB : Add Markers : +++++++++++++++++++++++++++++++++++++++++++ : circles\n");
     {tools::sg::markers* markers = new tools::sg::markers;
      markers->size = 10;
      markers->style = tools::sg::marker_circle_line;
      for (size_t i = 0; i < a_polymarker.size(); ++i) {
        markers->add(float(a_polymarker[i].x()),float(a_polymarker[i].y()),float(a_polymarker[i].z()));
      }
      currentNode->add(markers);}
     /*
      G4Circle circle (a_polymarker);
      G4double radius;
      if (circle.GetSizeType() == G4VMarker::world ) {
        radius = circle.GetWorldRadius();
      } else {  // Screen-size or none
        // Not figured out how to do screen-size, so use scene extent
        const G4double scale = 200.;  // Roughly pixles per scene
        radius = circle.GetScreenRadius()*fpScene->GetExtent().GetExtentRadius()/scale;
      }
      for (size_t iPoint = 0; iPoint < a_polymarker.size (); iPoint++) {
        auto position = fObjectTransformation*G4Translate3D(a_polymarker[iPoint]);
        auto node = new G4ToolsSGNode(currentNode);
        //node->addComponent(material);
        //node->addComponent(position);
        //node->addComponent(sphere);  // Sphere of given radius
      }
      */
    }break;
  case G4Polymarker::squares:{
    //::printf("debug : GB : Add Markers : +++++++++++++++++++++++++++++++++++++++++++ : square\n");
     {tools::sg::markers* markers = new tools::sg::markers;
      markers->size = 10;
      markers->style = tools::sg::marker_square_line;
      for (size_t i = 0; i < a_polymarker.size(); ++i) {
        markers->add(float(a_polymarker[i].x()),float(a_polymarker[i].y()),float(a_polymarker[i].z()));
      }
      currentNode->add(markers);}
      /*
      G4Square square (a_polymarker);
      G4double side;
      if (square.GetSizeType() == G4VMarker::world ) {
        side = square.GetWorldDiameter();
      } else {  // Screen-size or none
        // Not figured out how to do screen-size, so use scene extent
        const G4double scale = 200.;  // Roughly pixles per scene
        side = square.GetScreenDiameter()*fpScene->GetExtent().GetExtentRadius()/scale;
      }
      for (size_t iPoint = 0; iPoint < a_polymarker.size (); iPoint++) {
        auto position = fObjectTransformation*G4Translate3D(a_polymarker[iPoint]);
        auto node = new G4ToolsSGNode(currentNode);
        //node->addComponent(material);
        //node->addComponent(position);
        //node->addComponent(box);  // Cube of given side
      }
      */
  }break;
  }
}

void G4ToolsSGSceneHandler::AddPrimitive(const G4Text& a_text)
{
  //static G4bool first = true;
  //if (first) {
  //  first = false;
  //  G4Exception("G4ToolsSGSceneHandler::AddPrimitive(const G4Text& text)",
  //              "ToolsSG-0002", JustWarning,
  //              "Text drawing not yet implemented");
  //}
#ifdef G4TOOLSSG_DEBUG
  ::printf("debug : G4ToolsSGSceneHandler::AddPrimitive(const G4Text&) : \"%s\"\n",a_text.GetText().c_str());
#endif
  
  auto currentNode = GetOrCreateNode();
  if (!currentNode) return;  // Node not available

  fpVisAttribs = fpViewer->GetApplicableVisAttributes(a_text.GetVisAttributes());

  auto pos = a_text.GetPosition();
  //::printf("debug : Add Text : pos %g %g %g\n",pos.x(),pos.y(),pos.z());

  // Transformation
 {tools::sg::matrix* mtx = new tools::sg::matrix;
  auto elem = fObjectTransformation*G4Translate3D(a_text.GetPosition());
  //G4Transform3D& elem = fObjectTransformation*G4Translate3D(a_text.GetPosition());
  mtx->mtx.value().set_matrix(elem(0,0),elem(0,1),elem(0,2),elem(0,3),
                              elem(1,0),elem(1,1),elem(1,2),elem(1,3),
                              elem(2,0),elem(2,1),elem(2,2),elem(2,3),
                                      0,        0,        0,        1);
  currentNode->add(mtx);}

  // Material (colour, etc.)
 {const auto& colour = GetColour();
  tools::sg::rgba* mat = new tools::sg::rgba();
  mat->color =
    tools::colorf(float(colour.GetRed()),
                  float(colour.GetGreen()),
                  float(colour.GetBlue()),
                  float(colour.GetAlpha()));
  currentNode->add(mat);}

#ifdef TOOLS_USE_FREETYPE
  tools::sg::text_freetype_marker* text = new tools::sg::text_freetype_marker;
  text->font = tools::sg::font_lato_regular_ttf();
  text->front_face = tools::sg::winding_cw;
  text->modeling = tools::sg::font_pixmap;
#else
  tools::sg::text_hershey_marker* text = new tools::sg::text_hershey_marker;
//text->encoding.value(a_encoding);
#endif
  text->strings.add(a_text.GetText());
 {switch (a_text.GetLayout()) {
  default:
  case G4Text::left:
    text->hjust = tools::sg::left;
    break;
  case G4Text::centre:
    text->hjust = tools::sg::center;
    break;
  case G4Text::right:
    text->hjust = tools::sg::right;
    break;
  }}
//text->vjust.value(a_vjust);
  currentNode->add(text);
 
  /*
  //text.GetText(), text.GetXOffset(), etc.

  MarkerSizeType sizeType;
  G4double size = GetMarkerSize (text, sizeType);
  switch (sizeType) {
  default:
  case screen:
    // Draw in screen coordinates.  OK.
    break;
  case world:
    // Draw in world coordinates.   Not implemented.  Use size = 20.
    size = 20.;
    break;
  }
  */
}

void G4ToolsSGSceneHandler::AddPrimitive(const G4Circle& a_circle)
{
  G4Polymarker oneCircle(a_circle);
  oneCircle.push_back(a_circle.GetPosition());
  oneCircle.SetMarkerType(G4Polymarker::circles);
  // Call this AddPrimitive to avoid re-doing sub-class code.
  G4ToolsSGSceneHandler::AddPrimitive(oneCircle);
    
  /*
  auto currentNode = GetOrCreateNode();
  if (!currentNode) return;  // Node not available

  fpVisAttribs = fpViewer->GetApplicableVisAttributes(circle.GetVisAttributes());

  auto circleNode = new G4ToolsSGNode(currentNode);

  // Transformation
  auto position = fObjectTransformation*G4Translate3D(circle.GetPosition());
  //circleNode->addComponent(position);

  // Material (colour, etc.)
  const auto& colour = fpVisAttribs->GetColour();
  //circleNode->addComponent(colour);

  G4double radius;
  if (circle.GetSizeType() == G4VMarker::world ) {
    radius =circle.GetWorldRadius();
  } else {  // Screen-size or none
    // Not figured out how to do screen-size, so use scene extent
    const G4double scale = 200.;  // Roughly pixles per scene
    radius = circle.GetScreenRadius()*fpScene->GetExtent().GetExtentRadius()/scale;
  }
  //circleNode->addComponent(sphere);  // Of given radius
  */
}

void G4ToolsSGSceneHandler::AddPrimitive(const G4Square& a_square)
{
  G4Polymarker oneSquare(a_square);
  oneSquare.push_back(a_square.GetPosition());
  oneSquare.SetMarkerType(G4Polymarker::squares);
  // Call this AddPrimitive to avoid re-doing sub-class code.
  G4ToolsSGSceneHandler::AddPrimitive(oneSquare);

  /*
  auto currentNode = GetOrCreateNode();
  if (!currentNode) return;  // Node not available

  fpVisAttribs = fpViewer->GetApplicableVisAttributes(square.GetVisAttributes());

  auto squareNode = new G4ToolsSGNode(currentNode);

  // Transformation
  auto position = fObjectTransformation*G4Translate3D(square.GetPosition());
  //squareNode->addComponent(position);

  // Material (colour, etc.)
  const auto& colour = fpVisAttribs->GetColour();
  //squareNode->addComponent(colour);

  G4double side;
  if (square.GetSizeType() == G4VMarker::world ) {
    side = square.GetWorldDiameter();
  } else {  // Screen-size or none
    // Not figured out how to do screen-size, so use scene extent
    const G4double scale = 200.;  // Roughly pixles per scene
    side = square.GetScreenDiameter()*fpScene->GetExtent().GetExtentRadius()/scale;
  }
  //squareNode->addComponent(cube);  // Of given side
  */
}

void G4ToolsSGSceneHandler::AddPrimitive(const G4Polyhedron& a_polyhedron)
{
  if (a_polyhedron.GetNoFacets() == 0) return;

#ifdef G4TOOLSSG_DEBUG
  ::printf("debug : G4ToolsSGSceneHandler::AddPrimitive(const G4Polyhedron&) : %d\n",a_polyhedron.GetNoFacets());
#endif
  
  fpVisAttribs = fpViewer->GetApplicableVisAttributes(a_polyhedron.GetVisAttributes());

  // Roll out vertices and normals for the faces. Note that this means vertices
  // are duplicated. For example a box has 8 vertices, but to define 6 faces
  // you need 12 triangles and 36 vertices. If it was just a matter of vertices
  // we could restrict the number to 8 and use the indices to define the
  // triangles, but we also have to consider the normals. A vertex can be have
  // more than one normal, depending on which face it is being used to define.
  // So we roll out all the vertices and normals for each triangle.
  std::vector<G4Point3D> vertices;
  std::vector<G4Normal3D> normals;

  // Also roll out edges (as lines) for wireframe. Avoid duplicate lines,
  // including those that differ only in the order of vertices.
  typedef std::pair<G4Point3D,G4Point3D> Line;
  std::vector<Line> lines;
  auto insertIfNew = [&lines](const Line& newLine) {
    for (const auto& line: lines) {
      if ((newLine.first==line.first && newLine.second==line.second) ||
          (newLine.first==line.second && newLine.second==line.first))
      return;
    }
    lines.push_back(newLine);
  };

  G4bool isAuxilaryEdgeVisible = fpViewer->GetViewParameters().IsAuxEdgeVisible();
  G4bool notLastFace;
  do {
    G4int      nEdges;
    G4Point3D  vertex  [4];
    G4int      edgeFlag[4];
    G4Normal3D normal  [4];
    notLastFace = a_polyhedron.GetNextFacet(nEdges, vertex, edgeFlag, normal);
    vertices.push_back(vertex[0]);
    vertices.push_back(vertex[1]);
    vertices.push_back(vertex[2]);
    normals.push_back(normal[0]);
    normals.push_back(normal[1]);
    normals.push_back(normal[2]);
    if(isAuxilaryEdgeVisible||edgeFlag[0]>0)insertIfNew(Line(vertex[0],vertex[1]));
    if(isAuxilaryEdgeVisible||edgeFlag[1]>0)insertIfNew(Line(vertex[1],vertex[2]));
    if (nEdges == 3) {
      // Face is a triangle
      // One more line for wireframe, triangles for surfaces are complete
      if(isAuxilaryEdgeVisible||edgeFlag[2]>0)insertIfNew(Line(vertex[2],vertex[0]));
    } else if (nEdges == 4) {
      // Face is a quadrilateral
      // Create another triangle for surfaces, add two more lines for wireframe
      vertices.push_back(vertex[2]);
      vertices.push_back(vertex[3]);
      vertices.push_back(vertex[0]);
      normals.push_back(normal[2]);
      normals.push_back(normal[3]);
      normals.push_back(normal[0]);
      if(isAuxilaryEdgeVisible||edgeFlag[2]>0)insertIfNew(Line(vertex[2],vertex[3]));
      if(isAuxilaryEdgeVisible||edgeFlag[3]>0)insertIfNew(Line(vertex[3],vertex[0]));
    } else {
      G4cerr
      << "ERROR: polyhedron face with unexpected number of edges (" << nEdges << ')'
      << "\n  Tag: " << fpModel->GetCurrentTag()
      << G4endl;
      return;
    }
  } while (notLastFace);

  G4ViewParameters::DrawingStyle drawing_style = GetDrawingStyle (fpVisAttribs);
  switch (drawing_style) {
    case G4ViewParameters::wireframe:
      //vertices.clear();
      break;
    case G4ViewParameters::hlr:
      break;
    case G4ViewParameters::hsr:
      //lines.clear();
      break;
    case G4ViewParameters::hlhsr:
      break;
    case G4ViewParameters::cloud:
      // Shouldn't happen in this function (it's a polyhedron!) - ignore
      return;
  }

  auto currentNode = GetOrCreateNode();
  if (!currentNode) return;  // Node not available
  
  tools::sg::separator* sep = new tools::sg::separator;
  currentNode->add(sep);

  // Transformation
 {tools::sg::matrix* mtx = new tools::sg::matrix;
  G4Transform3D& elem = fObjectTransformation;
  mtx->mtx.value().set_matrix(elem(0,0),elem(0,1),elem(0,2),elem(0,3),
                              elem(1,0),elem(1,1),elem(1,2),elem(1,3),
                              elem(2,0),elem(2,1),elem(2,2),elem(2,3),
                                      0,        0,        0,        1);
  sep->add(mtx);}

  // Material (colour, etc.)
 {const auto& colour = fpVisAttribs->GetColour();
  tools::sg::rgba* mat = new tools::sg::rgba();
  mat->color =
    tools::colorf(float(colour.GetRed()),
                  float(colour.GetGreen()),
                  float(colour.GetBlue()),
                  float(colour.GetAlpha()));
  sep->add(mat);}

  if (drawing_style == G4ViewParameters::hlr ||
      drawing_style == G4ViewParameters::hsr ||
      drawing_style == G4ViewParameters::hlhsr) {

   {tools::sg::draw_style* ds = new tools::sg::draw_style;
    ds->style = tools::sg::draw_filled;
  //ds->cull_face = true;
    sep->add(ds);}

    tools::sg::atb_vertices* vtxs = new tools::sg::atb_vertices;
    vtxs->mode = tools::gl::triangles();
    sep->add(vtxs);
    
    const auto nVerts = vertices.size();
    for (size_t i = 0; i < nVerts; i++) {
      vtxs->add(float(vertices[i].x()),float(vertices[i].y()),float(vertices[i].z()));
      vtxs->add_normal(float(normals[i].x()),float(normals[i].y()),float(normals[i].z()));
    }
  }

  if (drawing_style == G4ViewParameters::wireframe ||
      drawing_style == G4ViewParameters::hlr ||
      drawing_style == G4ViewParameters::hlhsr) {

   {tools::sg::draw_style* ds = new tools::sg::draw_style;
    ds->style = tools::sg::draw_lines;
    ds->line_width = 1;
    sep->add(ds);}

    tools::sg::vertices* vtxs = new tools::sg::vertices;
    vtxs->mode = tools::gl::lines();  //segments
    sep->add(vtxs);
    
    for (const auto& line: lines) {
      vtxs->add(float(line.first.x()),float(line.first.y()),float(line.first.z()));
      vtxs->add(float(line.second.x()),float(line.second.y()),float(line.second.z()));
    }
  
  }
}

void G4ToolsSGSceneHandler::ClearStore ()
{
  fpToolsSGScene->clear();
  EstablishBaseNodes();
}

void G4ToolsSGSceneHandler::ClearTransientStore ()
{
  fpTransientObjects->clear();
}

#endif
