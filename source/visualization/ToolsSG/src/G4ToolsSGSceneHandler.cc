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

#include "G4ToolsSGSceneHandler.hh"

#include "G4ToolsSGNode.hh"

#include "G4TransportationManager.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Polyhedron.hh"
#include "G4Text.hh"
#include "G4Mesh.hh"
#include "G4PlotterManager.hh"

#include <tools/sg/separator>
#include <tools/sg/matrix>
#include <tools/sg/rgba>
#include <tools/sg/draw_style>
#include <tools/sg/atb_vertices>
#include <tools/sg/markers>
#ifdef TOOLS_USE_FREETYPE
#include <toolx/sg/text_freetype>
#include <tools/sg/strings>
#include <tools/font/lato_regular_ttf>
#include <tools/font/roboto_bold_ttf>
#include <toolx/sg/text_freetype_marker>
#else
#include <tools/sg/dummy_freetype>
#include <tools/sg/text_hershey_marker>
#endif

//for plotting:
#include <tools/sg/dummy_freetype>
#include <tools/sg/light_off>
#include <tools/sg/plots>
#include <tools/sg/h2plot_cp>
#include <tools/sg/plotter_style>
#include <tools/sg/event_dispatcher>
#include <tools/sg/path>
#include <tools/sg/search>
#include <tools/histo/h1d>
#include <tools/histo/h2d>
#include <tools/sg/plotter_some_styles>

#include <utility>

G4int G4ToolsSGSceneHandler::fSceneIdCount = 0;

G4ToolsSGSceneHandler::G4ToolsSGSceneHandler
(G4VGraphicsSystem& system, const G4String& name)
:parent(system, fSceneIdCount++, name)
,fFreetypeNode(0)
{
  //::printf("debug : G4ToolsSGSceneHandler : %lu, %s\n",this,name.c_str());
  EstablishBaseNodes();
#if defined(TOOLS_USE_FREETYPE)
  fFreetypeNode = new toolx::sg::text_freetype();
  fFreetypeNode->add_embedded_font(tools::sg::font_lato_regular_ttf(),tools::font::lato_regular_ttf);
  fFreetypeNode->add_embedded_font(tools::sg::font_roboto_bold_ttf(),tools::font::roboto_bold_ttf);
#else  
  fFreetypeNode = new tools::sg::dummy_freetype();
#endif
  Messenger::Create();
}

G4ToolsSGSceneHandler::~G4ToolsSGSceneHandler()
{
  //::printf("debug : ~G4ToolsSGSceneHandler : %lu\n",this);
  //WARNING : nodes may refer graphics managers (as tools/sg/[GL_manager,gl2ps_manager,zb_manager]
  //          used by viewers) to handle gstos (for GPU) or textures, then we have to delete them first.
  //          It is assumed that we pass here BEFORE the attached/managed viewers are deleted.
  fpTransient2DObjects.clear();
  fpPersistent2DObjects.clear();
  fpTransient3DObjects.clear();
  fpPersistent3DObjects.clear();
  delete fFreetypeNode;
}

void G4ToolsSGSceneHandler::EstablishBaseNodes()
{
  // Physical volume objects for each world hang from POs
  G4TransportationManager* transportationManager = G4TransportationManager::GetTransportationManager ();
  size_t nWorlds = transportationManager->GetNoWorlds();
  std::vector<G4VPhysicalVolume*>::iterator iterWorld = transportationManager->GetWorldsIterator();
  fpPhysicalVolumeObjects.resize(nWorlds);
  for (size_t i = 0; i < nWorlds; ++i, ++iterWorld) {
    G4VPhysicalVolume* _world = (*iterWorld);
    auto entity = new G4ToolsSGNode;
    fpPersistent3DObjects.add(entity);
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
    fpTransient3DObjects.add(sep);
    return sep;
  }

  auto* pPVModel = dynamic_cast<G4PhysicalVolumeModel*>(fpModel);

  if (!pPVModel) {  // Persistent objects (e.g., axes)
    tools::sg::separator* sep = new tools::sg::separator;
    fpPersistent3DObjects.add(sep);
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
    const G4int nChildren = (G4int)children.size();
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

void G4ToolsSGSceneHandler::ClearStore ()
{
  fpTransient2DObjects.clear();
  fpPersistent2DObjects.clear();
  fpTransient3DObjects.clear();
  fpPersistent3DObjects.clear();
  EstablishBaseNodes();
}

void G4ToolsSGSceneHandler::ClearTransientStore ()
{
  fpTransient2DObjects.clear();
  fpTransient3DObjects.clear();
}

void G4ToolsSGSceneHandler::AddPrimitive(const G4Polyline& a_polyline)
{
  //G4cout << "debug : G4ToolsSGSceneHandler::AddPrimitive(const G4Polyline&) : \n" << a_polyline << G4endl;
  if (a_polyline.size() == 0) return;

  tools::sg::separator* parentNode = 0;
  if(fProcessing2D) {
    parentNode = new tools::sg::separator;
    if (fReadyForTransients) {
      fpTransient2DObjects.add(parentNode);
    } else {
      fpPersistent2DObjects.add(parentNode);
    }

  } else {
    parentNode = GetOrCreateNode();
    if(!parentNode) return;

    tools::sg::matrix* mtx = new tools::sg::matrix;
    G4Transform3D& elem = fObjectTransformation;
    mtx->mtx.value().set_matrix(elem(0,0),elem(0,1),elem(0,2),elem(0,3),
                                elem(1,0),elem(1,1),elem(1,2),elem(1,3),
                                elem(2,0),elem(2,1),elem(2,2),elem(2,3),
                                        0,        0,        0,        1);
    parentNode->add(mtx);
  }

 {const auto& colour = GetColour(a_polyline);
  tools::sg::rgba* mat = new tools::sg::rgba();
  mat->color =
    tools::colorf(float(colour.GetRed()),
                  float(colour.GetGreen()),
                  float(colour.GetBlue()),
                  float(colour.GetAlpha()));
  parentNode->add(mat);}

 {tools::sg::draw_style* ds = new tools::sg::draw_style;
  ds->style = tools::sg::draw_lines;
  ds->line_width = 1;
  parentNode->add(ds);}

  tools::sg::vertices* vtxs = new tools::sg::vertices;
  vtxs->mode = tools::gl::line_strip();  //polyline
  parentNode->add(vtxs);
  
 {for (size_t i = 0; i < a_polyline.size(); ++i) {
    vtxs->add(float(a_polyline[i].x()),float(a_polyline[i].y()),float(a_polyline[i].z()));
  }}
  
}

void G4ToolsSGSceneHandler::AddPrimitive (const G4Polymarker& a_polymarker)
{
  //::printf("debug G4ToolsSGSceneHandler::AddPrimitive(const G4Polymarker&) : %lu, type %d\n",
  //	   a_polymarker.size(),a_polymarker.GetMarkerType());
  if (a_polymarker.size() == 0) return;
  auto currentNode = GetOrCreateNode();
  if (!currentNode) return;  // Node not available

  // Transformation
 {tools::sg::matrix* mtx = new tools::sg::matrix;
  G4Transform3D& elem = fObjectTransformation;
  mtx->mtx.value().set_matrix(elem(0,0),elem(0,1),elem(0,2),elem(0,3),
                              elem(1,0),elem(1,1),elem(1,2),elem(1,3),
                              elem(2,0),elem(2,1),elem(2,2),elem(2,3),
                                      0,        0,        0,        1);
  currentNode->add(mtx);}

 {const auto& colour = GetColour(a_polymarker);
  tools::sg::rgba* mat = new tools::sg::rgba();
  mat->color =
    tools::colorf(float(colour.GetRed()),
                  float(colour.GetGreen()),
                  float(colour.GetBlue()),
                  float(colour.GetAlpha()));
  currentNode->add(mat);}

  MarkerSizeType markerSizeType;
  G4double markerSize = GetMarkerSize(a_polymarker, markerSizeType);

  switch (a_polymarker.GetMarkerType()) {
    default:
    case G4Polymarker::dots:{
      //::printf("debug : GB : Add Markers : +++++++++++++++++++++++++++++++++++++++++++ : dots\n");
      tools::sg::draw_style* ds = new tools::sg::draw_style;
      ds->style = tools::sg::draw_points;
      ds->point_size = 1;
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
      G4double diameter = markerSize;  // OK for "screen-size" (the usual case)
      if (markerSizeType == G4VSceneHandler::world ) {
        const G4double scale = 200.;  // Roughly pixels per scene
        diameter *= fpScene->GetExtent().GetExtentRadius()/scale;
      }
      markers->size = diameter;
      markers->style = tools::sg::marker_circle_line;
      for (size_t i = 0; i < a_polymarker.size(); ++i) {
        markers->add(float(a_polymarker[i].x()),float(a_polymarker[i].y()),float(a_polymarker[i].z()));
      }
      currentNode->add(markers);}
    }break;
  case G4Polymarker::squares:{
    //::printf("debug : GB : Add Markers : +++++++++++++++++++++++++++++++++++++++++++ : square\n");
     {tools::sg::markers* markers = new tools::sg::markers;
      G4double side = markerSize;  // OK for "screen-size" (the usual case)
      if (markerSizeType == G4VSceneHandler::world ) {
        const G4double scale = 200.;  // Roughly pixels per scene
        side *= fpScene->GetExtent().GetExtentRadius()/scale;
      }
      markers->size = side;
      markers->style = tools::sg::marker_square_line;
      for (size_t i = 0; i < a_polymarker.size(); ++i) {
        markers->add(float(a_polymarker[i].x()),float(a_polymarker[i].y()),float(a_polymarker[i].z()));
      }
      currentNode->add(markers);}
  }break;
  }
}

void G4ToolsSGSceneHandler::AddPrimitive(const G4Text& a_text)
{
  //::printf("debug : G4ToolsSGSceneHandler::AddPrimitive(const G4Text&) : 000 : \"%s\"\n",a_text.GetText().c_str());
  //::printf("debug : G4ToolsSGSceneHandler::AddPrimitive(const G4Text&) : 2D ? %d\n",fProcessing2D);
  auto pos = a_text.GetPosition();
  //::printf("debug : Add Text : pos %g %g %g\n",pos.x(),pos.y(),pos.z());

  tools::sg::separator* parentNode = 0;
  if(fProcessing2D) {
    parentNode = new tools::sg::separator;
    if (fReadyForTransients) {
      fpTransient2DObjects.add(parentNode);
    } else {
      fpPersistent2DObjects.add(parentNode);
    }

    tools::sg::matrix* mtx = new tools::sg::matrix;
    mtx->set_translate(pos.x(),pos.y(),pos.z());
    parentNode->add(mtx);

  } else {
    parentNode = GetOrCreateNode();
    if (!parentNode) return;

    tools::sg::matrix* mtx = new tools::sg::matrix;
    auto elem = fObjectTransformation*G4Translate3D(pos);
    mtx->mtx.value().set_matrix(elem(0,0),elem(0,1),elem(0,2),elem(0,3),
                                elem(1,0),elem(1,1),elem(1,2),elem(1,3),
                                elem(2,0),elem(2,1),elem(2,2),elem(2,3),
                                        0,        0,        0,        1);
    parentNode->add(mtx);
  }

  MarkerSizeType sizeType;
  G4double size = GetMarkerSize(a_text, sizeType);
  
 {const auto& colour = GetTextColour(a_text);
  tools::sg::rgba* mat = new tools::sg::rgba();
  mat->color =
    tools::colorf(float(colour.GetRed()),
                  float(colour.GetGreen()),
                  float(colour.GetBlue()),
                  float(colour.GetAlpha()));
  parentNode->add(mat);}
 
#ifdef TOOLS_USE_FREETYPE
  toolx::sg::text_freetype_marker* text = new toolx::sg::text_freetype_marker;
  text->add_embedded_font(tools::sg::font_lato_regular_ttf(),tools::font::lato_regular_ttf);
  text->font = tools::sg::font_lato_regular_ttf();
  text->front_face = tools::sg::winding_cw;
//text->modeling = tools::sg::font_pixmap; //problem with Qt/GL. It slows rendering!
#else
  tools::sg::text_hershey_marker* text = new tools::sg::text_hershey_marker;
//text->encoding.value(a_encoding);
#endif
  text->height = float(size); //pixels
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
  parentNode->add(text);

}

void G4ToolsSGSceneHandler::AddPrimitive(const G4Circle& a_circle)
{
  G4Polymarker oneCircle(a_circle);
  oneCircle.push_back(a_circle.GetPosition());
  oneCircle.SetMarkerType(G4Polymarker::circles);
  // Call this AddPrimitive to avoid re-doing sub-class code.
  G4ToolsSGSceneHandler::AddPrimitive(oneCircle);
}

void G4ToolsSGSceneHandler::AddPrimitive(const G4Square& a_square)
{
  G4Polymarker oneSquare(a_square);
  oneSquare.push_back(a_square.GetPosition());
  oneSquare.SetMarkerType(G4Polymarker::squares);
  // Call this AddPrimitive to avoid re-doing sub-class code.
  G4ToolsSGSceneHandler::AddPrimitive(oneSquare);
}

void G4ToolsSGSceneHandler::AddPrimitive(const G4Polyhedron& a_polyhedron)
{
  if (a_polyhedron.GetNoFacets() == 0) return;

  //::printf("debug : G4ToolsSGSceneHandler::AddPrimitive(const G4Polyhedron&) : %d\n",a_polyhedron.GetNoFacets());
  
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
//    for (const auto& line: lines) {
//      if ((newLine.first==line.first && newLine.second==line.second) ||
//          (newLine.first==line.second && newLine.second==line.first))
//      return;
//    }
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

 {const auto& colour = GetColour(a_polyhedron);
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

void G4ToolsSGSceneHandler::AddCompound(const G4Mesh& mesh)
{
  StandardSpecialMeshRendering(mesh);
}

//plotting:
inline void SetRegionStyles(tools::xml::styles& a_styles,
		            tools::sg::plots& a_plots,
		            tools::sg::plotter& a_plotter,
                            const G4String& a_style) {
  if(a_style=="reset") {
    a_plotter.reset_style(true);
    a_plots.touch(); //to apply indirectly plots::set_plotter_layout() on _plotter.
  } else if( (a_style=="inlib_default")|| (a_style=="default")) {
    tools::sg::set_inlib_default_style(G4cout,a_styles.cmaps(),a_plotter,tools::sg::font_hershey());
  } else if(a_style=="ROOT_default") {
    tools::sg::set_ROOT_default_style(G4cout,a_styles.cmaps(),a_plotter,tools::sg::font_roboto_bold_ttf());
  } else if(a_style=="hippodraw") {
    tools::sg::set_hippodraw_style(G4cout,a_styles.cmaps(),a_plotter,tools::sg::font_lato_regular_ttf());
  } else {
    tools::sg::style_from_res(a_styles,a_style,a_plotter,false);
  }
}

inline tools::xml::styles::style_t* find_style(tools::xml::styles& a_styles,const std::string& a_name) {
  tools_vforit(tools::xml::styles::named_style_t,a_styles.named_styles(),it){
    if((*it).first==a_name) return &((*it).second);
  }
  return 0;
}

inline void SetPlotterStyles(tools::sg::plots& a_plots,
                             const std::vector<G4String>& a_plotter_styles,
                             const std::vector<G4Plotter::RegionStyle>& a_region_styles) {

  G4PlotterManager::Styles& _styles = G4PlotterManager::GetInstance().GetStyles();
 
  tools::xml::styles _tools_styles(G4cout);
  _tools_styles.add_colormap("default",tools::sg::style_default_colormap());
  _tools_styles.add_colormap("ROOT",tools::sg::style_ROOT_colormap());
  
 {tools_vforcit(G4PlotterManager::NamedStyle,_styles,it) {
    tools::xml::styles::style_t _tools_style;
    tools_vforcit(G4PlotterManager::StyleItem,(*it).second,its) {
      const G4String& param = (*its).first;
      if(param.find('.')==std::string::npos) {
        const G4String& value = (*its).second;
        _tools_style.push_back(tools::xml::styles::style_item_t(param,value));
      }
    }
    _tools_styles.add_style((*it).first,_tools_style);
  }}

  // sub styles:
 {tools_vforcit(G4PlotterManager::NamedStyle,_styles,it) {
    tools_vforcit(G4PlotterManager::StyleItem,(*it).second,its) {
      const G4String& param = (*its).first;
      std::string::size_type pos = param.rfind('.');
      if(pos!=std::string::npos) {
	std::string sub_style = (*it).first+"."+param.substr(0,pos);
        G4String parameter = param.substr(pos+1,param.size()-pos);
        const G4String& value = (*its).second;
        tools::xml::styles::style_t* _tools_style = find_style(_tools_styles,sub_style);
	if(_tools_style) {
          _tools_style->push_back(tools::xml::styles::style_item_t(parameter,value));
	} else {
          tools::xml::styles::style_t _tools_style_2;
          _tools_style_2.push_back(tools::xml::styles::style_item_t(parameter,value));
          _tools_styles.add_style(sub_style,_tools_style_2);
	}
      }
    }
  }}

 {unsigned int number = a_plots.number();
  for(unsigned int index=0;index<number;index++) {
    tools::sg::plotter* _plotter = a_plots.find_plotter(index);
    if(_plotter) {
      tools_vforcit(G4String,a_plotter_styles,it) {
        SetRegionStyles(_tools_styles,a_plots,*_plotter,*it);
      }
    }
  }}
 {tools_vforcit(G4Plotter::RegionStyle,a_region_styles,it) {
    tools::sg::plotter* _plotter = a_plots.find_plotter((*it).first);
    if(_plotter) {
      SetRegionStyles(_tools_styles,a_plots,*_plotter,(*it).second);
    }
  }}
}

inline void SetPlotterParameters(tools::sg::cmaps_t& a_cmaps,tools::sg::plots& a_plots,
                                 const std::vector<G4Plotter::RegionParameter>& a_region_parameters) {
  // parameter/field examples :
  //   title_automated
  //   title
  //   bins_style.0.color
  //   x_axis.divisions
  //   x_axis.line_style.color
  //   background_style.back_color
  tools_vforcit(G4Plotter::RegionParameter,a_region_parameters,it) {
    tools::sg::plotter* _plotter = a_plots.find_plotter((*it).first);
    if(_plotter) {
      const G4String& parameter = (*it).second.first;
      const G4String& value = (*it).second.second;
      tools::sg::field* fd = _plotter->find_field_by_name(parameter);
      if(!fd) fd = _plotter->find_field_by_name(_plotter->s_cls()+"."+parameter);
      if(fd) {if(fd->s2value(value)) continue;}
      // look for sf_enum for which value is given with a string, or
      // for sf<bool> for which value given with true/false, or
      // for a style, for example: bins_style.0.color:
      if(!_plotter->set_from_string(G4cout,a_cmaps,parameter,value)) {
        G4cout << "G4ToolsSGSceneHandler::SetPlotterParameters: plotter.set_from_string() failed for field "
               << tools::sout(parameter) << ", and value " << tools::sout(value) << "."
               << std::endl;
      }
    }
  }
}

#include "G4UImanager.hh"

void G4ToolsSGSceneHandler::SetPlotterHistograms(tools::sg::plots& a_plots) {
  a_plots.clear();
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI==NULL) return;
 {tools_vforcit(Region_h1,fRegionH1s,it) {
    tools::sg::plotter* _plotter = a_plots.find_plotter((*it).first);
    if(_plotter) {
      int hid = (*it).second;
      std::ostringstream os;
      os << hid;
      std::string cmd("/analysis/h1/get ");
      cmd += std::string(os.str());
      auto keepControlVerbose = UI->GetVerboseLevel();
      UI->SetVerboseLevel(0);
      G4int status = UI->ApplyCommand(cmd.c_str());
      UI->SetVerboseLevel(keepControlVerbose);
      if(status==G4UIcommandStatus::fCommandSucceeded) {
        G4String hexString = UI->GetCurrentValues("/analysis/h1/get");
        if(hexString.size()) {
          void* ptr;
          std::istringstream is(hexString);
          is >> ptr;
          tools::histo::h1d* _h = (tools::histo::h1d*)ptr;
          tools::sg::plottable* p = new tools::sg::h1d2plot_cp(*_h);
          _plotter->add_plottable(p); //give ownership of p to sg::plotter.
        }
      } else {
        G4cerr <<
        "G4ToolsSGSceneHandler::SetPlotterHistograms: cannot get histogram - maybe doesn't exist?"
        "\n  Maybe this app does not do analysis at all?"
        << G4endl;
      }
    }
  }}
 {tools_vforcit(Region_h2,fRegionH2s,it) {
    tools::sg::plotter* _plotter = a_plots.find_plotter((*it).first);
    if(_plotter) {
      int hid = (*it).second;
      std::ostringstream os;
      os << hid;
      std::string cmd("/analysis/h2/get ");
      cmd += std::string(os.str());
      auto keepControlVerbose = UI->GetVerboseLevel();
      UI->SetVerboseLevel(0);
      G4int status = UI->ApplyCommand(cmd.c_str());
      UI->SetVerboseLevel(keepControlVerbose);
      if(status==G4UIcommandStatus::fCommandSucceeded) {
        G4String hexString = UI->GetCurrentValues("/analysis/h2/get");
        if(hexString.size()) {
          void* ptr;
          std::istringstream is(hexString);
          is >> ptr;
          tools::histo::h2d* _h = (tools::histo::h2d*)ptr;
          tools::sg::plottable* p = new tools::sg::h2d2plot_cp(*_h);
          _plotter->add_plottable(p); //give ownership of p to sg::plotter.
        }
      } else {
        G4cerr <<
        "G4ToolsSGSceneHandler::SetPlotterHistograms: cannot get histogram - maybe doesn't exist?"
        "\n  Maybe this app does not do analysis at all?"
        << G4endl;
      }
    }
  }}
}

class plots_cbk : public tools::sg::ecbk {
  TOOLS_CBK(plots_cbk,plots_cbk,tools::sg::ecbk)
public:
  virtual tools::sg::return_action action() {
    if(const tools::sg::size_event* sz_evt = tools::sg::event_cast<tools::sg::event,tools::sg::size_event>(*m_event)){
      m_plots.adjust_size(sz_evt->width(),sz_evt->height());
      m_event_action->set_done(true);
      return tools::sg::return_to_render;
    }
    return tools::sg::return_none;
  }
public:
  plots_cbk(tools::sg::plots& a_plots)
  :parent()
  ,m_plots(a_plots)
  {}
  virtual ~plots_cbk(){}
public:
  plots_cbk(const plots_cbk& a_from)
  :parent(a_from)
  ,m_plots(a_from.m_plots)
  {}
  plots_cbk& operator=(const plots_cbk& a_from){
    parent::operator=(a_from);
    return *this;
  }
protected:
  tools::sg::plots& m_plots;
};

void G4ToolsSGSceneHandler::TouchPlotters(tools::sg::node& a_sg) {
  tools::sg::search_action sa(G4cout);
  const tools::sg::search_action::paths_t& paths = tools::sg::find_paths<tools::sg::plots>(sa,a_sg);
  tools_vforcit(tools::sg::path_t,paths,it) {
    tools::sg::plots* _plots = tools::sg::tail<tools::sg::plots>(*it);    
    if(_plots) {
      SetPlotterHistograms(*_plots);
    }
  }
}

void G4ToolsSGSceneHandler::AddPrimitive(const G4Plotter& a_plotter)
{
//G4cout << "debug : G4ToolsSGSceneHandler::AddPrimitive : 004" << std::endl;
  if(!fpViewer) return;
  
  auto currentNode = GetOrCreateNode();
  if (!currentNode) return;  // Node not available

  currentNode->add(new tools::sg::light_off());
  
  tools::sg::plots* _plots = new tools::sg::plots(*fFreetypeNode);
  currentNode->add(_plots);
  
  _plots->view_border = false;
  _plots->set_regions(a_plotter.GetColumns(),a_plotter.GetRows());

 {tools::sg::event_dispatcher* dpt = new tools::sg::event_dispatcher;
  dpt->add_callback(new plots_cbk(*_plots));
  currentNode->add(dpt);}

  SetPlotterStyles(*_plots,a_plotter.GetStyles(),a_plotter.GetRegionStyles());

  tools::sg::cmaps_t _cmaps;
  _cmaps["default"] = tools::sg::style_default_colormap();
  _cmaps["ROOT"] = tools::sg::style_ROOT_colormap();
  
  SetPlotterParameters(_cmaps,*_plots,a_plotter.GetRegionParameters());

  fRegionH1s = a_plotter.GetRegionH1s();
  fRegionH2s = a_plotter.GetRegionH2s();

  SetPlotterHistograms(*_plots);
}

void G4ToolsSGSceneHandler::Messenger::SetNewValue(G4UIcommand* a_cmd,G4String) {
  G4VSceneHandler* pSceneHandler = GetVisManager()->GetCurrentSceneHandler();
  if (!pSceneHandler) {
    G4cout << "G4ToolsSGSceneHandler::Messenger::SetNewValue: no current sceneHandler.  Please create one." << G4endl;
    return;
  }
  auto* tsg_scene_handler = dynamic_cast<G4ToolsSGSceneHandler*>(pSceneHandler);
  if(!tsg_scene_handler) {
    G4cout << "G4ToolsSGSceneHandler::Messenger::SetNewValue: current sceneHandler not a G4ToolsSGSceneHandler." << G4endl;
    return;
  }
  if(a_cmd==print_plotter_params) {
    tools::sg::dummy_freetype _ttf;
    tools::sg::plotter _plotter(_ttf);
    _plotter.print_available_customization(G4cout);
  }
}
