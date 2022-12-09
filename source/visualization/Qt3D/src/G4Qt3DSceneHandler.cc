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
// John Allison  17th June 2019

#include "G4Qt3DSceneHandler.hh"

#include "G4PhysicalVolumeModel.hh"
#include "G4LogicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4TransportationManager.hh"
#include "G4Box.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Polyhedron.hh"
#include "G4Scene.hh"
#include "G4Threading.hh"
#include "G4Mesh.hh"
#include "G4PseudoScene.hh"
#include "G4VisManager.hh"

#include "G4Qt3DViewer.hh"
#include "G4Qt3DUtils.hh"
#include "G4Qt3DQEntity.hh"

#include <Qt3DCore>
#include <Qt3DExtras>
#include <Qt3DRender>

#include <utility>

#define G4warn G4cout

// Qt3D seems to offer a choice of type - float or double. It would be nice
// to use double since it offers the prospect of higher precision, hopefully
// avoiding some issues that we see at high zoom. But it currently gives the
// following warning: "findBoundingVolumeComputeData: Position attribute not
// suited for bounding volume computation", so for now we use float.
#define PRECISION float
#define BASETYPE Qt3DRender::QAttribute::Float

G4int G4Qt3DSceneHandler::fSceneIdCount = 0;

G4Qt3DSceneHandler::G4Qt3DSceneHandler
 (G4VGraphicsSystem& system, const G4String& name)
: G4VSceneHandler(system, fSceneIdCount++, name)
{
#ifdef G4QT3DDEBUG
  G4cout << "G4Qt3DSceneHandler::G4Qt3DSceneHandler called" << G4endl;
#endif
  fpQt3DScene = new Qt3DCore::QEntity;
  fpQt3DScene->setObjectName("G4Qt3DSceneRoot");
  EstablishG4Qt3DQEntities();
}

G4Qt3DSceneHandler::~G4Qt3DSceneHandler()
{
  // Doesn't like this - it gives BAD_ACCESS in delete_entity_recursively.
  // Curiously the delete traceback shows three calls to this recursively:
  /*#1  0x0000000100411906 in (anonymous namespace)::delete_entity_recursively(Qt3DCore::QNode*) at /Users/johna/Geant4/geant4-dev/source/visualization/Qt3D/src/G4Qt3DSceneHandler.cc:60
  #2  0x0000000100411840 in G4Qt3DSceneHandler::~G4Qt3DSceneHandler() at /Users/johna/Geant4/geant4-dev/source/visualization/Qt3D/src/G4Qt3DSceneHandler.cc:169
  #3  0x0000000100411fc5 in G4Qt3DSceneHandler::~G4Qt3DSceneHandler() at /Users/johna/Geant4/geant4-dev/source/visualization/Qt3D/src/G4Qt3DSceneHandler.cc:168
  #4  0x0000000100411fe9 in G4Qt3DSceneHandler::~G4Qt3DSceneHandler() at /Users/johna/Geant4/geant4-dev/source/visualization/Qt3D/src/G4Qt3DSceneHandler.cc:168
  #5  0x0000000101032510 in G4VisManager::~G4VisManager() at /Users/johna/Geant4/geant4-dev/source/visualization/management/src/G4VisManager.cc:214
  #6  0x0000000100013885 in G4VisExecutive::~G4VisExecutive() at /Users/johna/Geant4/geant4-dev/source/visualization/management/include/G4VisExecutive.hh:119
  #7  0x00000001000119a5 in G4VisExecutive::~G4VisExecutive() at /Users/johna/Geant4/geant4-dev/source/visualization/management/include/G4VisExecutive.hh:119
  #8  0x00000001000119c9 in G4VisExecutive::~G4VisExecutive() at /Users/johna/Geant4/geant4-dev/source/visualization/management/include/G4VisExecutive.hh:119
  #9  0x00000001000117dd in main at /Users/johna/Geant4/geant4-dev/examples/basic/B1/exampleB1.cc:108
  */
  //if (fpQt3DScene) delete_entity_recursively(fpQt3DScene);
}

void G4Qt3DSceneHandler::EstablishG4Qt3DQEntities()
{
  fpTransientObjects  = new G4Qt3DQEntity(fpQt3DScene);  // Hangs from root
  fpTransientObjects  ->setObjectName("G4Qt3DTORoot");
  fpPersistentObjects = new G4Qt3DQEntity(fpQt3DScene);  // Hangs from root
  fpPersistentObjects ->setObjectName("G4Qt3DPORoot");

  // Physical volume objects for each world hang from POs
  G4TransportationManager* transportationManager
  = G4TransportationManager::GetTransportationManager ();
  std::size_t nWorlds = transportationManager->GetNoWorlds();
  std::vector<G4VPhysicalVolume*>::iterator iterWorld
  = transportationManager->GetWorldsIterator();
  fpPhysicalVolumeObjects.resize(nWorlds);
  for (std::size_t i = 0; i < nWorlds; ++i, ++iterWorld) {
    G4VPhysicalVolume* wrld = (*iterWorld);
    auto entity = new G4Qt3DQEntity(fpPersistentObjects);
    entity->setObjectName("G4Qt3DPORoot_"+QString(wrld->GetName()));
    entity->SetPVNodeID(G4PhysicalVolumeModel::G4PhysicalVolumeNodeID(wrld));
    fpPhysicalVolumeObjects[i] = entity;
  }
}

G4Qt3DQEntity* G4Qt3DSceneHandler::CreateNewNode()
{
  // Create a G4Qt3DQEntity node suitable for next solid or primitive

  G4Qt3DQEntity* newNode = nullptr;

  if (fReadyForTransients) {  // All transients hang from this node
    newNode = new G4Qt3DQEntity(fpTransientObjects);
    G4String name = fpModel? fpModel->GetGlobalTag(): "User";
    newNode->setObjectName(name.c_str());
    return newNode;
  }

  G4PhysicalVolumeModel* pPVModel =
    dynamic_cast<G4PhysicalVolumeModel*>(fpModel);

  if (!pPVModel) {  // Persistent objects (e.g., axes)
    newNode = new G4Qt3DQEntity(fpPersistentObjects);
    newNode->setObjectName(fpModel->GetGlobalTag().c_str());
    return newNode;
  }

  // So this is a G4PhysicalVolumeModel
  
  typedef G4PhysicalVolumeModel::G4PhysicalVolumeNodeID PVNodeID;
  typedef std::vector<PVNodeID> PVPath;
//  const PVPath& drawnPVPath = pPVModel->GetDrawnPVPath();
  const PVPath& fullPVPath  = pPVModel->GetFullPVPath();
  //G4int currentDepth = pPVModel->GetCurrentDepth();
  //G4VPhysicalVolume* pCurrentPV = pPVModel->GetCurrentPV();
  //G4LogicalVolume* pCurrentLV = pPVModel->GetCurrentLV();
  //G4Material* pCurrentMaterial = pPVModel->GetCurrentMaterial();
  // Note: pCurrentMaterial may be zero (parallel world).

#ifdef G4QTDEBUG
  G4cout << "A: " << fullPVPath << G4endl;  // DEBUG
#endif

  // Find appropriate root
  const std::size_t nWorlds = fpPhysicalVolumeObjects.size();
  std::size_t iWorld = 0;
  for (; iWorld < nWorlds; ++iWorld) {
    if (fullPVPath[0].GetPhysicalVolume() ==
        fpPhysicalVolumeObjects[iWorld]->GetPVNodeID().GetPhysicalVolume()) break;
  }
  if (iWorld == nWorlds) {
    G4Exception("G4Qt3DSceneHandler::CreateNewNode", "qt3D-0000", FatalException,
                "World mis-match - not possible(!?)");
  }

  // (Re-)establish pv path of root entity
  G4Qt3DQEntity* wrld = fpPhysicalVolumeObjects[iWorld];
  wrld->SetPVNodeID(fullPVPath[0]);

  // Create nodes as required
  G4Qt3DQEntity* node = wrld;
  newNode = node;
  const std::size_t depth = fullPVPath.size();
  std::size_t iDepth = 1;
  while (iDepth < depth) {
    const auto& children = node->children();
    const G4int nChildren = children.size();  // int size() (Qt covention?)
    G4int iChild = 0;
    G4Qt3DQEntity* child = nullptr;
    for (; iChild < nChildren; ++iChild) {
      child = static_cast<G4Qt3DQEntity*>(children[iChild]);
      if (child->GetPVNodeID() == fullPVPath[iDepth]) break;
    }
    if (iChild != nChildren) {  // Existing node found
      node = child;  // Must be the ancestor of new node (subsequent iteration)
    } else {
      // Add a new node as child of node
      newNode = new G4Qt3DQEntity(node);
      newNode->SetPVNodeID(fullPVPath[iDepth]);
      std::ostringstream oss;
      oss << newNode->GetPVNodeID().GetPhysicalVolume()->GetName()
      << ':' << newNode->GetPVNodeID().GetCopyNo();
      newNode->setObjectName(oss.str().c_str());
      node = newNode;
    }
    ++iDepth;
  }

  return node;
}

void G4Qt3DSceneHandler::PreAddSolid
 (const G4Transform3D& objectTransformation,
  const G4VisAttributes& visAttribs)
{  
  G4VSceneHandler::PreAddSolid(objectTransformation, visAttribs);
}

void G4Qt3DSceneHandler::PostAddSolid()
{
  G4VSceneHandler::PostAddSolid();
}

void G4Qt3DSceneHandler::BeginPrimitives2D(const G4Transform3D& objectTransformation)
{
// The x,y coordinates of the primitives passed to AddPrimitive are
// intrepreted as screen coordinates, -1 < x,y < 1.  The
// z-coordinate is ignored.
// IMPORTANT: invoke this from your polymorphic versions, e.g.:
// void MyXXXSceneHandler::BeginPrimitives2D
// (const G4Transform3D& objectTransformation) {
  static G4bool first = true;
  if (first) {
    first = false;
    G4Exception("G4Qt3DSceneHandler::BeginPrimitives2D", "qt3D-0001",
                JustWarning,
                "2D drawing not yet implemented");
  }
   G4VSceneHandler::BeginPrimitives2D (objectTransformation);
//   ...
}

void G4Qt3DSceneHandler::EndPrimitives2D ()
{
// IMPORTANT: invoke this from your polymorphic versions, e.g.:
// void MyXXXSceneHandler::EndPrimitives2D () {
//   ...
  G4VSceneHandler::EndPrimitives2D ();
}

void G4Qt3DSceneHandler::BeginPrimitives
 (const G4Transform3D& objectTransformation)
{
  G4VSceneHandler::BeginPrimitives(objectTransformation);
}

void G4Qt3DSceneHandler::EndPrimitives ()
{
  G4VSceneHandler::EndPrimitives ();
}

void G4Qt3DSceneHandler::AddPrimitive(const G4Polyline& polyline)
{
#ifdef G4QT3DDEBUG
  G4cout <<
  "G4Qt3DSceneHandler::AddPrimitive(const G4Polyline& polyline) called.\n"
  << polyline
  << G4endl;
#endif
  
  if (polyline.size() == 0) return;

  auto currentNode = CreateNewNode();
  if (!currentNode) {
    static G4bool first = true;
    if (first) {
      first = false;
      G4Exception("G4Qt3DSceneHandler::AddPrimitive(const G4Polyline&)",
		  "qt3d-0003", JustWarning,
		  "No available node!");
    }
    return;
  }

  fpVisAttribs = fpViewer->GetApplicableVisAttributes(polyline.GetVisAttributes());

  auto transform = G4Qt3DUtils::CreateQTransformFrom(fObjectTransformation);
  transform->setObjectName("transform");

  auto polylineEntity = new Qt3DCore::QEntity(currentNode);
  polylineEntity->addComponent(transform);

  const auto vertexByteSize  = 3*sizeof(PRECISION);

  const std::size_t nLines = polyline.size() - 1;
  QByteArray polylineByteArray;
  const auto polylineBufferByteSize = 2*nLines*vertexByteSize;
  polylineByteArray.resize((G4int)polylineBufferByteSize);
  auto polylineBufferArray = reinterpret_cast<PRECISION*>(polylineByteArray.data());
  G4int iLine = 0;
  for (std::size_t i = 0; i < nLines; ++i) {
    polylineBufferArray[iLine++] = polyline[i].x();
    polylineBufferArray[iLine++] = polyline[i].y();
    polylineBufferArray[iLine++] = polyline[i].z();
    polylineBufferArray[iLine++] = polyline[i+1].x();
    polylineBufferArray[iLine++] = polyline[i+1].y();
    polylineBufferArray[iLine++] = polyline[i+1].z();
  }
  auto polylineGeometry = new Qt3DRender::QGeometry();
  polylineGeometry->setObjectName("polylineGeometry");
  auto polylineBuffer = new Qt3DRender::QBuffer(polylineGeometry);
  polylineBuffer->setObjectName("Polyline buffer");
  polylineBuffer->setData(polylineByteArray);

  auto polylineAtt = new Qt3DRender::QAttribute;
  polylineAtt->setObjectName("Position attribute");
  polylineAtt->setName(Qt3DRender::QAttribute::defaultPositionAttributeName());
  polylineAtt->setBuffer(polylineBuffer);
  polylineAtt->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
  polylineAtt->setVertexBaseType(BASETYPE);
  polylineAtt->setVertexSize(3);
  polylineAtt->setCount((G4int)nLines);
  polylineAtt->setByteOffset(0);
  polylineAtt->setByteStride(vertexByteSize);

  const auto& colour = fpVisAttribs->GetColour();

  polylineGeometry->addAttribute(polylineAtt);

  auto material = new Qt3DExtras::QDiffuseSpecularMaterial();
  material->setObjectName("materialForPolyline");
  material->setAmbient(G4Qt3DUtils::ConvertToQColor(colour));
  material->setShininess(0.);
  material->setSpecular(0.);
  polylineEntity->addComponent(material);

  auto renderer = new Qt3DRender::QGeometryRenderer;
  renderer->setObjectName("polylineWireframeRenderer");
  renderer->setGeometry(polylineGeometry);
  renderer->setVertexCount(2*(G4int)nLines);
  renderer->setPrimitiveType(Qt3DRender::QGeometryRenderer::Lines);
  polylineEntity->addComponent(renderer);
}

void G4Qt3DSceneHandler::AddPrimitive (const G4Polymarker& polymarker)
{
  if (polymarker.size() == 0) return;

  auto currentNode = CreateNewNode();
  if (!currentNode) {
    static G4bool first = true;
    if (first) {
      first = false;
      G4Exception("G4Qt3DSceneHandler::AddPrimitive(const G4Polymarker&)",
		  "qt3d-0003", JustWarning,
		  "No available node!");
    }
    return;
  }

  fpVisAttribs = fpViewer->GetApplicableVisAttributes(polymarker.GetVisAttributes());

  MarkerSizeType markerSizeType;
  G4double markerSize = GetMarkerSize(polymarker, markerSizeType);

  switch (polymarker.GetMarkerType()) {
    default:
    case G4Polymarker::dots:
    {
      const std::size_t nDots = polymarker.size();

      auto transform = G4Qt3DUtils::CreateQTransformFrom(fObjectTransformation);
      transform->setObjectName("transform");

      auto polymarkerEntity = new Qt3DCore::QEntity(currentNode);
      polymarkerEntity->addComponent(transform);

      const auto vertexByteSize  = 3*sizeof(PRECISION);

      QByteArray polymarkerByteArray;
      const auto polymarkerBufferByteSize = nDots*vertexByteSize;
      polymarkerByteArray.resize((G4int)polymarkerBufferByteSize);
      auto polymarkerBufferArray = reinterpret_cast<PRECISION*>(polymarkerByteArray.data());
      G4int iMarker = 0;
      for (std::size_t i = 0; i < polymarker.size(); ++i) {
        polymarkerBufferArray[iMarker++] = polymarker[i].x();
        polymarkerBufferArray[iMarker++] = polymarker[i].y();
        polymarkerBufferArray[iMarker++] = polymarker[i].z();
      }
      auto polymarkerGeometry = new Qt3DRender::QGeometry();
      polymarkerGeometry->setObjectName("polymarkerGeometry");
      auto polymarkerBuffer = new Qt3DRender::QBuffer(polymarkerGeometry);
      polymarkerBuffer->setObjectName("Polymarker buffer");
      polymarkerBuffer->setData(polymarkerByteArray);

      auto polymarkerAtt = new Qt3DRender::QAttribute;
      polymarkerAtt->setObjectName("Position attribute");
      polymarkerAtt->setName(Qt3DRender::QAttribute::defaultPositionAttributeName());
      polymarkerAtt->setBuffer(polymarkerBuffer);
      polymarkerAtt->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
      polymarkerAtt->setVertexBaseType(BASETYPE);
      polymarkerAtt->setVertexSize(3);
      polymarkerAtt->setCount((G4int)nDots);
      polymarkerAtt->setByteOffset(0);
      polymarkerAtt->setByteStride(vertexByteSize);

      const auto& colour = fpVisAttribs->GetColour();

      polymarkerGeometry->addAttribute(polymarkerAtt);

      auto material = new Qt3DExtras::QDiffuseSpecularMaterial();
      material->setObjectName("materialForPolymarker");
      material->setAmbient(G4Qt3DUtils::ConvertToQColor(colour));
      material->setShininess(0.);
      material->setSpecular(0.);
      polymarkerEntity->addComponent(material);

      auto renderer = new Qt3DRender::QGeometryRenderer;
      renderer->setObjectName("polymarkerWireframeRenderer");
      renderer->setGeometry(polymarkerGeometry);
      renderer->setVertexCount((G4int)nDots);
      renderer->setPrimitiveType(Qt3DRender::QGeometryRenderer::Points);
      polymarkerEntity->addComponent(renderer);
    }
      break;
    case G4Polymarker::circles:
    {
      G4Circle circle (polymarker);  // Default circle

      const auto& colour = fpVisAttribs->GetColour();
      auto material = new Qt3DExtras::QDiffuseSpecularMaterial();
      material->setObjectName("materialForCircle");
      material->setAmbient(G4Qt3DUtils::ConvertToQColor(colour));
      if (colour.GetAlpha() < 1.) material->setAlphaBlendingEnabled(true);

      auto sphereMesh = new Qt3DExtras::QSphereMesh;
      sphereMesh->setObjectName("sphereMesh");
      G4double radius = markerSize/2.;
      if (markerSizeType == G4VSceneHandler::screen ) {
        // Not figured out how to do screen-size, so use scene extent
        const G4double scale = 200.;  // Roughly pixels per scene
        radius *= fpScene->GetExtent().GetExtentRadius()/scale;
      }
      sphereMesh->setRadius(radius);
//      sphereMesh->setInstanceCount(polymarker.size());  // Not undertood instancing yet

//      auto currentEntity = new Qt3DCore::QEntity(currentNode);  // Not undertood instancing yet
      for (std::size_t iPoint = 0; iPoint < polymarker.size(); ++iPoint) {
        auto position = fObjectTransformation*G4Translate3D(polymarker[iPoint]);
        auto transform = G4Qt3DUtils::CreateQTransformFrom(position);
	auto currentEntity = new Qt3DCore::QEntity(currentNode);  // Not undertood instancing yet
        currentEntity->addComponent(material);
        currentEntity->addComponent(transform);
        currentEntity->addComponent(sphereMesh);
      }
    }
      break;
    case G4Polymarker::squares:
    {
      G4Square square (polymarker);  // Default square

      const auto& colour = fpVisAttribs->GetColour();
      auto material = new Qt3DExtras::QDiffuseSpecularMaterial();
      material->setObjectName("materialForSquare");
      material->setAmbient(G4Qt3DUtils::ConvertToQColor(colour));
      if (colour.GetAlpha() < 1.) material->setAlphaBlendingEnabled(true);

      auto boxMesh = new Qt3DExtras::QCuboidMesh();
      boxMesh->setObjectName("boxMesh");
      G4double side = markerSize;
      if (markerSizeType == G4VSceneHandler::screen ) {
        // Not figured out how to do screen-size, so use scene extent
        const G4double scale = 200.;  // Roughly pixles per scene
        side *= fpScene->GetExtent().GetExtentRadius()/scale;
      }
      boxMesh->setXExtent(side);
      boxMesh->setYExtent(side);
      boxMesh->setZExtent(side);

      for (std::size_t iPoint = 0; iPoint < polymarker.size(); ++iPoint) {
        auto position = fObjectTransformation*G4Translate3D(polymarker[iPoint]);
        auto transform = G4Qt3DUtils::CreateQTransformFrom(position);
        auto currentEntity = new Qt3DCore::QEntity(currentNode);
        currentEntity->addComponent(material);
        currentEntity->addComponent(transform);
        currentEntity->addComponent(boxMesh);
      }
    }
      break;
  }
}

#ifdef G4QT3DDEBUG
void G4Qt3DSceneHandler::AddPrimitive(const G4Text& text) {
  G4cout <<
  "G4Qt3DSceneHandler::AddPrimitive(const G4Text& text) called.\n"
  << text
  << G4endl;
#else
void G4Qt3DSceneHandler::AddPrimitive(const G4Text& /*text*/) {
#endif

  static G4bool first = true;
  if (first) {
    first = false;
    G4Exception("G4Qt3DSceneHandler::AddPrimitive(const G4Text&)",
                "qt3D-0002", JustWarning,
                "Text drawing doesn't work yet");
  }  // OK. Not working, but let it execute, which it does without error.

  /* But it crashes after /vis/viewer/rebuild!!!
  auto currentNode = CreateNewNode();
   if (!currentNode) {
   static G4bool first = true;
   if (first) {
   first = false;
   G4Exception("G4Qt3DSceneHandler::AddPrimitive(const G4Text&)",
   "qt3d-0003", JustWarning,
   "No available node!");
   }
   return;
   }

  fpVisAttribs = fpViewer->GetApplicableVisAttributes(text.GetVisAttributes());

  auto position = fObjectTransformation*G4Translate3D(text.GetPosition());
  auto transform = G4Qt3DUtils::CreateQTransformFrom(position);
//  transform->setScale(10);
  transform->setScale(0.1);

//  auto currentEntity = new Qt3DCore::QEntity(currentNode);

  // This simply does not work
  auto qtext = new Qt3DExtras::QText2DEntity();
  qtext->setParent(currentNode);
//  qtext->setParent(currentEntity);  // ??  Doesn't help
  qtext->setText(text.GetText().c_str());
//  qtext->setHeight(100);
//  qtext->setWidth(1000);
  qtext->setHeight(20);
  qtext->setWidth(100);
  qtext->setColor(Qt::green);
  qtext->setFont(QFont("Courier New", 10));
  qtext->addComponent(transform);

  // This produces text in 3D facing +z - not what we want
//  const auto& colour = GetTextColour(text);
//  auto material = new Qt3DExtras::QDiffuseSpecularMaterial();
//  material->setObjectName("materialForText");
//  material->setAmbient(G4Qt3DUtils::ConvertToQColor(colour));
//  if (colour.GetAlpha() < 1.) material->setAlphaBlendingEnabled(true);
//
//  auto textMesh = new Qt3DExtras::QExtrudedTextMesh();
//  textMesh->setText(text.GetText().c_str());
//  textMesh->setFont(QFont("Courier New", 10));
//  textMesh->setDepth(.01f);
//
//  currentNode->addComponent(material);
//  currentNode->addComponent(transform);
//  currentNode->addComponent(textMesh);
   */
}

void G4Qt3DSceneHandler::AddPrimitive(const G4Circle& circle)
{
#ifdef G4QT3DDEBUG
  G4cout <<
  "G4Qt3DSceneHandler::AddPrimitive(const G4Circle& circle) called.\n"
  << circle
  << G4endl;
#endif

#ifdef G4QT3DDEBUG
  MarkerSizeType sizeType;
  G4double size = GetMarkerSize (circle, sizeType);
  switch (sizeType) {
    default:
    case screen:
      // Draw in screen coordinates.
      G4cout << "screen";
      break;
    case world:
      // Draw in world coordinates.
      G4cout << "world";
      break;
  }
  G4cout << " size: " << size << G4endl;
#endif

  auto currentNode = CreateNewNode();
  if (!currentNode) {
    static G4bool first = true;
    if (first) {
      first = false;
      G4Exception("G4Qt3DSceneHandler::AddPrimitive(const G4Circle&)",
		  "qt3d-0003", JustWarning,
		  "No available node!");
    }
    return;
  }

  fpVisAttribs = fpViewer->GetApplicableVisAttributes(circle.GetVisAttributes());

  auto position = fObjectTransformation*G4Translate3D(circle.GetPosition());
  auto transform = G4Qt3DUtils::CreateQTransformFrom(position);

  const auto& colour = fpVisAttribs->GetColour();
  auto material = new Qt3DExtras::QDiffuseSpecularMaterial();
  material->setObjectName("materialForCircle");
  material->setAmbient(G4Qt3DUtils::ConvertToQColor(colour));
  if (colour.GetAlpha() < 1.) material->setAlphaBlendingEnabled(true);

  auto sphereMesh = new Qt3DExtras::QSphereMesh;
  sphereMesh->setObjectName("sphereMesh");
  G4double radius;
  if (circle.GetSizeType() == G4VMarker::world ) {
    radius =circle.GetWorldRadius();
  } else {  // Screen-size or none
    // Not figured out how to do screen-size, so use scene extent
    const G4double scale = 200.;  // Roughly pixles per scene
    radius = circle.GetScreenRadius()*fpScene->GetExtent().GetExtentRadius()/scale;
  }
  sphereMesh->setRadius(radius);

  auto currentEntity = new Qt3DCore::QEntity(currentNode);
  currentEntity->addComponent(material);
  currentEntity->addComponent(transform);
  currentEntity->addComponent(sphereMesh);
}

void G4Qt3DSceneHandler::AddPrimitive(const G4Square& square)
{
#ifdef G4QT3DDEBUG
  G4cout <<
  "G4Qt3DSceneHandler::AddPrimitive(const G4Square& square) called.\n"
  << square
  << G4endl;
#endif

#ifdef G4QT3DDEBUG
  MarkerSizeType sizeType;
  G4double size = GetMarkerSize (square, sizeType);
  switch (sizeType) {
    default:
    case screen:
      // Draw in screen coordinates.
      G4cout << "screen";
      break;
    case world:
      // Draw in world coordinates.
      G4cout << "world";
      break;
  }
  G4cout << " size: " << size << G4endl;
#endif

  auto currentNode = CreateNewNode();
  if (!currentNode) {
    static G4bool first = true;
    if (first) {
      first = false;
      G4Exception("G4Qt3DSceneHandler::AddPrimitive(const G4Square&)",
		  "qt3d-0003", JustWarning,
		  "No available node!");
    }
    return;
  }

  fpVisAttribs = fpViewer->GetApplicableVisAttributes(square.GetVisAttributes());

  auto position = fObjectTransformation*G4Translate3D(square.GetPosition());
  auto transform = G4Qt3DUtils::CreateQTransformFrom(position);

  const auto& colour = fpVisAttribs->GetColour();
  auto material = new Qt3DExtras::QDiffuseSpecularMaterial();
  material->setObjectName("materialForSquare");
  material->setAmbient(G4Qt3DUtils::ConvertToQColor(colour));
  if (colour.GetAlpha() < 1.) material->setAlphaBlendingEnabled(true);

  auto boxMesh = new Qt3DExtras::QCuboidMesh();
  boxMesh->setObjectName("boxMesh");
  G4double side;
  if (square.GetSizeType() == G4VMarker::world ) {
    side = square.GetWorldDiameter();
  } else {  // Screen-size or none
    // Not figured out how to do screen-size, so use scene extent
    const G4double scale = 200.;  // Roughly pixles per scene
    side = square.GetScreenDiameter()*fpScene->GetExtent().GetExtentRadius()/scale;
  }
  boxMesh->setXExtent(side);
  boxMesh->setYExtent(side);
  boxMesh->setZExtent(side);

  auto currentEntity = new Qt3DCore::QEntity(currentNode);
  currentEntity->addComponent(material);
  currentEntity->addComponent(transform);
  currentEntity->addComponent(boxMesh);
}

void G4Qt3DSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron)
{
  auto currentNode = CreateNewNode();
  if (!currentNode) {
    static G4bool first = true;
    if (first) {
      first = false;
      G4Exception("G4Qt3DSceneHandler::AddPrimitive(const G4Polyhedron&)",
		  "qt3d-0003", JustWarning,
		  "No available node!");
    }
    return;
  }

  if (polyhedron.GetNoFacets() == 0) return;

  fpVisAttribs = fpViewer->GetApplicableVisAttributes(polyhedron.GetVisAttributes());

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
    // For a large polyhedron, eliminating lines like this is prohibitively
    // expensive. Comment out for now, and maybe unwind altogether in future.
    // Allow the graphics-reps utilities to optimise things like this.
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
    notLastFace = polyhedron.GetNextFacet(nEdges, vertex, edgeFlag, normal);
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
      G4warn
      << "ERROR: polyhedron face with unexpected number of edges (" << nEdges << ')'
      << "\n  Tag: " << fpModel->GetCurrentTag()
      << G4endl;
      return;
    }
  } while (notLastFace);
  const auto nVerts = vertices.size();
  const auto nLines = lines.size();

  // Now put stuff into Qt objects

  auto transform = G4Qt3DUtils::CreateQTransformFrom(fObjectTransformation);
  transform->setObjectName("transform");

  Qt3DCore::QEntity* wireframeEntity = nullptr;
  Qt3DCore::QEntity* surfaceEntity   = nullptr;
  static G4int errorCount = 0;
  G4ViewParameters::DrawingStyle drawing_style = GetDrawingStyle (fpVisAttribs);
  switch (drawing_style) {
    case G4ViewParameters::wireframe:
      wireframeEntity = new Qt3DCore::QEntity(currentNode);
      wireframeEntity->addComponent(transform);
      break;
    case G4ViewParameters::hlr:
      wireframeEntity = new Qt3DCore::QEntity(currentNode);
      wireframeEntity->addComponent(transform);
      surfaceEntity = new Qt3DCore::QEntity(currentNode);
      surfaceEntity->addComponent(transform);
      break;
    case G4ViewParameters::hsr:
      surfaceEntity = new Qt3DCore::QEntity(currentNode);
      surfaceEntity->addComponent(transform);
      break;
    case G4ViewParameters::hlhsr:
      wireframeEntity = new Qt3DCore::QEntity(currentNode);
      wireframeEntity->addComponent(transform);
      surfaceEntity = new Qt3DCore::QEntity(currentNode);
      surfaceEntity->addComponent(transform);
      break;
    case G4ViewParameters::cloud:
      // Shouldn't happen in this function (it's a polyhedron!)
      if (errorCount == 0) {
        ++errorCount;
        G4warn << "WARNING: Qt3D: cloud drawing not implemented" << G4endl;
      }
      return;
      break;
  }

  const auto vertexByteSize  = 3*sizeof(PRECISION);

  Qt3DRender::QGeometry* vertexGeometry = nullptr;
  Qt3DRender::QGeometry* lineGeometry   = nullptr;

  Qt3DRender::QAttribute* positionAtt = nullptr;
  Qt3DRender::QAttribute* normalAtt   = nullptr;
  Qt3DRender::QAttribute* lineAtt     = nullptr;

  Qt3DRender::QBuffer* vertexBuffer = nullptr;
  if (drawing_style == G4ViewParameters::hlr ||
      drawing_style == G4ViewParameters::hsr ||
      drawing_style == G4ViewParameters::hlhsr) {

    // Put vertices, normals into  QByteArray
    // Accomodates both vertices and normals - hence 2*
    QByteArray vertexByteArray;
    const auto vertexBufferByteSize = 2*nVerts*vertexByteSize;
    vertexByteArray.resize((G4int)vertexBufferByteSize);
    auto vertexBufferArray = reinterpret_cast<PRECISION*>(vertexByteArray.data());
    G4int i1 = 0;
    for (std::size_t i = 0; i < nVerts; ++i) {
      vertexBufferArray[i1++] = vertices[i].x();
      vertexBufferArray[i1++] = vertices[i].y();
      vertexBufferArray[i1++] = vertices[i].z();
      vertexBufferArray[i1++] = normals[i].x();
      vertexBufferArray[i1++] = normals[i].y();
      vertexBufferArray[i1++] = normals[i].z();
    }
    // Vertex buffer (vertices and normals)
    vertexGeometry = new Qt3DRender::QGeometry();
    vertexGeometry->setObjectName("vertexGeometry");
    vertexBuffer = new Qt3DRender::QBuffer(vertexGeometry);
    vertexBuffer->setObjectName("Vertex buffer");
    vertexBuffer->setData(vertexByteArray);

    // Position attribute
    positionAtt = new Qt3DRender::QAttribute;
    positionAtt->setObjectName("Position attribute");
    positionAtt->setName(Qt3DRender::QAttribute::defaultPositionAttributeName());
    positionAtt->setBuffer(vertexBuffer);
    positionAtt->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
    positionAtt->setVertexBaseType(BASETYPE);
    positionAtt->setVertexSize(3);
    positionAtt->setCount((G4int)nVerts);
    positionAtt->setByteOffset(0);
    positionAtt->setByteStride(2*vertexByteSize);

    // Normal attribute
    normalAtt = new Qt3DRender::QAttribute;
    normalAtt->setObjectName("Normal attribute");
    normalAtt->setName(Qt3DRender::QAttribute::defaultNormalAttributeName());
    normalAtt->setBuffer(vertexBuffer);
    normalAtt->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
    normalAtt->setVertexBaseType(BASETYPE);
    normalAtt->setVertexSize(3);
    normalAtt->setCount((G4int)nVerts);
    normalAtt->setByteOffset(vertexByteSize);
    normalAtt->setByteStride(2*vertexByteSize);
  }

  Qt3DRender::QBuffer* lineBuffer = nullptr;
  if (drawing_style == G4ViewParameters::wireframe ||
      drawing_style == G4ViewParameters::hlr ||
      drawing_style == G4ViewParameters::hlhsr) {

    // Put lines into a QByteArray
    QByteArray lineByteArray;
    const auto lineBufferByteSize = 2*nLines*vertexByteSize;
    lineByteArray.resize((G4int)lineBufferByteSize);
    auto lineBufferArray = reinterpret_cast<PRECISION*>(lineByteArray.data());
    G4int i2 = 0;
    for (const auto& line: lines) {
      lineBufferArray[i2++] = line.first.x();
      lineBufferArray[i2++] = line.first.y();
      lineBufferArray[i2++] = line.first.z();
      lineBufferArray[i2++] = line.second.x();
      lineBufferArray[i2++] = line.second.y();
      lineBufferArray[i2++] = line.second.z();
    }
    // Line loop buffer
    lineGeometry = new Qt3DRender::QGeometry();
    lineGeometry->setObjectName("lineGeometry");
    lineBuffer = new Qt3DRender::QBuffer(lineGeometry);
    lineBuffer->setObjectName("Line buffer");
    lineBuffer->setData(lineByteArray);

    // Line attribute
    lineAtt = new Qt3DRender::QAttribute;
    lineAtt->setObjectName("Position attribute");
    lineAtt->setName(Qt3DRender::QAttribute::defaultPositionAttributeName());
    lineAtt->setBuffer(lineBuffer);
    lineAtt->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
    lineAtt->setVertexBaseType(BASETYPE);
    lineAtt->setVertexSize(3);
    lineAtt->setCount((G4int)nLines);
    lineAtt->setByteOffset(0);
    lineAtt->setByteStride(vertexByteSize);
  }

  // Create material and renderer(s)...

  const auto& colour = fpVisAttribs->GetColour();
  Qt3DExtras::QDiffuseSpecularMaterial* material;
  Qt3DRender::QGeometryRenderer* renderer;
  switch (drawing_style) {
      
    case G4ViewParameters::wireframe:

      lineGeometry->addAttribute(lineAtt);

      material = new Qt3DExtras::QDiffuseSpecularMaterial();
      material->setObjectName("materialForWireframe");
      material->setAmbient(G4Qt3DUtils::ConvertToQColor(colour));
      material->setShininess(0.);
      material->setSpecular(0.);
      wireframeEntity->addComponent(material);

      renderer = new Qt3DRender::QGeometryRenderer;
      renderer->setObjectName("polyhedronWireframeRenderer");
      renderer->setGeometry(lineGeometry);
      renderer->setVertexCount(2*(G4int)nLines);
      renderer->setPrimitiveType(Qt3DRender::QGeometryRenderer::Lines);
      wireframeEntity->addComponent(renderer);

      break;

    case G4ViewParameters::hlr:

      // Surfaces with background colour to hide the edges

      vertexGeometry->addAttribute(positionAtt);
      vertexGeometry->addAttribute(normalAtt);

      material = new Qt3DExtras::QDiffuseSpecularMaterial();
      material->setObjectName("materialForHiddenLines");
      material->setAmbient(Qt::white);  // White for now (should be from fVP)
      material->setShininess(0.);
      material->setSpecular(0.);
      surfaceEntity->addComponent(material);

      renderer = new Qt3DRender::QGeometryRenderer;
      renderer->setObjectName("polyhedronSurfaceRenderer");
      renderer->setGeometry(vertexGeometry);
      renderer->setVertexCount((G4int)nVerts);
      renderer->setPrimitiveType(Qt3DRender::QGeometryRenderer::Triangles);
      surfaceEntity->addComponent(renderer);

      // Edges

      lineGeometry->addAttribute(lineAtt);

      material = new Qt3DExtras::QDiffuseSpecularMaterial();
      material->setObjectName("materialForWireFrame");
      material->setAmbient(G4Qt3DUtils::ConvertToQColor(colour));
      material->setShininess(0.);
      material->setSpecular(0.);
      wireframeEntity->addComponent(material);

      renderer = new Qt3DRender::QGeometryRenderer;
      renderer->setObjectName("polyhedronWireframeRenderer");
      renderer->setGeometry(lineGeometry);
      renderer->setVertexCount(2*(G4int)nLines);
      renderer->setPrimitiveType(Qt3DRender::QGeometryRenderer::Lines);
      wireframeEntity->addComponent(renderer);

      break;

    case G4ViewParameters::hsr:

      vertexGeometry->addAttribute(positionAtt);
      vertexGeometry->addAttribute(normalAtt);

      material = new Qt3DExtras::QDiffuseSpecularMaterial();
      material->setObjectName("materialForSurface");
      material->setAmbient(G4Qt3DUtils::ConvertToQColor(colour));
      if (colour.GetAlpha() < 1.) material->setAlphaBlendingEnabled(true);
      surfaceEntity->addComponent(material);

      renderer = new Qt3DRender::QGeometryRenderer;
      renderer->setObjectName("polyhedronSurfaceRenderer");
      renderer->setGeometry(vertexGeometry);
      renderer->setVertexCount((G4int)nVerts);
      renderer->setPrimitiveType(Qt3DRender::QGeometryRenderer::Triangles);
      surfaceEntity->addComponent(renderer);

      break;

    case G4ViewParameters::hlhsr:

      // Surfaces

      vertexGeometry->addAttribute(positionAtt);
      vertexGeometry->addAttribute(normalAtt);

      material = new Qt3DExtras::QDiffuseSpecularMaterial();
      material->setObjectName("materialForSurface");
      material->setAmbient(G4Qt3DUtils::ConvertToQColor(colour));
      if (colour.GetAlpha() < 1.) material->setAlphaBlendingEnabled(true);
      surfaceEntity->addComponent(material);

      renderer = new Qt3DRender::QGeometryRenderer;
      renderer->setObjectName("polyhedronSurfaceRenderer");
      renderer->setGeometry(vertexGeometry);
      renderer->setVertexCount((G4int)nVerts);
      renderer->setPrimitiveType(Qt3DRender::QGeometryRenderer::Triangles);
      surfaceEntity->addComponent(renderer);

      // Edges

      lineGeometry->addAttribute(lineAtt);

      material = new Qt3DExtras::QDiffuseSpecularMaterial();
      material->setObjectName("materialForWireframe");
      material->setAmbient(G4Qt3DUtils::ConvertToQColor(colour));
      material->setShininess(0.);
      material->setSpecular(0.);
      wireframeEntity->addComponent(material);

      renderer = new Qt3DRender::QGeometryRenderer;
      renderer->setObjectName("polyhedronSurfaceRenderer");
      renderer->setGeometry(lineGeometry);
      renderer->setVertexCount(2*(G4int)nLines);
      renderer->setPrimitiveType(Qt3DRender::QGeometryRenderer::Lines);
      wireframeEntity->addComponent(renderer);

      break;

    case G4ViewParameters::cloud:
      // Case trapped at start of function, so no need to implement
      break;
  }
}

void G4Qt3DSceneHandler::AddCompound(const G4Mesh& mesh)
{
  StandardSpecialMeshRendering(mesh);
}

void G4Qt3DSceneHandler::ClearStore ()
{
  G4Qt3DUtils::delete_components_and_children_of_entity_recursively(fpQt3DScene);
  EstablishG4Qt3DQEntities();
}

void G4Qt3DSceneHandler::ClearTransientStore ()
{
  G4Qt3DUtils::delete_components_and_children_of_entity_recursively(fpTransientObjects);
}
