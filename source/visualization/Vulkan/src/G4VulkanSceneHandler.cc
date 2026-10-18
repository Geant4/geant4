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
// Alex Vergara Gil  17th October 2026
// Ported to Vulkan from the OpenGL work by Andrew Walkden (27th March 1996).
// Vulkan stored scene - creates Vulkan display lists.
// Vulkan immediate scene - draws immediately to buffer
//                           (saving space on server).

#include "G4VulkanSceneHandler.hh"

#include "G4AttHolder.hh"
#include "G4Circle.hh"
#include "G4LogicalVolume.hh"
#include "G4Mesh.hh"
#include "G4Normal3D.hh"
#include "G4VulkanTransform3D.hh"
#include "G4VulkanViewer.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4Point3D.hh"
#include "G4Polyhedron.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4PseudoScene.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4Scene.hh"
#include "G4Square.hh"
#include "G4Text.hh"
#include "G4Transform3D.hh"
#include "G4VMarker.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4VisAttributes.hh"
#include "G4VisExtent.hh"
#include "G4VisManager.hh"

const VKubyte G4VulkanSceneHandler::fStippleMaskHashed[128] = {
  0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55,
  0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55,
  0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55,
  0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55,
  0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55,
  0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55,
  0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55,
  0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55};

G4VulkanSceneHandler::G4VulkanSceneHandler(G4VGraphicsSystem& system, G4int id,
                                           const G4String& name)
  : G4VSceneHandler(system, id, name),
    fPickName(0),
    fThreePassCapable(false),
    fSecondPassForTransparencyRequested(false),
    fSecondPassForTransparency(false),
    fThirdPassForNonHiddenMarkersRequested(false),
    fThirdPassForNonHiddenMarkers(false),
    fEdgeFlag(true)
{}

G4VulkanSceneHandler::~G4VulkanSceneHandler()
{
  ClearStore();
}

void G4VulkanSceneHandler::ClearAndDestroyAtts()
{
  std::map<VKuint, G4AttHolder*>::iterator i;
  for (i = fPickMap.begin(); i != fPickMap.end(); ++i)
    delete i->second;
  fPickMap.clear();
}

G4int G4VulkanSceneHandler::fEntitiesFlushInterval = 100;
G4VulkanSceneHandler::FlushAction G4VulkanSceneHandler::fFlushAction =
  G4VulkanSceneHandler::NthEvent;

void G4VulkanSceneHandler::ScaledFlush()
{
  if (fReadyForTransients)
  {
    // Drawing transients, e.g., trajectories.

    if (!fpScene)
    {
      // No scene - shouldn't happen
      vkFlush();
      return;
    }
    // Get event from modeling parameters
    if (!fpModel)
    {
      // No model - shouldn't happen
      vkFlush();
      return;
    }
    const G4ModelingParameters* modelingParameters = fpModel->GetModelingParameters();
    if (!modelingParameters)
    {
      // No modeling parameters - shouldn't happen
      vkFlush();
      return;
    }
    const G4Event* thisEvent = modelingParameters->GetEvent();
    if (!thisEvent)
    {
      // No event, so not in event loop.
      if (fFlushAction == endOfEvent)
      {
        fFlushAction = endOfRun;
      }
      else if (fFlushAction == NthEvent)
      {
        fFlushAction = NthPrimitive;
      }
    }
    G4RunManager* runMan = G4RunManagerFactory::GetMasterRunManager();
    if (!runMan)
    {
      // No run manager - shouldn't happen
      vkFlush();
      return;
    }
    const G4Run* thisRun = runMan->GetCurrentRun();
    if (!thisRun)
    {
      // No run, so not in event loop.
      if (fFlushAction == endOfRun)
      {
        fFlushAction = NthPrimitive;
      }
      else if (fFlushAction == NthEvent)
      {
        fFlushAction = NthPrimitive;
      }
    }

    switch (fFlushAction)
    {
      case endOfEvent:
        // If "/vis/scene/endOfEventAction refresh", primitives are flushed at
        // end of run anyway, so only scale if false.
        if (!fpScene->GetRefreshAtEndOfEvent())
        {
          // But if "/vis/scene/endOfEventAction accumulate", ShowView is not
          // called until end of run, so we have to watch for a new event.
          // Get event from modeling parameters
          G4int thisEventID = thisEvent->GetEventID();
          static G4int lastEventID = 0;
          if (thisEventID != lastEventID)
          {
            vkFlush();
            lastEventID = thisEventID;
          }
        }
        break;
      case endOfRun:
        // If "/vis/scene/endOfRunAction refresh", primitives are flushed at
        // end of run anyway, so only scale if false.
        if (!fpScene->GetRefreshAtEndOfRun())
        {
          // If "/vis/scene/endOfRunAction accumulate", ShowView is never called
          // so we have to watch for a new run.
          G4int thisRunID = thisRun->GetRunID();
          static G4int lastRunID = 0;
          if (thisRunID != lastRunID)
          {
            vkFlush();
            lastRunID = thisRunID;
          }
        }
        break;
      case eachPrimitive:
        // This is equivalent to numeric with fEntitiesFlushInterval == 1.
        fEntitiesFlushInterval = 1;
        [[fallthrough]];  // Fall through to NthPrimitive.
      case NthPrimitive: {  // Encapsulate in scope {} brackets to satisfy Windows.
        static G4int primitivesWaitingToBeFlushed = 0;
        primitivesWaitingToBeFlushed++;
        if (primitivesWaitingToBeFlushed < fEntitiesFlushInterval) return;
        vkFlush();
        primitivesWaitingToBeFlushed = 0;
        break;
      }
      case NthEvent:
        // If "/vis/scene/endOfEventAction refresh", primitives are flushed at
        // end of event anyway, so only scale if false.
        if (!fpScene->GetRefreshAtEndOfEvent())
        {
          G4int thisEventID = thisEvent->GetEventID();
          static G4int lastEventID = 0;
          if (thisEventID != lastEventID)
          {
            static G4int eventsWaitingToBeFlushed = 0;
            eventsWaitingToBeFlushed++;
            if (eventsWaitingToBeFlushed < fEntitiesFlushInterval) return;
            vkFlush();
            eventsWaitingToBeFlushed = 0;
            lastEventID = thisEventID;
          }
        }
        break;
      case never:
        break;
      default:
        break;
    }
  }

  else

  {
    // For run duration model drawing (detector drawing):
    // Immediate mode: a huge speed up is obtained if flushes are scaled.
    // Stored mode: no discernable difference since drawing is done to the
    //   back buffer and then swapped.
    // So eachPrimitive and NthPrimitive make sense.  But endOfEvent and
    // endOfRun are treated as "no action", i.e., a flush will only be issued,
    // as happens anyway, when drawing is complete.

    switch (fFlushAction)
    {
      case endOfEvent:
        break;
      case endOfRun:
        break;
      case eachPrimitive:
        // This is equivalent to NthPrimitive with fEntitiesFlushInterval == 1.
        fEntitiesFlushInterval = 1;
        [[fallthrough]];  // Fall through to NthPrimitive.
      case NthPrimitive: {  // Encapsulate in scope {} brackets to satisfy Windows.
        static G4int primitivesWaitingToBeFlushed = 0;
        primitivesWaitingToBeFlushed++;
        if (primitivesWaitingToBeFlushed < fEntitiesFlushInterval) return;
        vkFlush();
        primitivesWaitingToBeFlushed = 0;
        break;
      }
      case NthEvent:
        break;
      case never:
        break;
      default:
        break;
    }
  }
}

void G4VulkanSceneHandler::ProcessScene()
{
  fThreePassCapable = true;

  G4VSceneHandler::ProcessScene();

  // Repeat if required...
  if (fSecondPassForTransparencyRequested)
  {
    fSecondPassForTransparency = true;
    G4VSceneHandler::ProcessScene();
    fSecondPassForTransparency = false;
    fSecondPassForTransparencyRequested = false;
  }

  // And again if required...
  if (fThirdPassForNonHiddenMarkersRequested)
  {
    fThirdPassForNonHiddenMarkers = true;
    G4VSceneHandler::ProcessScene();
    fThirdPassForNonHiddenMarkers = false;
    fThirdPassForNonHiddenMarkersRequested = false;
  }

  fThreePassCapable = false;
}

void G4VulkanSceneHandler::PreAddSolid(const G4Transform3D& objectTransformation,
                                       const G4VisAttributes& visAttribs)
{
  G4VSceneHandler::PreAddSolid(objectTransformation, visAttribs);
}

void G4VulkanSceneHandler::BeginPrimitives(const G4Transform3D& objectTransformation)
{
  vkClearColor(0.0, 0.0, 0.0, 0.0);
  vkClearDepth(1.0);
  vkDisable(VK_LINE_SMOOTH);
  vkDisable(VK_POLYGON_SMOOTH);
  vkDisable(VK_POINT_SMOOTH);

  G4VSceneHandler::BeginPrimitives(objectTransformation);
}

void G4VulkanSceneHandler::EndPrimitives()
{
  G4VSceneHandler::EndPrimitives();
}

void G4VulkanSceneHandler::BeginPrimitives2D(const G4Transform3D& objectTransformation)
{
  vkClearColor(0.0, 0.0, 0.0, 0.0);
  vkClearDepth(1.0);
  vkDisable(VK_LINE_SMOOTH);
  vkDisable(VK_POLYGON_SMOOTH);
  vkDisable(VK_POINT_SMOOTH);

  G4VSceneHandler::BeginPrimitives2D(objectTransformation);
}

void G4VulkanSceneHandler::EndPrimitives2D()
{
  G4VSceneHandler::EndPrimitives2D();
}

G4DisplacedSolid* G4VulkanSceneHandler::CreateSectionSolid()
{
  return G4VSceneHandler::CreateSectionSolid();
  // If clipping done in G4VulkanViewer::SetView
  // return 0;
  // Note: if you change this, you must also change
  // G4VulkanStoredViewer::CompareForKernelVisit
}

G4DisplacedSolid* G4VulkanSceneHandler::CreateCutawaySolid()
{
  // return G4VSceneHandler::CreateCutawaySolid();
  // If cutaway done in G4VulkanViewer::SetView.
  return 0;
  // Note: if you change this, you must also change
  // G4VulkanStoredViewer::CompareForKernelVisit
}

void G4VulkanSceneHandler::AddPrimitive(const G4Polyline& line)
{
  std::size_t nPoints = line.size();
  if (nPoints <= 0) return;

  // Note: colour and depth test treated in sub-class.

  vkDisable(VK_LIGHTING);

  G4double lineWidth = GetLineWidth(fpVisAttribs);
  // Need access to method in G4VulkanViewer.  static_cast doesn't
  // work with a virtual base class, so use dynamic_cast.  No need to
  // test the outcome since viewer is guaranteed to be a
  // G4VulkanViewer, but test it anyway to keep Coverity happy.
  G4VulkanViewer* pVKViewer = dynamic_cast<G4VulkanViewer*>(fpViewer);
  if (pVKViewer) pVKViewer->ChangeLineWidth(lineWidth);

  fEdgeFlag = true;
  vkBegin(VK_LINE_STRIP);
  // No ned vkEdgeFlag for lines :
  // Boundary and nonboundary edge flags on vertices are significant only if VK_POLYGON_MODE is set
  // to VK_POINT or VK_LINE.  See vkPolygonMode.

  //  vkEdgeFlag (VK_TRUE);
  for (std::size_t iPoint = 0; iPoint < nPoints; ++iPoint)
  {
    G4double x, y, z;
    x = line[iPoint].x();
    y = line[iPoint].y();
    z = line[iPoint].z();
    vkVertex3d(x, y, z);
  }
  vkEnd();
}

void G4VulkanSceneHandler::AddPrimitive(const G4Polymarker& polymarker)
{
  if (polymarker.size() == 0)
  {
    return;
  }

  // Note: colour and depth test treated in sub-class.

  vkDisable(VK_LIGHTING);

  MarkerSizeType sizeType;
  G4double size = GetMarkerSize(polymarker, sizeType);

  // Need access to method in G4VulkanViewer.  static_cast doesn't
  // work with a virtual base class, so use dynamic_cast.  No need to
  // test the outcome since viewer is guaranteed to be a
  // G4VulkanViewer, but test it anyway to keep Coverity happy.
  G4VulkanViewer* pVKViewer = dynamic_cast<G4VulkanViewer*>(fpViewer);
  if (!pVKViewer) return;

  if (sizeType == world)
  {  // Size specified in world coordinates.
    G4double lineWidth = GetLineWidth(fpVisAttribs);
    pVKViewer->ChangeLineWidth(lineWidth);

    G4VMarker::FillStyle style = polymarker.GetFillStyle();

    // G4bool filled = false;  Not actually used - comment out to prevent compiler warnings (JA).
    static G4bool hashedWarned = false;

    switch (style)
    {
      case G4VMarker::noFill:
        vkPolygonMode(VK_FRONT_AND_BACK, VK_LINE);
        vkEdgeFlag(VK_TRUE);
        // filled = false;
        break;
      case G4VMarker::hashed:
        if (!hashedWarned)
        {
          G4cout << "Hashed fill style in G4VulkanSceneHandler."
                 << "\n  Not implemented.  Using G4VMarker::filled." << G4endl;
          hashedWarned = true;
        }
        // Maybe use
        // vkPolygonStipple (fStippleMaskHashed);
        [[fallthrough]];  // Drop through to filled...
      case G4VMarker::filled:
        vkPolygonMode(VK_FRONT_AND_BACK, VK_FILL);
        // filled = true;
        break;
    }
  }

  // Draw...
  if (sizeType == world)
  {  // Size specified in world coordinates.

    G4int nSides;
    G4double startPhi;
    switch (polymarker.GetMarkerType())
    {
      default:
      case G4Polymarker::dots:
        [[fallthrough]];  // Fall through to circles
      case G4Polymarker::circles:
        nSides = GetNoOfSides(fpVisAttribs);
        startPhi = 0.;
        break;
      case G4Polymarker::squares:
        nSides = 4;
        startPhi = -pi / 4.;
        break;
    }

    const G4Vector3D& viewpointDirection = fpViewer->GetViewParameters().GetViewpointDirection();
    const G4Vector3D& up = fpViewer->GetViewParameters().GetUpVector();
    const G4double dPhi = twopi / nSides;
    const G4double radius = size / 2.;
    G4Vector3D start = radius * (up.cross(viewpointDirection)).unit();
    G4double phi;
    G4int i;
    for (size_t iPoint = 0; iPoint < polymarker.size(); iPoint++)
    {
      fEdgeFlag = true;
      vkBegin(VK_POLYGON);
      for (i = 0, phi = startPhi; i < nSides; i++, phi += dPhi)
      {
        G4Vector3D r = start;
        r.rotate(phi, viewpointDirection);
        G4Vector3D p = polymarker[iPoint] + r;
        vkVertex3d(p.x(), p.y(), p.z());
      }
      vkEnd();
    }
  }
  else
  {  // Size specified in screen (window) coordinates.

    pVKViewer->ChangePointSize(size);

    switch (polymarker.GetMarkerType())
    {
      default:
      case G4Polymarker::dots:
        [[fallthrough]];  // Fall through to circles
      case G4Polymarker::circles:
        if (fpViewer->GetViewParameters().IsDotsSmooth())
        {
          vkEnable(VK_POINT_SMOOTH);
        }
        break;
      case G4Polymarker::squares:
        break;
    }
    vkBegin(VK_POINTS);
    for (size_t iPoint = 0; iPoint < polymarker.size(); iPoint++)
    {
      G4Point3D centre = polymarker[iPoint];
      vkVertex3d(centre.x(), centre.y(), centre.z());
    }
    vkEnd();
  }
}

void G4VulkanSceneHandler::AddPrimitive(const G4Text& text)
{
  // Pass to specific viewer via virtual function DrawText.
  G4VulkanViewer* pVKViewer = dynamic_cast<G4VulkanViewer*>(fpViewer);
  if (pVKViewer) pVKViewer->DrawText(text);
}

void G4VulkanSceneHandler::AddPrimitive(const G4Circle& circle)
{
  G4Polymarker oneCircle(circle);
  oneCircle.push_back(circle.GetPosition());
  oneCircle.SetMarkerType(G4Polymarker::circles);
  // Call this AddPrimitive to avoid re-doing sub-class code.
  G4VulkanSceneHandler::AddPrimitive(oneCircle);
}

void G4VulkanSceneHandler::AddPrimitive(const G4Square& square)
{
  G4Polymarker oneSquare(square);
  oneSquare.push_back(square.GetPosition());
  oneSquare.SetMarkerType(G4Polymarker::squares);
  // Call this AddPrimitive to avoid re-doing sub-class code.
  G4VulkanSceneHandler::AddPrimitive(oneSquare);
}

// Method for handling G4Polyhedron objects for drawing solids.
void G4VulkanSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron)
{
  // Assume all facets are planar convex quadrilaterals.
  // Draw each facet individually

  if (polyhedron.GetNoFacets() == 0) return;

  // Need access to data in G4VulkanViewer.  static_cast doesn't work
  // with a virtual base class, so use dynamic_cast.
  G4VulkanViewer* pVKViewer = dynamic_cast<G4VulkanViewer*>(fpViewer);
  if (!pVKViewer) return;

  // Get view parameters that the user can force through the vis
  // attributes, thereby over-riding the current view parameter.
  G4ViewParameters::DrawingStyle drawing_style = GetDrawingStyle(fpVisAttribs);

  // Note that in stored mode, because this call gets embedded in a display
  //  list, it is the colour _at the time of_ creation of the display list, so
  //  even if the colour is changed, for example, by interaction with a Qt
  //  window, current_colour does not change.
  VKfloat* painting_colour;
  VKfloat clear_colour[4];
  VKfloat current_colour[4];
  vkGetFloatv(VK_CURRENT_COLOR, current_colour);

  G4bool isTransparent = false;
  if (current_colour[3] < 1.)
  {  // This object is transparent
    isTransparent = true;
  }

  if (drawing_style == G4ViewParameters::hlr)
  {
    // This is the colour used to paint surfaces in hlr mode.
    vkGetFloatv(VK_COLOR_CLEAR_VALUE, clear_colour);
    painting_colour = clear_colour;
  }
  else
  {  // drawing_style == G4ViewParameters::hlhsr
    painting_colour = current_colour;
  }

  G4double lineWidth = GetLineWidth(fpVisAttribs);
  pVKViewer->ChangeLineWidth(lineWidth);

  G4bool isAuxEdgeVisible = GetAuxEdgeVisible(fpVisAttribs);

  G4bool clipping = pVKViewer->fVP.IsSection() || pVKViewer->fVP.IsCutaway();

  // Lighting disabled unless otherwise requested
  vkDisable(VK_LIGHTING);

  switch (drawing_style)
  {
    case (G4ViewParameters::hlhsr):
      // Set up as for hidden line removal but paint polygon faces later...
    case (G4ViewParameters::hlr):
      vkEnable(VK_STENCIL_TEST);
      // The stencil buffer is cleared in G4VulkanViewer::ClearView.
      // The procedure below leaves it clear.
      vkStencilFunc(VK_ALWAYS, 0, 1);
      vkStencilOp(VK_INVERT, VK_INVERT, VK_INVERT);
      vkEnable(VK_DEPTH_TEST);
      vkDepthFunc(VK_LEQUAL);
      if (isTransparent)
      {
        // Transparent...
        vkColorMaterial(VK_FRONT_AND_BACK, VK_AMBIENT_AND_DIFFUSE);
        vkEnable(VK_COLOR_MATERIAL);
        // vkDisable (VK_CULL_FACE);
        vkPolygonMode(VK_FRONT_AND_BACK, VK_LINE);
      }
      else
      {
        // Opaque...
        if (clipping)
        {
          vkColorMaterial(VK_FRONT_AND_BACK, VK_AMBIENT_AND_DIFFUSE);
          vkEnable(VK_COLOR_MATERIAL);
          // vkDisable (VK_CULL_FACE);
          vkPolygonMode(VK_FRONT_AND_BACK, VK_LINE);
        }
        else
        {
          vkColorMaterial(VK_FRONT, VK_AMBIENT_AND_DIFFUSE);
          vkEnable(VK_COLOR_MATERIAL);
          // vkEnable (VK_CULL_FACE);
          // vkCullFace (VK_BACK);
          vkPolygonMode(VK_FRONT, VK_LINE);
        }
      }
      break;
    case (G4ViewParameters::hsr):
      vkEnable(VK_DEPTH_TEST);
      vkDepthFunc(VK_LEQUAL);
      if (isTransparent)
      {
        // Transparent...
        vkDepthMask(VK_FALSE);  // Make depth buffer read-only.
        vkColorMaterial(VK_FRONT_AND_BACK, VK_AMBIENT_AND_DIFFUSE);
        vkEnable(VK_COLOR_MATERIAL);
        // vkDisable (VK_CULL_FACE);
        vkPolygonMode(VK_FRONT_AND_BACK, VK_FILL);
      }
      else
      {
        // Opaque...
        vkDepthMask(VK_TRUE);  // Make depth buffer writable (default).
        if (clipping)
        {
          vkColorMaterial(VK_FRONT_AND_BACK, VK_AMBIENT_AND_DIFFUSE);
          vkEnable(VK_COLOR_MATERIAL);
          // vkDisable (VK_CULL_FACE);
          vkPolygonMode(VK_FRONT_AND_BACK, VK_FILL);
        }
        else
        {
          vkColorMaterial(VK_FRONT, VK_AMBIENT_AND_DIFFUSE);
          vkEnable(VK_COLOR_MATERIAL);
          // vkEnable (VK_CULL_FACE);
          // vkCullFace (VK_BACK);
          vkPolygonMode(VK_FRONT, VK_FILL);
        }
      }
      if (!fProcessing2D) vkEnable(VK_LIGHTING);
      break;
    case (G4ViewParameters::wireframe):
    default:
      vkEnable(VK_DEPTH_TEST);
      vkDepthFunc(VK_LEQUAL);  //??? was VK_ALWAYS
      // vkDisable (VK_CULL_FACE);
      vkPolygonMode(VK_FRONT_AND_BACK, VK_LINE);
      break;
  }

  // Loop through all the facets...
  fEdgeFlag = true;
  vkBegin(VK_QUADS);
  vkEdgeFlag(VK_TRUE);
  G4bool notLastFace;
  do
  {
    // First, find vertices, edgeflags and normals and note "not last facet"...
    G4Point3D vertex[4];
    G4int edgeFlag[4];
    G4Normal3D normals[4];
    G4int nEdges;
    notLastFace = polyhedron.GetNextFacet(nEdges, vertex, edgeFlag, normals);

    // Loop through the four edges of each G4Facet...
    for (G4int edgeCount = 0; edgeCount < nEdges; ++edgeCount)
    {
      // Check to see if edge is visible or not...
      if (isAuxEdgeVisible)
      {
        edgeFlag[edgeCount] = 1;
      }
      if (edgeFlag[edgeCount] > 0)
      {
        if (fEdgeFlag != true)
        {
          vkEdgeFlag(VK_TRUE);
          fEdgeFlag = true;
        }
      }
      else
      {
        if (fEdgeFlag != false)
        {
          vkEdgeFlag(VK_FALSE);
          fEdgeFlag = false;
        }
      }
      vkNormal3d(normals[edgeCount].x(), normals[edgeCount].y(), normals[edgeCount].z());
      vkVertex3d(vertex[edgeCount].x(), vertex[edgeCount].y(), vertex[edgeCount].z());
    }

    // HepPolyhedron produces triangles too; in that case add an extra
    // vertex identical to first...
    if (nEdges == 3)
    {
      G4int edgeCount = 3;
      normals[edgeCount] = normals[0];
      vertex[edgeCount] = vertex[0];
      edgeFlag[edgeCount] = -1;
      if (fEdgeFlag != false)
      {
        vkEdgeFlag(VK_FALSE);
        fEdgeFlag = false;
      }

      vkNormal3d(normals[edgeCount].x(), normals[edgeCount].y(), normals[edgeCount].z());
      vkVertex3d(vertex[edgeCount].x(), vertex[edgeCount].y(), vertex[edgeCount].z());
    }
    // Trap situation where number of edges is > 4...
    if (nEdges > 4)
    {
      G4cerr << "G4VulkanSceneHandler::AddPrimitive(G4Polyhedron): WARNING"
                "\n   G4Polyhedron facet with "
             << nEdges << " edges" << G4endl;
    }

    // Do it all over again (twice) for hlr...
    if (drawing_style == G4ViewParameters::hlr || drawing_style == G4ViewParameters::hlhsr)
    {
      vkDisable(VK_COLOR_MATERIAL);  // Revert to vkMaterial for hlr/sr.
      vkEnd();  // Placed here to balance vkBegin above, allowing GL

      // state changes below, then vkBegin again.  Avoids
      // having vkBegin/End pairs *inside* loop in the more
      // usual case of no hidden line removal.

      // Lighting disabled unless otherwise requested
      vkDisable(VK_LIGHTING);

      // Draw through stencil...
      vkStencilFunc(VK_EQUAL, 0, 1);
      vkStencilOp(VK_KEEP, VK_KEEP, VK_KEEP);
      if (drawing_style == G4ViewParameters::hlhsr)
      {
        if (!fProcessing2D) vkEnable(VK_LIGHTING);
      }
      vkEnable(VK_DEPTH_TEST);
      vkDepthFunc(VK_LEQUAL);
      if (isTransparent)
      {
        // Transparent...
        vkDepthMask(VK_FALSE);  // Make depth buffer read-only.
        // vkDisable (VK_CULL_FACE);
        vkPolygonMode(VK_FRONT_AND_BACK, VK_FILL);
      }
      else
      {
        // Opaque...
        vkDepthMask(VK_TRUE);  // Make depth buffer writable (default).
        if (clipping)
        {
          // vkDisable (VK_CULL_FACE);
          vkPolygonMode(VK_FRONT_AND_BACK, VK_FILL);
        }
        else
        {
          // vkEnable (VK_CULL_FACE);
          // vkCullFace (VK_BACK);
          vkPolygonMode(VK_FRONT, VK_FILL);
        }
      }
      if (drawing_style == G4ViewParameters::hlr)
      {
        if (isTransparent)
        {
          // Transparent - don't paint...
          goto end_of_drawing_through_stencil;
        }
      }
      if (isTransparent)
      {
        // Transparent...
        vkMaterialfv(VK_FRONT_AND_BACK, VK_AMBIENT_AND_DIFFUSE, painting_colour);
      }
      else
      {
        // Opaque...
        vkMaterialfv(VK_FRONT, VK_AMBIENT_AND_DIFFUSE, painting_colour);
      }
      vkColor4fv(painting_colour);
      vkBegin(VK_QUADS);
      vkEdgeFlag(VK_TRUE);
      fEdgeFlag = true;

      for (int edgeCount = 0; edgeCount < 4; ++edgeCount)
      {
        if (edgeFlag[edgeCount] > 0)
        {
          if (fEdgeFlag != true)
          {
            vkEdgeFlag(VK_TRUE);
            fEdgeFlag = true;
          }
        }
        else
        {
          if (fEdgeFlag != false)
          {
            vkEdgeFlag(VK_FALSE);
            fEdgeFlag = false;
          }
        }
        vkNormal3d(normals[edgeCount].x(), normals[edgeCount].y(), normals[edgeCount].z());
        vkVertex3d(vertex[edgeCount].x(), vertex[edgeCount].y(), vertex[edgeCount].z());
      }
      vkEnd();
    end_of_drawing_through_stencil:

      // and once more to reset the stencil bits...
      vkStencilFunc(VK_ALWAYS, 0, 1);
      vkStencilOp(VK_INVERT, VK_INVERT, VK_INVERT);
      vkDepthFunc(VK_LEQUAL);  // to make sure line gets drawn.
      if (isTransparent)
      {
        // Transparent...
        // vkDisable (VK_CULL_FACE);
        vkPolygonMode(VK_FRONT_AND_BACK, VK_LINE);
      }
      else
      {
        // Opaque...
        if (clipping)
        {
          // vkDisable (VK_CULL_FACE);
          vkPolygonMode(VK_FRONT_AND_BACK, VK_LINE);
        }
        else
        {
          // vkEnable (VK_CULL_FACE);
          // vkCullFace (VK_BACK);
          vkPolygonMode(VK_FRONT, VK_LINE);
        }
      }
      vkDisable(VK_LIGHTING);
      vkColor4fv(current_colour);
      fEdgeFlag = true;
      vkBegin(VK_QUADS);
      vkEdgeFlag(VK_TRUE);
      fEdgeFlag = true;
      for (int edgeCount = 0; edgeCount < 4; ++edgeCount)
      {
        if (edgeFlag[edgeCount] > 0)
        {
          if (fEdgeFlag != true)
          {
            vkEdgeFlag(VK_TRUE);
            fEdgeFlag = true;
          }
        }
        else
        {
          if (fEdgeFlag != false)
          {
            vkEdgeFlag(VK_FALSE);
            fEdgeFlag = false;
          }
        }
        vkNormal3d(normals[edgeCount].x(), normals[edgeCount].y(), normals[edgeCount].z());
        vkVertex3d(vertex[edgeCount].x(), vertex[edgeCount].y(), vertex[edgeCount].z());
      }
      vkEnd();

      vkDepthFunc(VK_LEQUAL);  // Revert for next facet.
      fEdgeFlag = true;
      vkBegin(VK_QUADS);  // Ready for next facet.  GL
      vkEdgeFlag(VK_TRUE);
      fEdgeFlag = true;
      // says it ignores incomplete
      // quadrilaterals, so final empty
      // vkBegin/End sequence should be OK.
    }
  } while (notLastFace);

  vkEnd();
  vkDisable(VK_STENCIL_TEST);  // Revert to default for next primitive.
  vkDepthMask(VK_TRUE);  // Revert to default for next primitive.
  vkDisable(VK_LIGHTING);  // Revert to default for next primitive.
}

void G4VulkanSceneHandler::AddCompound(const G4VTrajectory& traj)
{
  G4VSceneHandler::AddCompound(traj);  // For now.
}

void G4VulkanSceneHandler::AddCompound(const G4VHit& hit)
{
  G4VSceneHandler::AddCompound(hit);  // For now.
}

void G4VulkanSceneHandler::AddCompound(const G4VDigi& digi)
{
  G4VSceneHandler::AddCompound(digi);  // For now.
}

void G4VulkanSceneHandler::AddCompound(const G4THitsMap<G4double>& hits)
{
  G4VSceneHandler::AddCompound(hits);  // For now.
}

void G4VulkanSceneHandler::AddCompound(const G4THitsMap<G4StatDouble>& hits)
{
  G4VSceneHandler::AddCompound(hits);  // For now.
}

void G4VulkanSceneHandler::AddCompound(const G4Mesh& mesh)
{
  StandardSpecialMeshRendering(mesh);
}
