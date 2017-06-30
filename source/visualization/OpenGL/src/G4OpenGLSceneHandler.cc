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
// $Id: G4OpenGLSceneHandler.cc 104288 2017-05-23 13:23:23Z gcosmo $
//
// 
// Andrew Walkden  27th March 1996
// OpenGL stored scene - creates OpenGL display lists.
// OpenGL immediate scene - draws immediately to buffer
//                           (saving space on server).

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#include "G4OpenGLSceneHandler.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLTransform3D.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"
#include "G4Transform3D.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4VMarker.hh"
#include "G4Polyhedron.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Scene.hh"
#include "G4VisExtent.hh"
#include "G4AttHolder.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#endif
#include "G4Run.hh"

const GLubyte G4OpenGLSceneHandler::fStippleMaskHashed [128] = {
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55
};

G4OpenGLSceneHandler::G4OpenGLSceneHandler (G4VGraphicsSystem& system,
                                            G4int id,
                                            const G4String& name):
G4VSceneHandler (system, id, name),
#ifdef G4OPENGL_VERSION_2
fEmulate_GL_QUADS(false),
#endif
fPickName(0),
fThreePassCapable(false),
fSecondPassForTransparencyRequested(false),
fSecondPassForTransparency(false),
fThirdPassForNonHiddenMarkersRequested(false),
fThirdPassForNonHiddenMarkers(false),
fEdgeFlag(true)
{
}

G4OpenGLSceneHandler::~G4OpenGLSceneHandler ()
{
  ClearStore ();
}

void G4OpenGLSceneHandler::ClearAndDestroyAtts()
{
  std::map<GLuint, G4AttHolder*>::iterator i;
  for (i = fPickMap.begin(); i != fPickMap.end(); ++i) delete i->second;
  fPickMap.clear();
}

G4int G4OpenGLSceneHandler::fEntitiesFlushInterval = 100;
G4OpenGLSceneHandler::FlushAction
G4OpenGLSceneHandler::fFlushAction = G4OpenGLSceneHandler::NthEvent;

void G4OpenGLSceneHandler::ScaledFlush()
{
  if (fReadyForTransients) {

    // Drawing transients, e.g., trajectories.

    if (!fpScene) {
      // No scene, so probably not in event loop.
      glFlush();
      return;
    }
    // Get event from modeling parameters
    if (!fpModel) {
      // No model, so probably not in event loop
      glFlush();
      return;
    }
    const G4ModelingParameters* modelingParameters =
    fpModel->GetModelingParameters();
    if (!modelingParameters) {
      // No modeling parameters, so probably not in event loop.
      glFlush();
      return;
    }
    const G4Event* thisEvent = modelingParameters->GetEvent();
    if (!thisEvent) {
      // No event, so probably not in event loop.
      glFlush();
      return;
    }
    G4RunManager* runMan = G4RunManager::GetRunManager();
#ifdef G4MULTITHREADED
    if (G4Threading::IsMultithreadedApplication()) {
      runMan = G4MTRunManager::GetMasterRunManager();
    }
#endif
    if (!runMan) {
      glFlush();
      return;
    }
    const G4Run* thisRun = runMan->GetCurrentRun();
    if (!thisRun) {
      glFlush();
      return;
    }

    switch (fFlushAction) {
      case endOfEvent:
        // If "/vis/scene/endOfEventAction refresh", primitives are flushed at
        // end of run anyway, so only scale if false.
        if (!fpScene->GetRefreshAtEndOfEvent()) {
          // But if "/vis/scene/endOfEventAction accumulate", ShowView is not
          // called until end of run, so we have to watch for a new event.
          // Get event from modeling parameters
          G4int thisEventID = thisEvent->GetEventID();
          static G4int lastEventID = 0;
          if (thisEventID != lastEventID) {
            glFlush();
            lastEventID = thisEventID;
          }
        }
        break;
      case endOfRun:
        // If "/vis/scene/endOfRunAction refresh", primitives are flushed at
        // end of run anyway, so only scale if false.
        if (!fpScene->GetRefreshAtEndOfRun()) {
          // If "/vis/scene/endOfRunAction accumulate", ShowView is never called
          // so we have to watch for a new run.
          G4int thisRunID = thisRun->GetRunID();
          static G4int lastRunID = 0;
          if (thisRunID != lastRunID) {
            glFlush();
            lastRunID = thisRunID;
          }
        }
        break;
      case eachPrimitive:
        // This is equivalent to numeric with fEntitiesFlushInterval == 1.
        fEntitiesFlushInterval = 1;  // fallthrough
        // Fall through to NthPrimitive.
      case NthPrimitive:
      { // Encapsulate in scope {} brackets to satisfy Windows.
        static G4int primitivesWaitingToBeFlushed = 0;
        primitivesWaitingToBeFlushed++;
        if (primitivesWaitingToBeFlushed < fEntitiesFlushInterval) return;
        glFlush();
        primitivesWaitingToBeFlushed = 0;
        break;
      }
      case NthEvent:
        // If "/vis/scene/endOfEventAction refresh", primitives are flushed at
        // end of event anyway, so only scale if false.
        if (!fpScene->GetRefreshAtEndOfEvent()) {
          G4int thisEventID = thisEvent->GetEventID();
          static G4int lastEventID = 0;
          if (thisEventID != lastEventID) {
            static G4int eventsWaitingToBeFlushed = 0;
            eventsWaitingToBeFlushed++;
            if (eventsWaitingToBeFlushed < fEntitiesFlushInterval) return;
            glFlush();
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

    switch (fFlushAction) {
      case endOfEvent:
        break;
      case endOfRun:
        break;
      case eachPrimitive:
        // This is equivalent to NthPrimitive with fEntitiesFlushInterval == 1.
        fEntitiesFlushInterval = 1;  // fallthrough
        // Fall through to NthPrimitive.
      case NthPrimitive:
      { // Encapsulate in scope {} brackets to satisfy Windows.
        static G4int primitivesWaitingToBeFlushed = 0;
        primitivesWaitingToBeFlushed++;
        if (primitivesWaitingToBeFlushed < fEntitiesFlushInterval) return;
        glFlush();
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

void G4OpenGLSceneHandler::ProcessScene()
{
  fThreePassCapable = true;
  
  G4VSceneHandler::ProcessScene();

  // Repeat if required...
  if (fSecondPassForTransparencyRequested) {
    fSecondPassForTransparency = true;
    G4VSceneHandler::ProcessScene();
    fSecondPassForTransparency = false;
    fSecondPassForTransparencyRequested = false;
  }

  // And again if required...
  if (fThirdPassForNonHiddenMarkersRequested) {
    fThirdPassForNonHiddenMarkers = true;
    G4VSceneHandler::ProcessScene();
    fThirdPassForNonHiddenMarkers = false;
    fThirdPassForNonHiddenMarkersRequested = false;
  }
  
  fThreePassCapable = false;
}

void G4OpenGLSceneHandler::PreAddSolid
(const G4Transform3D& objectTransformation,
 const G4VisAttributes& visAttribs)
{
  G4VSceneHandler::PreAddSolid (objectTransformation, visAttribs);
}

void G4OpenGLSceneHandler::BeginPrimitives
(const G4Transform3D& objectTransformation)
{
  G4VSceneHandler::BeginPrimitives (objectTransformation);
}

void G4OpenGLSceneHandler::EndPrimitives ()
{
  G4VSceneHandler::EndPrimitives ();
}

void G4OpenGLSceneHandler::BeginPrimitives2D
(const G4Transform3D& objectTransformation)
{
  G4VSceneHandler::BeginPrimitives2D (objectTransformation);
}

void G4OpenGLSceneHandler::EndPrimitives2D ()
{
  G4VSceneHandler::EndPrimitives2D ();
}

G4VSolid* G4OpenGLSceneHandler::CreateSectionSolid ()
{
  return G4VSceneHandler::CreateSectionSolid();
  // If clipping done in G4OpenGLViewer::SetView
  // return 0;
}

G4VSolid* G4OpenGLSceneHandler::CreateCutawaySolid ()
{
  // Cutaway done in G4OpenGLViewer::SetView.
  return 0;
  // Else
  // return G4VSceneHandler::CreateCutawaySolid();
}

void G4OpenGLSceneHandler::AddPrimitive (const G4Polyline& line)
{
  G4int nPoints = line.size ();
  if (nPoints <= 0) return;

  // Note: colour and depth test treated in sub-class.

#ifndef G4OPENGL_VERSION_2
  glDisable (GL_LIGHTING);
#endif
  
  G4double lineWidth = GetLineWidth(fpVisAttribs);
  // Need access to method in G4OpenGLViewer.  static_cast doesn't
  // work with a virtual base class, so use dynamic_cast.  No need to
  // test the outcome since viewer is guaranteed to be a
  // G4OpenGLViewer, but test it anyway to keep Coverity happy.
  G4OpenGLViewer* pGLViewer = dynamic_cast<G4OpenGLViewer*>(fpViewer);
  if (pGLViewer) pGLViewer->ChangeLineWidth(lineWidth);

  fEdgeFlag = true;
#ifndef G4OPENGL_VERSION_2
  glBegin (GL_LINE_STRIP);
  // No ned glEdgeFlag for lines :
  // Boundary and nonboundary edge flags on vertices are significant only if GL_POLYGON_MODE is set to GL_POINT or GL_LINE.  See glPolygonMode.
  
  //  glEdgeFlag (GL_TRUE);
  for (G4int iPoint = 0; iPoint < nPoints; iPoint++) {
  G4double x, y, z;
    x = line[iPoint].x(); 
    y = line[iPoint].y();
    z = line[iPoint].z();
    glVertex3d (x, y, z);
  }
  glEnd ();
#else
  glBeginVBO(GL_LINE_STRIP);

  for (G4int iPoint = 0; iPoint < nPoints; iPoint++) {
    fOglVertex.push_back(line[iPoint].x());
    fOglVertex.push_back(line[iPoint].y());
    fOglVertex.push_back(line[iPoint].z());
    // normal
    fOglVertex.push_back(0);
    fOglVertex.push_back(0);
    fOglVertex.push_back(1);
  }
  
  glEndVBO();
#endif
}

void G4OpenGLSceneHandler::AddPrimitive (const G4Polymarker& polymarker)
{
  if (polymarker.size() == 0) {
    return;
  }

  // Note: colour and depth test treated in sub-class.

#ifndef G4OPENGL_VERSION_2
  glDisable (GL_LIGHTING);
#endif
  
  MarkerSizeType sizeType;
  G4double size = GetMarkerSize(polymarker, sizeType);

  // Need access to method in G4OpenGLViewer.  static_cast doesn't
  // work with a virtual base class, so use dynamic_cast.  No need to
  // test the outcome since viewer is guaranteed to be a
  // G4OpenGLViewer, but test it anyway to keep Coverity happy.
  G4OpenGLViewer* pGLViewer = dynamic_cast<G4OpenGLViewer*>(fpViewer);
  if (!pGLViewer) return;

  if (sizeType == world) {  // Size specified in world coordinates.
    G4double lineWidth = GetLineWidth(fpVisAttribs);
    pGLViewer->ChangeLineWidth(lineWidth);
  
    G4VMarker::FillStyle style = polymarker.GetFillStyle();
    
    // G4bool filled = false;  Not actually used - comment out to prevent compiler warnings (JA).
    static G4bool hashedWarned = false;
    
    switch (style) {
      case G4VMarker::noFill:
        glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
        glEdgeFlag (GL_TRUE);
        //filled = false;
        break;
      case G4VMarker::hashed:
        if (!hashedWarned) {  // fallthrough
          G4cout << "Hashed fill style in G4OpenGLSceneHandler."
          << "\n  Not implemented.  Using G4VMarker::filled."
          << G4endl;
          hashedWarned = true;
        }  // fallthrough
        // Maybe use
        //glPolygonStipple (fStippleMaskHashed);
        // Drop through to filled...
      case G4VMarker::filled:
        glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
        //filled = true;
        break;
    }
  }
  
  // Draw...
  if (sizeType == world) {  // Size specified in world coordinates.

    G4int nSides;
    G4double startPhi;
    switch (polymarker.GetMarkerType()) {
    default:
    case G4Polymarker::dots:
        size = 1.;  // fallthrough
      // Drop through to circles
    case G4Polymarker::circles:
      nSides = GetNoOfSides(fpVisAttribs);
      startPhi = 0.;
      break;
    case G4Polymarker::squares:
      nSides = 4;
      startPhi = -pi / 4.;
      break;
    }

    const G4Vector3D& viewpointDirection =
      fpViewer -> GetViewParameters().GetViewpointDirection();
    const G4Vector3D& up = fpViewer->GetViewParameters().GetUpVector();
    const G4double dPhi = twopi / nSides;
    const G4double radius = size / 2.;
    G4Vector3D start = radius * (up.cross(viewpointDirection)).unit();
    G4double phi;
    G4int i;
    for (size_t iPoint = 0; iPoint < polymarker.size (); iPoint++) {
      fEdgeFlag = true;
#ifndef G4OPENGL_VERSION_2
      glBegin (GL_POLYGON);
      for (i = 0, phi = startPhi; i < nSides; i++, phi += dPhi) {
	G4Vector3D r = start; r.rotate(phi, viewpointDirection);
	G4Vector3D p = polymarker[iPoint] + r;
	glVertex3d (p.x(), p.y(), p.z());
      }
      glEnd ();
#else
      glBeginVBO (GL_TRIANGLE_STRIP);
      for (i = 0, phi = startPhi; i < nSides; i++, phi += dPhi) {
        G4Vector3D r = start; r.rotate(phi, viewpointDirection);
        G4Vector3D p = polymarker[iPoint] + r;

        fOglVertex.push_back(p.x());
        fOglVertex.push_back(p.y());
        fOglVertex.push_back(p.z());
        // normal
        fOglVertex.push_back(0);
        fOglVertex.push_back(0);
        fOglVertex.push_back(1);
      }
      glEndVBO ();
#endif
    }

  } else { // Size specified in screen (window) coordinates.

    pGLViewer->ChangePointSize(size);

    //Antialiasing only for circles
#ifndef G4OPENGL_VERSION_2
    switch (polymarker.GetMarkerType()) {
    default:
    case G4Polymarker::dots:
    case G4Polymarker::circles:
      glEnable (GL_POINT_SMOOTH); break;
    case G4Polymarker::squares:
      glDisable (GL_POINT_SMOOTH); break;
    }
#endif
#ifndef G4OPENGL_VERSION_2
    glBegin (GL_POINTS);
    for (size_t iPoint = 0; iPoint < polymarker.size (); iPoint++) {
      G4Point3D centre = polymarker[iPoint];
      glVertex3d(centre.x(),centre.y(),centre.z());
    }
    glEnd();     
#else
    glBeginVBO(GL_POINTS);

    for (size_t iPoint = 0; iPoint < polymarker.size (); iPoint++) {
      fOglVertex.push_back(polymarker[iPoint].x());
      fOglVertex.push_back(polymarker[iPoint].y());
      fOglVertex.push_back(polymarker[iPoint].z());
      fOglVertex.push_back(0);
      fOglVertex.push_back(0);
      fOglVertex.push_back(1);
    }
    glEndVBO();
#endif
  }
}

void G4OpenGLSceneHandler::AddPrimitive (const G4Text& text) {
  // Pass to specific viewer via virtual function DrawText.
  // FIXME : Not ready for OPENGL2 for the moment
#ifdef G4OPENGL_VERSION_2
  return;
#endif
  G4OpenGLViewer* pGLViewer = dynamic_cast<G4OpenGLViewer*>(fpViewer);
  if (pGLViewer) pGLViewer->DrawText(text);
}

void G4OpenGLSceneHandler::AddPrimitive (const G4Circle& circle) {
  G4Polymarker oneCircle(circle);
  oneCircle.push_back(circle.GetPosition());
  oneCircle.SetMarkerType(G4Polymarker::circles);
  // Call this AddPrimitive to avoid re-doing sub-class code.
  G4OpenGLSceneHandler::AddPrimitive(oneCircle);
}

void G4OpenGLSceneHandler::AddPrimitive (const G4Square& square) {
  G4Polymarker oneSquare(square);
  oneSquare.push_back(square.GetPosition());
  oneSquare.SetMarkerType(G4Polymarker::squares);
  // Call this AddPrimitive to avoid re-doing sub-class code.
  G4OpenGLSceneHandler::AddPrimitive(oneSquare);
}

void G4OpenGLSceneHandler::AddPrimitive (const G4Scale& scale)
{
  G4VSceneHandler::AddPrimitive(scale);
}

//Method for handling G4Polyhedron objects for drawing solids.
void G4OpenGLSceneHandler::AddPrimitive (const G4Polyhedron& polyhedron) {

  // Assume all facets are planar convex quadrilaterals.
  // Draw each facet individually
  
  if (polyhedron.GetNoFacets() == 0) return;

  // Need access to data in G4OpenGLViewer.  static_cast doesn't work
  // with a virtual base class, so use dynamic_cast.
  G4OpenGLViewer* pGLViewer = dynamic_cast<G4OpenGLViewer*>(fpViewer);
  if (!pGLViewer) return;
  
  // Get view parameters that the user can force through the vis
  // attributes, thereby over-riding the current view parameter.
  G4ViewParameters::DrawingStyle drawing_style = GetDrawingStyle (fpVisAttribs);

  // Note that in stored mode, because this call gets embedded in a display
  //  list, it is the colour _at the time of_ creation of the display list, so
  //  even if the colour is changed, for example, by interaction with a Qt
  //  window, current_colour does not change.
  GLfloat* painting_colour;
  GLfloat current_colour [4];
  glGetFloatv (GL_CURRENT_COLOR, current_colour);
  
  G4bool isTransparent = false;
  if (current_colour[3] < 1.) {  // This object is transparent
    isTransparent = true;
  }

  
  if  (drawing_style == G4ViewParameters::hlr) {
    // This is the colour used to paint surfaces in hlr mode.
    GLfloat clear_colour[4];
    glGetFloatv (GL_COLOR_CLEAR_VALUE, clear_colour);
    painting_colour = clear_colour;
  } else {  // drawing_style == G4ViewParameters::hlhsr
    painting_colour = current_colour;
  }

  G4double lineWidth = GetLineWidth(fpVisAttribs);
  pGLViewer->ChangeLineWidth(lineWidth);

  G4bool isAuxEdgeVisible = GetAuxEdgeVisible (fpVisAttribs);

  G4bool clipping = pGLViewer->fVP.IsSection() || pGLViewer->fVP.IsCutaway();

  // Lighting disabled unless otherwise requested
#ifndef G4OPENGL_VERSION_2
  glDisable (GL_LIGHTING);
#endif

  switch (drawing_style) {
  case (G4ViewParameters::hlhsr):
    // Set up as for hidden line removal but paint polygon faces later...
  case (G4ViewParameters::hlr):
    glEnable (GL_STENCIL_TEST);
    // The stencil buffer is cleared in G4OpenGLViewer::ClearView.
    // The procedure below leaves it clear.
    glStencilFunc (GL_ALWAYS, 0, 1);
    glStencilOp (GL_INVERT, GL_INVERT, GL_INVERT);
    glEnable (GL_DEPTH_TEST);
    glDepthFunc (GL_LEQUAL);
    if (isTransparent) {
      // Transparent...
      glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
      glEnable(GL_COLOR_MATERIAL);
      glDisable (GL_CULL_FACE);
      glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
    } else {
      // Opaque...
      if (clipping) {
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
        glEnable(GL_COLOR_MATERIAL);
	glDisable (GL_CULL_FACE);
	glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
      } else {
        glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
        glEnable(GL_COLOR_MATERIAL);
	glEnable (GL_CULL_FACE);
	glCullFace (GL_BACK);
	glPolygonMode (GL_FRONT, GL_LINE);
      }
    }
    break;
  case (G4ViewParameters::hsr):
    glEnable (GL_DEPTH_TEST);
    glDepthFunc (GL_LEQUAL);    
    if (isTransparent) {
      // Transparent...
      glDepthMask (GL_FALSE);  // Make depth buffer read-only.
      glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
#ifndef G4OPENGL_VERSION_2
      glEnable(GL_COLOR_MATERIAL);
#endif
      glDisable (GL_CULL_FACE);
      glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    } else {
      // Opaque...
      glDepthMask (GL_TRUE);  // Make depth buffer writable (default).
      if (clipping) {
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
        glEnable(GL_COLOR_MATERIAL);
	glDisable (GL_CULL_FACE);
	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
      } else {
        glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
#ifndef G4OPENGL_VERSION_2
        glEnable(GL_COLOR_MATERIAL);
#endif
  glEnable (GL_CULL_FACE);
	glCullFace (GL_BACK);
	glPolygonMode (GL_FRONT, GL_FILL);
      }
    }
#ifndef G4OPENGL_VERSION_2
    if (!fProcessing2D) glEnable (GL_LIGHTING);
#endif
      break;
  case (G4ViewParameters::wireframe):
  default:
    glEnable (GL_DEPTH_TEST);
    glDepthFunc (GL_LEQUAL);    //??? was GL_ALWAYS
    glDisable (GL_CULL_FACE);
    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
    break;
  }

  //Loop through all the facets...
  fEdgeFlag = true;
#ifndef G4OPENGL_VERSION_2
  glBegin (GL_QUADS);
  glEdgeFlag (GL_TRUE);
#else
  fEmulate_GL_QUADS = true;
  glBeginVBO(GL_TRIANGLE_STRIP);
#endif
  G4bool notLastFace;
  do {

    //First, find vertices, edgeflags and normals and note "not last facet"...
    G4Point3D vertex[4];
    G4int edgeFlag[4];
    G4Normal3D normals[4];
    G4int nEdges;
    notLastFace = polyhedron.GetNextFacet(nEdges, vertex, edgeFlag, normals);

    //Loop through the four edges of each G4Facet...
    for(G4int edgeCount = 0; edgeCount < nEdges; ++edgeCount) {
      // Check to see if edge is visible or not...
      if (isAuxEdgeVisible) {
	edgeFlag[edgeCount] = 1;
      }
#ifndef G4OPENGL_VERSION_2
      if (edgeFlag[edgeCount] > 0) {
        if (fEdgeFlag != true) {
          glEdgeFlag (GL_TRUE);
          fEdgeFlag = true;
        }
      } else {
        if (fEdgeFlag != false) {
          glEdgeFlag (GL_FALSE);
          fEdgeFlag = false;
        }
      }
      glNormal3d (normals[edgeCount].x(),
		  normals[edgeCount].y(),
		  normals[edgeCount].z());
      glVertex3d (vertex[edgeCount].x(), 
		  vertex[edgeCount].y(),
		  vertex[edgeCount].z());
#else

      fOglVertex.push_back(vertex[edgeCount].x());
      fOglVertex.push_back(vertex[edgeCount].y());
      fOglVertex.push_back(vertex[edgeCount].z());
      
      fOglVertex.push_back(normals[edgeCount].x());
      fOglVertex.push_back(normals[edgeCount].y());
      fOglVertex.push_back(normals[edgeCount].z());

#endif
      
    }
   
    // HepPolyhedron produces triangles too; in that case add an extra
    // vertex identical to first...
    if (nEdges == 3) {
      G4int edgeCount = 3;
      normals[edgeCount] = normals[0];
      vertex[edgeCount] = vertex[0];
#ifndef G4OPENGL_VERSION_2
      edgeFlag[edgeCount] = -1;
      if (fEdgeFlag != false) {
        glEdgeFlag (GL_FALSE);
        fEdgeFlag = false;
      }

      glNormal3d (normals[edgeCount].x(),
		  normals[edgeCount].y(), 
		  normals[edgeCount].z());
      glVertex3d (vertex[edgeCount].x(),
		  vertex[edgeCount].y(), 
		  vertex[edgeCount].z());
#else
      fOglVertex.push_back(vertex[edgeCount].x());
      fOglVertex.push_back(vertex[edgeCount].y());
      fOglVertex.push_back(vertex[edgeCount].z());
      
      fOglVertex.push_back(normals[edgeCount].x());
      fOglVertex.push_back(normals[edgeCount].y());
      fOglVertex.push_back(normals[edgeCount].z());

#endif
    }
    // Trap situation where number of edges is > 4...
    if (nEdges > 4) {
      G4cerr <<
	"G4OpenGLSceneHandler::AddPrimitive(G4Polyhedron): WARNING"
	"\n   G4Polyhedron facet with " << nEdges << " edges" << G4endl;
    }

    
    // Do it all over again (twice) for hlr...
    if  (drawing_style == G4ViewParameters::hlr ||
	 drawing_style == G4ViewParameters::hlhsr) {

#ifndef G4OPENGL_VERSION_2
      glDisable(GL_COLOR_MATERIAL); // Revert to glMaterial for hlr/sr.
#endif

#ifndef G4OPENGL_VERSION_2
      glEnd ();  // Placed here to balance glBegin above, allowing GL
#else
      glEndVBO();
#endif
      // state changes below, then glBegin again.  Avoids
		 // having glBegin/End pairs *inside* loop in the more
		 // usual case of no hidden line removal.

      // Lighting disabled unless otherwise requested
      glDisable (GL_LIGHTING);

      // Draw through stencil...
      glStencilFunc (GL_EQUAL, 0, 1);
      glStencilOp (GL_KEEP, GL_KEEP, GL_KEEP);
      if (drawing_style == G4ViewParameters::hlhsr) {
	if (!fProcessing2D) glEnable (GL_LIGHTING);
      }
      glEnable (GL_DEPTH_TEST);
      glDepthFunc (GL_LEQUAL);    
      if (isTransparent) {
	// Transparent...
	glDepthMask (GL_FALSE);  // Make depth buffer read-only.
	glDisable (GL_CULL_FACE);
	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
      } else {
	// Opaque...
	glDepthMask (GL_TRUE);  // Make depth buffer writable (default).
	if (clipping) {
	  glDisable (GL_CULL_FACE);
	  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
	} else {
	  glEnable (GL_CULL_FACE);
	  glCullFace (GL_BACK);
	  glPolygonMode (GL_FRONT, GL_FILL);
	}
      }
      if  (drawing_style == G4ViewParameters::hlr) {
	if (isTransparent) {
	  // Transparent - don't paint...
	  goto end_of_drawing_through_stencil;
	}
      }
      if (isTransparent) {
	// Transparent...
	glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, painting_colour);
      } else {
	// Opaque...
	glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE, painting_colour);
      }
      glColor4fv (painting_colour);
#ifndef G4OPENGL_VERSION_2
      glBegin (GL_QUADS);
      glEdgeFlag (GL_TRUE);
      fEdgeFlag = true;
#else
      fEmulate_GL_QUADS = true;
      glBeginVBO(GL_TRIANGLE_STRIP);
#endif

      for (int edgeCount = 0; edgeCount < 4; ++edgeCount) {
#ifndef G4OPENGL_VERSION_2
        if (edgeFlag[edgeCount] > 0) {
          if (fEdgeFlag != true) {
            glEdgeFlag (GL_TRUE);
            fEdgeFlag = true;
          }
        } else {
          if (fEdgeFlag != false) {
            glEdgeFlag (GL_FALSE);
            fEdgeFlag = false;
          }
        }
				glNormal3d (normals[edgeCount].x(),
							normals[edgeCount].y(),
							normals[edgeCount].z());
        glVertex3d (vertex[edgeCount].x(),
                    vertex[edgeCount].y(),
                    vertex[edgeCount].z());
#else
				fOglVertex.push_back(vertex[edgeCount].x());
        fOglVertex.push_back(vertex[edgeCount].y());
        fOglVertex.push_back(vertex[edgeCount].z());
        
        fOglVertex.push_back(normals[edgeCount].x());
        fOglVertex.push_back(normals[edgeCount].y());
        fOglVertex.push_back(normals[edgeCount].z());

#endif
      }
#ifndef G4OPENGL_VERSION_2
      glEnd ();
#else
      glEndVBO();
#endif
    end_of_drawing_through_stencil:

      // and once more to reset the stencil bits...
      glStencilFunc (GL_ALWAYS, 0, 1);
      glStencilOp (GL_INVERT, GL_INVERT, GL_INVERT);
      glDepthFunc (GL_LEQUAL);  // to make sure line gets drawn.  
      if (isTransparent) {
	// Transparent...
	glDisable (GL_CULL_FACE);
	glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
      } else {
	// Opaque...
	if (clipping) {
	  glDisable (GL_CULL_FACE);
	  glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
	} else {
	  glEnable (GL_CULL_FACE);
	  glCullFace (GL_BACK);
	  glPolygonMode (GL_FRONT, GL_LINE);
	}
      }
      glDisable (GL_LIGHTING);
      glColor4fv (current_colour);
      fEdgeFlag = true;
#ifndef G4OPENGL_VERSION_2
      glBegin (GL_QUADS);
      glEdgeFlag (GL_TRUE);
      fEdgeFlag = true;
#else
      fEmulate_GL_QUADS = true;
      glBeginVBO(GL_TRIANGLE_STRIP);
#endif
      for (int edgeCount = 0; edgeCount < 4; ++edgeCount) {
#ifndef G4OPENGL_VERSION_2
        if (edgeFlag[edgeCount] > 0) {
          if (fEdgeFlag != true) {
            glEdgeFlag (GL_TRUE);
            fEdgeFlag = true;
          }
        } else {
          if (fEdgeFlag != false) {
            glEdgeFlag (GL_FALSE);
            fEdgeFlag = false;
          }
        }
        glNormal3d (normals[edgeCount].x(),
                    normals[edgeCount].y(),
                    normals[edgeCount].z());
        glVertex3d (vertex[edgeCount].x(),
                    vertex[edgeCount].y(),
                    vertex[edgeCount].z());
#else
				fOglVertex.push_back(vertex[edgeCount].x());
        fOglVertex.push_back(vertex[edgeCount].y());
        fOglVertex.push_back(vertex[edgeCount].z());
        
        fOglVertex.push_back(normals[edgeCount].x());
        fOglVertex.push_back(normals[edgeCount].y());
        fOglVertex.push_back(normals[edgeCount].z());

#endif
      }
#ifndef G4OPENGL_VERSION_2
      glEnd ();
#else
      glEndVBO();
#endif

      glDepthFunc (GL_LEQUAL);   // Revert for next facet.
      fEdgeFlag = true;
#ifndef G4OPENGL_VERSION_2
      glBegin (GL_QUADS);      // Ready for next facet.  GL
      glEdgeFlag (GL_TRUE);
      fEdgeFlag = true;
      // says it ignores incomplete
      // quadrilaterals, so final empty
      // glBegin/End sequence should be OK.
#else
      fEmulate_GL_QUADS = true;
      glBeginVBO(GL_TRIANGLE_STRIP);
#endif
    }
  } while (notLastFace);  
  
#ifndef G4OPENGL_VERSION_2
  glEnd ();
#else
	
// FIXME: du grand n'importe quoi en test
// Cube optimization
  
  // store old DrawType because in case of optimization it could be changed
  GLenum oldDrawArrayType = fDrawArrayType;

  if (dynamic_cast<const G4PolyhedronTrd2*>(&polyhedron)) {
//    OptimizeVBOForTrd();
  } else if (dynamic_cast<const G4PolyhedronCons*>(&polyhedron)) {
//    OptimizeVBOForCons((polyhedron.GetNoVertices()-2)/2 ); // top + bottom + all faces
  }

  glEndVBO();
  fDrawArrayType = oldDrawArrayType;
#endif
  
  glDisable (GL_STENCIL_TEST);  // Revert to default for next primitive.
  glDepthMask (GL_TRUE);        // Revert to default for next primitive.
#ifndef G4OPENGL_VERSION_2
  glDisable (GL_LIGHTING);      // Revert to default for next primitive.
#endif
}

void G4OpenGLSceneHandler::AddCompound(const G4VTrajectory& traj) {
  G4VSceneHandler::AddCompound(traj);  // For now.
}

void G4OpenGLSceneHandler::AddCompound(const G4VHit& hit) {
  G4VSceneHandler::AddCompound(hit);  // For now.
}

void G4OpenGLSceneHandler::AddCompound(const G4VDigi& digi) {
  G4VSceneHandler::AddCompound(digi);  // For now.
}

void G4OpenGLSceneHandler::AddCompound(const G4THitsMap<G4double>& hits) {
  G4VSceneHandler::AddCompound(hits);  // For now.
}

void G4OpenGLSceneHandler::AddCompound(const G4THitsMap<G4StatDouble>& hits) {
  G4VSceneHandler::AddCompound(hits);  // For now.
}


#ifdef G4OPENGL_VERSION_2

// Optimize vertex and indices in order to render less vertex in OpenGL VBO/IBO
void G4OpenGLSceneHandler::OptimizeVBOForTrd(){
  
  /* HOW IT IS BUILD (as we receive it from fOglVertex : 
   */

   std::vector<double> vertices;
  vertices.insert (vertices.end(),fOglVertex.begin(),fOglVertex.begin()+6*6); // ABCDEF
  vertices.insert (vertices.end(),fOglVertex.begin()+9*6,fOglVertex.begin()+9*6+6); // G
  vertices.insert (vertices.end(),fOglVertex.begin()+13*6,fOglVertex.begin()+13*6+6); // H
  fOglVertex = vertices;
  
  int myarray [] = {
    3,2,0,1,4,5,7,6, 6,0,4,3,7,2,6,1,5
  };
  fOglIndices.insert(fOglIndices.begin(), myarray, myarray+17/*36*/);

  fDrawArrayType = GL_TRIANGLE_STRIP;
}

// Optimize vertex and indices in order to render less vertex in OpenGL VBO/IBO
void G4OpenGLSceneHandler::OptimizeVBOForCons(G4int aNoFaces){
  // Optimized, 1st level : 10f/15sec with 1000 cones
  //    DrawElements:208 vertex and 605 (2*100+2*100+2*100+5) indices for a 100 face cone

  /* surface of polycone : could be optimized
   for 100 faces : 
   - 100*4 = 400 points
   - 100*2+2 = 202 points with TRIANGLE_STRIP
   Total :
     n*4+n*4+n*4 = n*12
    optimize : n*2+2+1+n+1 = n*3+3 (factor 4)
    but could do better : n faces should give = n*2+2
   */
  
  /*
         0
        / \
       2---4   6 ....2
       |   |
       3---5   7 ....3
        \ /
         1
   */
  // First, faces
  std::vector<double> vertices;

  // Add bottom and top vertex
  // aNoFaces*4*6+6 : nb Faces * 4 points per face * 6 vertex by point + 1 point offset
  vertices.insert (vertices.end(),fOglVertex.begin()+ (aNoFaces*4)*6,fOglVertex.begin()+(aNoFaces*4)*6+6); // 0
  vertices.insert (vertices.end(),fOglVertex.begin()+ (aNoFaces*8+1)*6,fOglVertex.begin()+(aNoFaces*8+1)*6+6); // 1
  
  // Add facets points
  G4int posInVertice;
  for (G4int a = 0; a<aNoFaces; a++) {
    posInVertice = a*4*6;
    vertices.insert (vertices.end(),fOglVertex.begin()+posInVertice,fOglVertex.begin()+posInVertice+1*6+6); // AB
  }
  vertices.insert (vertices.end(),fOglVertex.begin(),fOglVertex.begin()+1*6*6); // AB
  fOglVertex = vertices;

  // Add indices for top :
  // simple version :    0-2-0-4-0-6-0-8-0-10..
  // optimized version : 2-0-4-6- 6-0-8-10.. but we have to deal with odd faces numbers
  for (G4int a=0; a<aNoFaces; a++) {
    fOglIndices.push_back(0);
    fOglIndices.push_back(a*2+2);
  }
  // close strip
  fOglIndices.push_back(0);
  fOglIndices.push_back(2);

  // Add indices for faces
  for (G4int a = 0; a<aNoFaces; a++) {
    fOglIndices.push_back(a*2+2);
    fOglIndices.push_back(a*2+1+2);
  }
  fOglIndices.push_back(2);
  fOglIndices.push_back(2+1);
  
  // Second : top
  // 3-1-5-1-7-1-9-1..
  for (G4int a=0; a<aNoFaces; a++) {
    fOglIndices.push_back(a*2+3);
    fOglIndices.push_back(1);
  }
  // close strip
  fOglIndices.push_back(0+3);
  
  fDrawArrayType = GL_TRIANGLE_STRIP;
  fEmulate_GL_QUADS = false;
}

void G4OpenGLSceneHandler::glBeginVBO(GLenum type)  {
  fDrawArrayType = type;
#ifndef G4VIS_BUILD_OPENGLWT_DRIVER
  glGenBuffers(1,&fVertexBufferObject);
  glGenBuffers(1,&fIndicesBufferObject);
#else
  fVertexBufferObject = glCreateBuffer(); //glGenBuffer(1,fVertexBufferObject_2)
  fIndicesBufferObject = glCreateBuffer(); //glGenBuffer(1,fIndicesBufferObject_2)
#endif

  // clear data and indices for OpenGL
  fOglVertex.clear();
  fOglIndices.clear();
}

// 2 cases :
/*
 glDrawArray : if there is no vertex indices : fOglIndices.size() == 0
 glDrawElements : if there is vertex indices : fOglIndices.size() != 0
 
 */
void G4OpenGLSceneHandler::glEndVBO()  {
  if (fOglIndices.size() == 0) {


    std::vector<double> vertices;
    // check if it is a GL_QUADS emulation
    if (fEmulate_GL_QUADS == true) {
      fEmulate_GL_QUADS = false;
      // A point has 6 double : Vx Vy Vz Nx Ny Nz
      // A QUAD should be like this
      /*
           0   3/4  7/8   ..
       
           1   2/5  6/9   ..
       */
      // And if 3==4 and 2==5, we should do it like this for a TRIANGLES_STRIP
      /*
       0   4   8   ..
       | / | / |
       1   5   9   ..
       // Optimized, 1st level : 24f/15sec with 10 cones
       // non Optimized, 1st level : 12f/15sec with 10 cones
       */
      // should be 4 points
      for (unsigned int a=0; a<fOglVertex.size(); a+=6*4) {
        vertices.insert (vertices.end(),fOglVertex.begin()+a,fOglVertex.begin()+a+1*6+6); // 0-1
        // if 2-3 == 4-5, do not add them
        // if differents, we are obliged to create a new GL_TRIANGLE_STRIP
        if (a+4*6+5 < fOglVertex.size()) {
          if ((fOglVertex[a+2*6+0] != fOglVertex[a+5*6+0]) || //Vx for 2 and 5
              (fOglVertex[a+2*6+1] != fOglVertex[a+5*6+1]) || //Vy for 2 and 5
              (fOglVertex[a+2*6+2] != fOglVertex[a+5*6+2]) || //Vz for 2 and 5
              (fOglVertex[a+2*6+3] != fOglVertex[a+5*6+3]) || //Px for 2 and 5
              (fOglVertex[a+2*6+4] != fOglVertex[a+5*6+4]) || //Py for 2 and 5
              (fOglVertex[a+2*6+5] != fOglVertex[a+5*6+5]) || //Pz for 2 and 5
              
              (fOglVertex[a+3*6+0] != fOglVertex[a+4*6+0]) || //Vx for 3 and 4
              (fOglVertex[a+3*6+1] != fOglVertex[a+4*6+1]) || //Vy for 3 and 4
              (fOglVertex[a+3*6+2] != fOglVertex[a+4*6+2]) || //Vz for 3 and 4
              (fOglVertex[a+3*6+3] != fOglVertex[a+4*6+3]) || //Px for 3 and 4
              (fOglVertex[a+3*6+4] != fOglVertex[a+4*6+4]) || //Py for 3 and 4
              (fOglVertex[a+3*6+5] != fOglVertex[a+4*6+5])) { //Pz for 3 and 4
            // add last points
            vertices.insert (vertices.end(),fOglVertex.begin()+a+3*6,fOglVertex.begin()+a+3*6+6); // 3
            vertices.insert (vertices.end(),fOglVertex.begin()+a+2*6,fOglVertex.begin()+a+2*6+6); // 2
            // build and send the GL_TRIANGLE_STRIP
            drawVBOArray(vertices);
            vertices.clear();
          }
        } else { // end of volume
          vertices.insert (vertices.end(),fOglVertex.begin()+a+3*6,fOglVertex.begin()+a+3*6+6); // 3
          vertices.insert (vertices.end(),fOglVertex.begin()+a+2*6,fOglVertex.begin()+a+2*6+6); // 2
        }
      }
      fOglVertex = vertices;
    }

    drawVBOArray(fOglVertex);

  } else {
  
    // Bind VBO
    glBindBuffer(GL_ARRAY_BUFFER, fVertexBufferObject);
    
    // Load fOglVertex into VBO
#ifndef G4VIS_BUILD_OPENGLWT_DRIVER
    int sizeV = fOglVertex.size();
    // FIXME : perhaps a problem withBufferData in OpenGL other than WebGL ?
//    void glBufferData(	GLenum target, GLsizeiptr size, const GLvoid * data, GLenum usage);
    glBufferData(GL_ARRAY_BUFFER, sizeof(double)*sizeV, &fOglVertex[0], GL_STATIC_DRAW);
#else
    glBufferDatafv(GL_ARRAY_BUFFER, fOglVertex.begin(), fOglVertex.end(), GL_STATIC_DRAW);
#endif
    
    // Bind IBO
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, fIndicesBufferObject);
    
    // Load fOglVertex into VBO
#ifndef G4VIS_BUILD_OPENGLWT_DRIVER
    int sizeI = fOglIndices.size();
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(int)*sizeI, &fOglIndices[0], GL_STATIC_DRAW);
#else
    glBufferDataiv(GL_ELEMENT_ARRAY_BUFFER, fOglIndices.begin(), fOglIndices.end(), GL_STATIC_DRAW, GL_UNSIGNED_BYTE);
#endif
    
    //----------------------------
    // Draw VBO
    //----------------------------
    glBindBuffer(GL_ARRAY_BUFFER, fVertexBufferObject);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, fIndicesBufferObject);
        
    // the fVertexPositionAttribute is inside the G4OpenGLViewer
    G4OpenGLViewer* pGLViewer = dynamic_cast<G4OpenGLViewer*>(fpViewer);
    if (pGLViewer) {
      glEnableVertexAttribArray(pGLViewer->fVertexPositionAttribute);

      glVertexAttribPointer(pGLViewer->fVertexPositionAttribute,
                            3,     // size: Every vertex has an X, Y anc Z component
                            GL_FLOAT, // type: They are floats
                            GL_FALSE, // normalized: Please, do NOT normalize the vertices
                            2*3*4, // stride: The first byte of the next vertex is located this
                            //         amount of bytes further. The format of the VBO is
                            //         vx, vy, vz, nx, ny, nz and every element is a
                            //         Float32, hence 4 bytes large
                            0);    // offset: The byte position of the first vertex in the buffer
    }
    
    
    glBindBuffer(GL_ARRAY_BUFFER, fVertexBufferObject);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, fIndicesBufferObject);
//    glDrawElements(fDrawArrayType, fOglIndices.size(), GL_UNSIGNED_SHORT, 0);
    glDrawElements(fDrawArrayType, fOglIndices.size(), GL_UNSIGNED_SHORT, 0);
    
    if (pGLViewer) {
      glDisableVertexAttribArray(pGLViewer->fVertexPositionAttribute);
    }

    // delete the buffer
#ifndef G4VIS_BUILD_OPENGLWT_DRIVER
    glDeleteBuffers(1,&fVertexBufferObject);
#else
    glDeleteBuffer(fVertexBufferObject);
#endif
  }
}
          
void G4OpenGLSceneHandler::drawVBOArray(std::vector<double> vertices)  {
#ifndef G4VIS_BUILD_OPENGLWT_DRIVER
  glGenBuffers(1,&fVertexBufferObject);
  glGenBuffers(1,&fIndicesBufferObject);
#else
  fVertexBufferObject = glCreateBuffer(); //glGenBuffer(1,fVertexBufferObject_2)
  fIndicesBufferObject = glCreateBuffer(); //glGenBuffer(1,fIndicesBufferObject_2)
#endif

  // Bind this buffer
  glBindBuffer(GL_ARRAY_BUFFER, fVertexBufferObject);
  // Load oglData into VBO
#ifndef G4VIS_BUILD_OPENGLWT_DRIVER
  int s = vertices.size();
  glBufferData(GL_ARRAY_BUFFER, sizeof(double)*s, &vertices[0], GL_STATIC_DRAW);
#else
  glBufferDatafv(GL_ARRAY_BUFFER, vertices.begin(), vertices.end(), GL_STATIC_DRAW);
#endif
  
  //----------------------------
  // Draw VBO
  //----------------------------
  glBindBuffer(GL_ARRAY_BUFFER, fVertexBufferObject);
  
  // the fVertexPositionAttribute is inside the G4OpenGLViewer
  G4OpenGLViewer* pGLViewer = dynamic_cast<G4OpenGLViewer*>(fpViewer);
  if (pGLViewer) {
    glEnableVertexAttribArray(pGLViewer->fVertexPositionAttribute);

//    glVertexAttribPointer(GLuint index, GLint size, GLenum type, GLboolean normalized, GLsizei stride, const GLvoid *pointer)

/*
 GL_DOUBLE
 Warning: This section describes legacy OpenGL APIs that have been removed from core OpenGL 3.1 and above (they are only deprecated in OpenGL 3.0). It is recommended that you not use this functionality in your programs.
 
 glLoadMatrixd, glRotated and any other function that have to do with the double type. Most GPUs don't support GL_DOUBLE (double) so the driver will convert the data to GL_FLOAT (float) and send to the GPU. If you put GL_DOUBLE data in a VBO, the performance might even be much worst than immediate mode (immediate mode means glBegin, glVertex, glEnd). GL doesn't offer any better way to know what the GPU prefers.
 */
#ifndef G4VIS_BUILD_OPENGLWT_DRIVER
    glVertexAttribPointer(pGLViewer->fVertexPositionAttribute,
                          3,     // size: Every vertex has an X, Y anc Z component
                          GL_DOUBLE, // type: They are double
                          GL_FALSE, // normalized: Please, do NOT normalize the vertices
                          6*sizeof(double), // stride: The first byte of the next vertex is located this
                          //         amount of bytes further. The format of the VBO is
                          //         vx, vy, vz, nx, ny, nz and every element is a
                          //         Float32, hence 4 bytes large
                          0);    // offset: The byte position of the first vertex in the buffer
#else
    glVertexAttribPointer(pGLViewer->fVertexPositionAttribute,
                          3,     // size: Every vertex has an X, Y anc Z component
                          GL_FLOAT, // type: They are floats
                          GL_FALSE, // normalized: Please, do NOT normalize the vertices
                          2*3*4,    // stride: The first byte of the next vertex is located this
                          //         amount of bytes further. The format of the VBO is
                          //         vx, vy, vz, nx, ny, nz and every element is a
                          //         Float32, hence 4 bytes large
                          0);    // offset: The byte position of the first vertex in the buffer
#endif
  }
  
  glDrawArrays(fDrawArrayType, // GL_POINTS, GL_LINE_STRIP, GL_LINE_LOOP, GL_LINES, GL_TRIANGLE_FAN, GL_TRIANGLE_STRIP, and GL_TRIANGLES
               0, vertices.size()/6);
  if (pGLViewer) {
#ifndef G4VIS_BUILD_OPENGLWT_DRIVER
    glDisableClientState( GL_VERTEX_ARRAY );
#else
    glDisableVertexAttribArray(pGLViewer->fVertexPositionAttribute);
#endif
  }
  
  // delete the buffer
#ifndef G4VIS_BUILD_OPENGLWT_DRIVER
  glDeleteBuffers(1,&fVertexBufferObject);
#else
  glDeleteBuffer(fVertexBufferObject);
#endif
}
#endif

#endif
