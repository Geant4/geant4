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
// Andrew Walkden  27th March 1996
// OpenGL stored scene - creates OpenGL display lists.
// OpenGL immediate scene - draws immediately to buffer
//                           (saving space on server).


#  include "G4OpenGLSceneHandler.hh"
#  include "G4OpenGLViewer.hh"
#  include "G4OpenGLTransform3D.hh"
#  include "G4Point3D.hh"
#  include "G4Normal3D.hh"
#  include "G4Transform3D.hh"
#  include "G4Polyline.hh"
#  include "G4Polymarker.hh"
#  include "G4Text.hh"
#  include "G4Circle.hh"
#  include "G4Square.hh"
#  include "G4VMarker.hh"
#  include "G4Polyhedron.hh"
#  include "G4VisAttributes.hh"
#  include "G4PhysicalVolumeModel.hh"
#  include "G4VPhysicalVolume.hh"
#  include "G4LogicalVolume.hh"
#  include "G4VSolid.hh"
#  include "G4Scene.hh"
#  include "G4VisExtent.hh"
#  include "G4AttHolder.hh"
#  include "G4PhysicalConstants.hh"
#  include "G4RunManager.hh"
#  include "G4Run.hh"
#  include "G4RunManagerFactory.hh"
#  include "G4Mesh.hh"
#  include "G4PseudoScene.hh"
#  include "G4VisManager.hh"

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
      // No scene - shouldn't happen
      glFlush();
      return;
    }
    // Get event from modeling parameters
    if (!fpModel) {
      // No model - shouldn't happen
      glFlush();
      return;
    }
    const G4ModelingParameters* modelingParameters =
    fpModel->GetModelingParameters();
    if (!modelingParameters) {
      // No modeling parameters - shouldn't happen
      glFlush();
      return;
    }
    const G4Event* thisEvent = modelingParameters->GetEvent();
    if (!thisEvent) {
      // No event, so not in event loop.
      if (fFlushAction == endOfEvent) {
        fFlushAction = endOfRun;
      } else if (fFlushAction == NthEvent) {
        fFlushAction = NthPrimitive;
      }
    }
    G4RunManager* runMan = G4RunManagerFactory::GetMasterRunManager();
    if (!runMan) {
      // No run manager - shouldn't happen
      glFlush();
      return;
    }
    const G4Run* thisRun = runMan->GetCurrentRun();
    if (!thisRun) {
      // No run, so not in event loop.
      if (fFlushAction == endOfRun) {
        fFlushAction = NthPrimitive;
      } else if (fFlushAction == NthEvent) {
        fFlushAction = NthPrimitive;
      }
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
        fEntitiesFlushInterval = 1;
	[[fallthrough]];  // Fall through to NthPrimitive.
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
        fEntitiesFlushInterval = 1;
	[[fallthrough]];  // Fall through to NthPrimitive.
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

G4DisplacedSolid* G4OpenGLSceneHandler::CreateSectionSolid ()
{
  return G4VSceneHandler::CreateSectionSolid();
  // If clipping done in G4OpenGLViewer::SetView
  // return 0;
  // Note: if you change this, you must also change
  // G4OpenGLStoredViewer::CompareForKernelVisit
}

G4DisplacedSolid* G4OpenGLSceneHandler::CreateCutawaySolid ()
{
  // return G4VSceneHandler::CreateCutawaySolid();
  // If cutaway done in G4OpenGLViewer::SetView.
  return 0;
  // Note: if you change this, you must also change
  // G4OpenGLStoredViewer::CompareForKernelVisit
}

void G4OpenGLSceneHandler::AddPrimitive (const G4Polyline& line)
{
  std::size_t nPoints = line.size ();
  if (nPoints <= 0) return;

  // Note: colour and depth test treated in sub-class.

  glDisable (GL_LIGHTING);
  
  G4double lineWidth = GetLineWidth(fpVisAttribs);
  // Need access to method in G4OpenGLViewer.  static_cast doesn't
  // work with a virtual base class, so use dynamic_cast.  No need to
  // test the outcome since viewer is guaranteed to be a
  // G4OpenGLViewer, but test it anyway to keep Coverity happy.
  G4OpenGLViewer* pGLViewer = dynamic_cast<G4OpenGLViewer*>(fpViewer);
  if (pGLViewer) pGLViewer->ChangeLineWidth(lineWidth);

  fEdgeFlag = true;
  glBegin (GL_LINE_STRIP);
  // No ned glEdgeFlag for lines :
  // Boundary and nonboundary edge flags on vertices are significant only if GL_POLYGON_MODE is set to GL_POINT or GL_LINE.  See glPolygonMode.
  
  //  glEdgeFlag (GL_TRUE);
  for (std::size_t iPoint = 0; iPoint < nPoints; ++iPoint) {
  G4double x, y, z;
    x = line[iPoint].x(); 
    y = line[iPoint].y();
    z = line[iPoint].z();
    glVertex3d (x, y, z);
  }
  glEnd ();
}

void G4OpenGLSceneHandler::AddPrimitive (const G4Polymarker& polymarker)
{
  if (polymarker.size() == 0) {
    return;
  }

  // Note: colour and depth test treated in sub-class.

  glDisable (GL_LIGHTING);
  
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
        if (!hashedWarned) {
          G4cout << "Hashed fill style in G4OpenGLSceneHandler."
          << "\n  Not implemented.  Using G4VMarker::filled."
          << G4endl;
          hashedWarned = true;
        }
        // Maybe use
        //glPolygonStipple (fStippleMaskHashed);
	[[fallthrough]];   // Drop through to filled...
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
        size = 1.;
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
      glBegin (GL_POLYGON);
      for (i = 0, phi = startPhi; i < nSides; i++, phi += dPhi) {
	G4Vector3D r = start; r.rotate(phi, viewpointDirection);
	G4Vector3D p = polymarker[iPoint] + r;
	glVertex3d (p.x(), p.y(), p.z());
      }
      glEnd ();
    }

  } else { // Size specified in screen (window) coordinates.

    pGLViewer->ChangePointSize(size);

    //Antialiasing only for circles
    switch (polymarker.GetMarkerType()) {
    default:
    case G4Polymarker::dots:
    case G4Polymarker::circles:
      glEnable (GL_POINT_SMOOTH); break;
    case G4Polymarker::squares:
      glDisable (GL_POINT_SMOOTH); break;
    }
    glBegin (GL_POINTS);
    for (size_t iPoint = 0; iPoint < polymarker.size (); iPoint++) {
      G4Point3D centre = polymarker[iPoint];
      glVertex3d(centre.x(),centre.y(),centre.z());
    }
    glEnd();     
  }
}

void G4OpenGLSceneHandler::AddPrimitive (const G4Text& text) {
  // Pass to specific viewer via virtual function DrawText.
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
  GLfloat clear_colour[4];
  GLfloat current_colour[4];
  glGetFloatv (GL_CURRENT_COLOR, current_colour);
  
  G4bool isTransparent = false;
  if (current_colour[3] < 1.) {  // This object is transparent
    isTransparent = true;
  }

  if  (drawing_style == G4ViewParameters::hlr) {
    // This is the colour used to paint surfaces in hlr mode.
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
  glDisable (GL_LIGHTING);

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
      //glDisable (GL_CULL_FACE);
      glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
    } else {
      // Opaque...
      if (clipping) {
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
        glEnable(GL_COLOR_MATERIAL);
	//glDisable (GL_CULL_FACE);
	glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
      } else {
        glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
        glEnable(GL_COLOR_MATERIAL);
	//glEnable (GL_CULL_FACE);
	//glCullFace (GL_BACK);
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
      glEnable(GL_COLOR_MATERIAL);
      //glDisable (GL_CULL_FACE);
      glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    } else {
      // Opaque...
      glDepthMask (GL_TRUE);  // Make depth buffer writable (default).
      if (clipping) {
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
        glEnable(GL_COLOR_MATERIAL);
	//glDisable (GL_CULL_FACE);
	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
      } else {
        glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
        glEnable(GL_COLOR_MATERIAL);
        //glEnable (GL_CULL_FACE);
	//glCullFace (GL_BACK);
	glPolygonMode (GL_FRONT, GL_FILL);
      }
    }
    if (!fProcessing2D) glEnable (GL_LIGHTING);
      break;
  case (G4ViewParameters::wireframe):
  default:
    glEnable (GL_DEPTH_TEST);
    glDepthFunc (GL_LEQUAL);    //??? was GL_ALWAYS
    //glDisable (GL_CULL_FACE);
    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
    break;
  }

  //Loop through all the facets...
  fEdgeFlag = true;
  glBegin (GL_QUADS);
  glEdgeFlag (GL_TRUE);
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
    }
   
    // HepPolyhedron produces triangles too; in that case add an extra
    // vertex identical to first...
    if (nEdges == 3) {
      G4int edgeCount = 3;
      normals[edgeCount] = normals[0];
      vertex[edgeCount] = vertex[0];
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

      glDisable(GL_COLOR_MATERIAL); // Revert to glMaterial for hlr/sr.
      glEnd ();  // Placed here to balance glBegin above, allowing GL
      
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
	//glDisable (GL_CULL_FACE);
	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
      } else {
	// Opaque...
	glDepthMask (GL_TRUE);  // Make depth buffer writable (default).
	if (clipping) {
	  //glDisable (GL_CULL_FACE);
	  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
	} else {
	  //glEnable (GL_CULL_FACE);
	  //glCullFace (GL_BACK);
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
      glBegin (GL_QUADS);
      glEdgeFlag (GL_TRUE);
      fEdgeFlag = true;

      for (int edgeCount = 0; edgeCount < 4; ++edgeCount) {
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
      }
      glEnd ();
    end_of_drawing_through_stencil:

      // and once more to reset the stencil bits...
      glStencilFunc (GL_ALWAYS, 0, 1);
      glStencilOp (GL_INVERT, GL_INVERT, GL_INVERT);
      glDepthFunc (GL_LEQUAL);  // to make sure line gets drawn.  
      if (isTransparent) {
	// Transparent...
	//glDisable (GL_CULL_FACE);
	glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
      } else {
	// Opaque...
	if (clipping) {
	  //glDisable (GL_CULL_FACE);
	  glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
	} else {
	  //glEnable (GL_CULL_FACE);
	  //glCullFace (GL_BACK);
	  glPolygonMode (GL_FRONT, GL_LINE);
	}
      }
      glDisable (GL_LIGHTING);
      glColor4fv (current_colour);
      fEdgeFlag = true;
      glBegin (GL_QUADS);
      glEdgeFlag (GL_TRUE);
      fEdgeFlag = true;
      for (int edgeCount = 0; edgeCount < 4; ++edgeCount) {
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
      }
      glEnd ();

      glDepthFunc (GL_LEQUAL);   // Revert for next facet.
      fEdgeFlag = true;
      glBegin (GL_QUADS);      // Ready for next facet.  GL
      glEdgeFlag (GL_TRUE);
      fEdgeFlag = true;
      // says it ignores incomplete
      // quadrilaterals, so final empty
      // glBegin/End sequence should be OK.
    }
  } while (notLastFace);  
  
  glEnd ();
  glDisable (GL_STENCIL_TEST);  // Revert to default for next primitive.
  glDepthMask (GL_TRUE);        // Revert to default for next primitive.
  glDisable (GL_LIGHTING);      // Revert to default for next primitive.
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

void G4OpenGLSceneHandler::AddCompound(const G4Mesh& mesh) {
  StandardSpecialMeshRendering(mesh);
}
