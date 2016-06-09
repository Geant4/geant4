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
// $Id$
//
// 
// Andrew Walkden  27th March 1996
// OpenGL stored scene - creates OpenGL display lists.
// OpenGL immediate scene - draws immediately to buffer
//                           (saving space on server).

#ifdef G4VIS_BUILD_OPENGL_DRIVER

// Included here - problems with HP compiler if not before other includes?
#include "G4NURBS.hh"

// Here follows a special for Mesa, the OpenGL emulator.  Does not affect
// other OpenGL's, as far as I'm aware.   John Allison 18/9/96.
#define CENTERLINE_CLPP  /* CenterLine C++ workaround: */
// Also seems to be required for HP's CC and AIX xlC, at least.

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

G4OpenGLSceneHandler::G4OpenGLSceneHandler (G4VGraphicsSystem& system,
                                            G4int id,
                                            const G4String& name):
G4VSceneHandler (system, id, name),
fPickName(0),
// glFlush take about 90% time.  Dividing glFlush number by 100 will
// change the first vis time from 100% to 10+90/100 = 10,9%.
fEventsDrawInterval(1),
fEventsWaitingToBeFlushed(0),
fThreePassCapable(false),
fSecondPassForTransparencyRequested(false),
fSecondPassForTransparency(false),
fThirdPassForNonHiddenMarkersRequested(false),
fThirdPassForNonHiddenMarkers(false)
{}

G4OpenGLSceneHandler::~G4OpenGLSceneHandler ()
{
  ClearStore ();
}

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

void G4OpenGLSceneHandler::ClearAndDestroyAtts()
{
  std::map<GLuint, G4AttHolder*>::iterator i;
  for (i = fPickMap.begin(); i != fPickMap.end(); ++i) delete i->second;
  fPickMap.clear();
}

void G4OpenGLSceneHandler::ScaledFlush()
{
  fEventsWaitingToBeFlushed++;
  if (fEventsWaitingToBeFlushed < fEventsDrawInterval) return;
  glFlush();
  fEventsWaitingToBeFlushed = 0;
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
  // Clipping done in G4OpenGLViewer::SetView.
  // return 0;

  // But...OpenGL no longer seems to reconstruct clipped edges, so,
  // when the BooleanProcessor is up to it, abandon this and use
  // generic clipping in G4VSceneHandler::CreateSectionSolid...
  return G4VSceneHandler::CreateSectionSolid();
}

G4VSolid* G4OpenGLSceneHandler::CreateCutawaySolid ()
{
  // Cutaway done in G4OpenGLViewer::SetView.
  // return 0;

  // But...if not, when the BooleanProcessor is up to it...
  return G4VSceneHandler::CreateCutawaySolid();
}

void G4OpenGLSceneHandler::AddPrimitive (const G4Polyline& line)
{
  G4int nPoints = line.size ();
  if (nPoints <= 0) return;

  // Note: colour and depth test treated in sub-class.

  glDisable (GL_LIGHTING);

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA =
    fpViewer -> GetApplicableVisAttributes (line.GetVisAttributes ());

  G4double lineWidth = GetLineWidth(pVA);
  // Need access to method in G4OpenGLViewer.  static_cast doesn't
  // work with a virtual base class, so use dynamic_cast.  No need to
  // test the outcome since viewer is guaranteed to be a
  // G4OpenGLViewer, but test it anyway to keep Coverity happy.
  G4OpenGLViewer* pGLViewer = dynamic_cast<G4OpenGLViewer*>(fpViewer);
  if (pGLViewer) pGLViewer->ChangeLineWidth(lineWidth);

  glBegin (GL_LINE_STRIP);
  for (G4int iPoint = 0; iPoint < nPoints; iPoint++) {
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
  
  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA =
    fpViewer -> GetApplicableVisAttributes (polymarker.GetVisAttributes ());

  G4double lineWidth = GetLineWidth(pVA);
  // Need access to method in G4OpenGLViewer.  static_cast doesn't
  // work with a virtual base class, so use dynamic_cast.  No need to
  // test the outcome since viewer is guaranteed to be a
  // G4OpenGLViewer, but test it anyway to keep Coverity happy.
  G4OpenGLViewer* pGLViewer = dynamic_cast<G4OpenGLViewer*>(fpViewer);
  if (pGLViewer) pGLViewer->ChangeLineWidth(lineWidth);

  G4VMarker::FillStyle style = polymarker.GetFillStyle();

  // G4bool filled = false;  Not actually used - comment out to prevent compiler warnings (JA).
  static G4bool hashedWarned = false;
  
  switch (style) {
  case G4VMarker::noFill: 
    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
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
    // Drop through to filled...  
  case G4VMarker::filled:
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    //filled = true;
    break;
  }

  MarkerSizeType sizeType;
  G4double size = GetMarkerSize(polymarker, sizeType);

  // Draw...
  if (sizeType == world) {  // Size specified in world coordinates.

    G4int nSides;
    G4double startPhi;
    switch (polymarker.GetMarkerType()) {
    default:
    case G4Polymarker::dots:
      size = 1.;
      // Drop through to circles
    case G4Polymarker::circles:
      nSides = GetNoOfSides(pVA);
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
  
  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA =
    fpViewer -> GetApplicableVisAttributes (polyhedron.GetVisAttributes ());

  // Get view parameters that the user can force through the vis
  // attributes, thereby over-riding the current view parameter.
  G4ViewParameters::DrawingStyle drawing_style = GetDrawingStyle (pVA);

  // Note that in stored mode, because this call gets embedded in a display
  //  list, it is the colour _at the time of_ creation of the display list, so
  //  even if the colour is changed, for example, by interaction with a Qt
  //  window, current_colour does not change.
  GLfloat current_colour [4];
  glGetFloatv (GL_CURRENT_COLOR, current_colour);
  
  G4bool isTransparent = false;
  if (current_colour[3] < 1.) {  // This object is transparent
    isTransparent = true;
  }

  // This is the colour used to paint surfaces in hlr mode.
  GLfloat clear_colour[4];
  glGetFloatv (GL_COLOR_CLEAR_VALUE, clear_colour);
  
  G4double lineWidth = GetLineWidth(pVA);
  pGLViewer->ChangeLineWidth(lineWidth);

  G4bool isAuxEdgeVisible = GetAuxEdgeVisible (pVA);

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
      glEnable(GL_COLOR_MATERIAL);
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
        glEnable(GL_COLOR_MATERIAL);
	glEnable (GL_CULL_FACE);
	glCullFace (GL_BACK);
	glPolygonMode (GL_FRONT, GL_FILL);
      }
    }
    if (!fProcessing2D) glEnable (GL_LIGHTING);
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
  glBegin (GL_QUADS);
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
	glEdgeFlag (GL_TRUE);
      } else {
	glEdgeFlag (GL_FALSE);
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
      glEdgeFlag (GL_FALSE);
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

    glDisable(GL_COLOR_MATERIAL); // Revert to glMaterial for hlr/sr.

    // Do it all over again (twice) for hlr...
    if  (drawing_style == G4ViewParameters::hlr ||
	 drawing_style == G4ViewParameters::hlhsr) {

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
      GLfloat* painting_colour;
      if  (drawing_style == G4ViewParameters::hlr) {
	if (isTransparent) {
	  // Transparent - don't paint...
	  goto end_of_drawing_through_stencil;
	}
	painting_colour = clear_colour;
      } else {  // drawing_style == G4ViewParameters::hlhsr
	painting_colour = current_colour;
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
      for (int edgeCount = 0; edgeCount < 4; ++edgeCount) {
	if (edgeFlag[edgeCount] > 0) {
	  glEdgeFlag (GL_TRUE);
	} else {
	  glEdgeFlag (GL_FALSE);
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
      glBegin (GL_QUADS);
      for (int edgeCount = 0; edgeCount < 4; ++edgeCount) {
	if (edgeFlag[edgeCount] > 0) {
	  glEdgeFlag (GL_TRUE);
	} else {
	  glEdgeFlag (GL_FALSE);
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
      glBegin (GL_QUADS);      // Ready for next facet.  GL
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

//Method for handling G4NURBS objects for drawing solids.
//Knots and Ctrl Pnts MUST be arrays of GLfloats.
void G4OpenGLSceneHandler::AddPrimitive (const G4NURBS& nurb) {

  GLUnurbsObj *gl_nurb;
  gl_nurb = gluNewNurbsRenderer ();

  GLfloat *u_knot_array, *u_knot_array_ptr;
  u_knot_array = u_knot_array_ptr = new GLfloat [nurb.GetnbrKnots(G4NURBS::U)];
  G4NURBS::KnotsIterator u_iterator (nurb, G4NURBS::U);
  while (u_iterator.pick (u_knot_array_ptr++)){}

  GLfloat *v_knot_array, *v_knot_array_ptr;
  v_knot_array = v_knot_array_ptr = new GLfloat [nurb.GetnbrKnots(G4NURBS::V)];
  G4NURBS::KnotsIterator v_iterator (nurb, G4NURBS::V);
  while (v_iterator.pick (v_knot_array_ptr++)){}

  GLfloat *ctrl_pnt_array, *ctrl_pnt_array_ptr;
  ctrl_pnt_array = ctrl_pnt_array_ptr =
    new GLfloat [nurb.GettotalnbrCtrlPts () * G4NURBS::NofC];
  G4NURBS::CtrlPtsCoordsIterator c_p_iterator (nurb);
  while (c_p_iterator.pick (ctrl_pnt_array_ptr++)){}

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA =
    fpViewer -> GetApplicableVisAttributes (nurb.GetVisAttributes ());

  // Get view parameters that the user can force through the vis
  // attributes, thereby over-riding the current view parameter.
  G4ViewParameters::DrawingStyle drawing_style = GetDrawingStyle (pVA);
  //G4bool isAuxEdgeVisible = GetAuxEdgeVisible (pVA);
  
  //Get colour, etc..
  const G4Colour& c = pVA -> GetColour ();

  switch (drawing_style) {

  case (G4ViewParameters::hlhsr):
    //    G4cout << "Hidden line removal not implememented in G4OpenGL.\n"
    // << "Using hidden surface removal." << G4endl;
  case (G4ViewParameters::hsr):
    {
      if (!fProcessing2D) glEnable (GL_LIGHTING);
      glEnable (GL_DEPTH_TEST);
      glEnable (GL_AUTO_NORMAL);
      glEnable (GL_NORMALIZE);
      gluNurbsProperty (gl_nurb, GLU_DISPLAY_MODE, GLU_FILL);
      gluNurbsProperty (gl_nurb, GLU_SAMPLING_TOLERANCE, 50.0);
      GLfloat materialColour [4];
      materialColour [0] = c.GetRed ();
      materialColour [1] = c.GetGreen ();
      materialColour [2] = c.GetBlue ();
      materialColour [3] = 1.0;  // = c.GetAlpha () for transparency -
				 // but see complication in
				 // AddPrimitive(const G4Polyhedron&).
      glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE, materialColour);
      break;
    }
  case (G4ViewParameters::hlr):
    //    G4cout << "Hidden line removal not implememented in G4OpenGL.\n"
    // << "Using wireframe." << G4endl;
  case (G4ViewParameters::wireframe):
  default:
    glDisable (GL_LIGHTING);
//    glDisable (GL_DEPTH_TEST);
    glEnable (GL_DEPTH_TEST);
    glDisable (GL_AUTO_NORMAL);
    glDisable (GL_NORMALIZE);
    gluNurbsProperty (gl_nurb, GLU_DISPLAY_MODE, GLU_OUTLINE_POLYGON);
    gluNurbsProperty (gl_nurb, GLU_SAMPLING_TOLERANCE, 50.0);
    glColor4d(c.GetRed(), c.GetGreen(), c.GetBlue(),c.GetAlpha());
    break;
  }	

  gluBeginSurface (gl_nurb);
  G4int u_stride = 4;
  G4int v_stride = nurb.GetnbrCtrlPts(G4NURBS::U) * 4;

  gluNurbsSurface (gl_nurb, 
		   nurb.GetnbrKnots (G4NURBS::U), (GLfloat*)u_knot_array, 
		   nurb.GetnbrKnots (G4NURBS::V), (GLfloat*)v_knot_array,
		   u_stride,
		   v_stride,  
		   ctrl_pnt_array,
		   nurb.GetUorder (),
		   nurb.GetVorder (), 
		   GL_MAP2_VERTEX_4);
  
  gluEndSurface (gl_nurb);

  delete [] u_knot_array;  // These should be allocated with smart allocators
  delete [] v_knot_array;  // to avoid memory explosion.
  delete [] ctrl_pnt_array;

  gluDeleteNurbsRenderer (gl_nurb);
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

#endif
