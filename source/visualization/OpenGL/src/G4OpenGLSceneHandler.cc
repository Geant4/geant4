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
// $Id: G4OpenGLSceneHandler.cc,v 1.59 2010-05-30 09:53:05 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4OpenGLFontBaseStore.hh"
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

G4OpenGLSceneHandler::G4OpenGLSceneHandler (G4VGraphicsSystem& system,
			      G4int id,
			      const G4String& name):
  G4VSceneHandler (system, id, name),
  fPickName(0),
  fProcessing2D (false),
  fProcessingPolymarker(false)
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
  fProcessing2D = true;
}

void G4OpenGLSceneHandler::EndPrimitives2D ()
{
  fProcessing2D = false;
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

  // Loads G4Atts for picking...
  if (fpViewer->GetViewParameters().IsPicking()) {
    G4AttHolder* holder = new G4AttHolder;
    LoadAtts(line, holder);
    fPickMap[fPickName] = holder;
  }

  // Note: colour treated in sub-class.

  if (fpViewer -> GetViewParameters ().IsMarkerNotHidden ())
    glDisable (GL_DEPTH_TEST);
  else {glEnable (GL_DEPTH_TEST); glDepthFunc (GL_LESS);}

  glDisable (GL_LIGHTING);

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA =
    fpViewer -> GetApplicableVisAttributes (line.GetVisAttributes ());

  G4double lineWidth = GetLineWidth(pVA);
  glLineWidth(lineWidth);

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
  G4int nPoints = polymarker.size ();
  if (nPoints <= 0) return;

  fProcessingPolymarker = true;

  // Loads G4Atts for picking...
  if (fpViewer->GetViewParameters().IsPicking()) {
    G4AttHolder* holder = new G4AttHolder;
    LoadAtts(polymarker, holder);
    fPickMap[fPickName] = holder;
  }

  switch (polymarker.GetMarkerType()) {
  default:
  case G4Polymarker::dots:
    {
      for (size_t iPoint = 0; iPoint < polymarker.size (); iPoint++) {
        G4Circle dot (polymarker);
        dot.SetPosition (polymarker[iPoint]);
        dot.SetWorldSize  (0.);
        dot.SetScreenSize (0.1);  // Very small circle.
        G4OpenGLSceneHandler::AddPrimitive (dot);
      }
    }
    break;
  case G4Polymarker::circles:
    {
      std::vector <G4VMarker> circleV; 
      for (size_t iPoint = 0; iPoint < polymarker.size (); iPoint++) {
        G4Circle circle (polymarker);
        // If not already drawn
        circle.SetPosition (polymarker[iPoint]);
        circleV.push_back(circle);
        //      G4OpenGLSceneHandler::AddPrimitive (circle);
      }
      G4OpenGLSceneHandler::AddPrimitives (circleV);
    }
    break;
  case G4Polymarker::squares:
    {
      std::vector <G4VMarker> squareV; 
      for (size_t iPoint = 0; iPoint < polymarker.size (); iPoint++) {
        G4Square square (polymarker);
        square.SetPosition (polymarker[iPoint]);
        squareV.push_back(square);
        //      G4OpenGLSceneHandler::AddPrimitive (square);
      }
      G4OpenGLSceneHandler::AddPrimitives (squareV);
    }
    break;
  }

  fProcessingPolymarker = false;
}

void G4OpenGLSceneHandler::AddPrimitive (const G4Text& text) {

  // Loads G4Atts for picking...
  if (fpViewer->GetViewParameters().IsPicking()) {
    G4AttHolder* holder = new G4AttHolder;
    LoadAtts(text, holder);
    fPickMap[fPickName] = holder;
  }

  const G4Colour& c = GetTextColour (text);  // Picks up default if none.
  MarkerSizeType sizeType;
  G4double size = GetMarkerSize (text, sizeType);
  G4ThreeVector position (text.GetPosition ());
  G4String textString = text.GetText();

  G4int font_base = G4OpenGLFontBaseStore::GetFontBase(fpViewer,size);
  if (font_base < 0) {
    static G4int callCount = 0;
    ++callCount;
    if (callCount <= 10 || callCount%100 == 0) {
      G4cout <<
	"G4OpenGLSceneHandler::AddPrimitive (const G4Text&) call count "
	     << callCount <<
	"\n  No fonts available."
	"\n  Called with text \""
	     << text.GetText ()
	     << "\"\n  at " << position
	     << ", size " << size
	     << ", offsets " << text.GetXOffset () << ", " << text.GetYOffset ()
	     << ", type " << G4int(sizeType)
	     << ", colour " << c
	     << G4endl;
    }
    return;
  }
  const char* textCString = textString.c_str();
  glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
  glDisable (GL_DEPTH_TEST);
  glDisable (GL_LIGHTING);
  
  glRasterPos3d(position.x(),position.y(),position.z());
  // No action on offset or layout at present.
   glPushAttrib(GL_LIST_BIT);
   glListBase(font_base);
   glCallLists(strlen(textCString), GL_UNSIGNED_BYTE, (GLubyte *)textCString);
   glPopAttrib();
}

void G4OpenGLSceneHandler::AddPrimitive (const G4Circle& circle) {
  glEnable (GL_POINT_SMOOTH);
  AddCircleSquare (circle, G4OpenGLBitMapStore::circle);
}

void G4OpenGLSceneHandler::AddPrimitive (const G4Square& square) {
  glDisable (GL_POINT_SMOOTH);
  AddCircleSquare (square, G4OpenGLBitMapStore::square);
}

void G4OpenGLSceneHandler::AddPrimitives (std::vector <G4VMarker> square) {
  glDisable (GL_POINT_SMOOTH);
  AddCircleSquareVector (square, G4OpenGLBitMapStore::square);
}

void G4OpenGLSceneHandler::AddCircleSquare
(const G4VMarker& marker,
 G4OpenGLBitMapStore::Shape shape) {

  std::vector <G4VMarker> circleVector;
  circleVector.push_back(marker);
  AddCircleSquareVector(circleVector,shape);
}
 
void G4OpenGLSceneHandler::AddCircleSquareVector
(std::vector <G4VMarker> marker,
 G4OpenGLBitMapStore::Shape shape) {

  if (marker.size() == 0) {
    return;
  }

  if (!fProcessingPolymarker) {  // Polymarker has already loaded atts.
    // Loads G4Atts for picking...
    if (fpViewer->GetViewParameters().IsPicking()) {
      G4AttHolder* holder = new G4AttHolder;
      LoadAtts(marker[0], holder);
      fPickMap[fPickName] = holder;
    }
  }

  // Note: colour treated in sub-class.

  if (fpViewer -> GetViewParameters ().IsMarkerNotHidden ()) {
    glDisable (GL_DEPTH_TEST);
  } else {
    glEnable (GL_DEPTH_TEST); glDepthFunc (GL_LESS);
  }
  
  glDisable (GL_LIGHTING);
  
  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA =
    fpViewer -> GetApplicableVisAttributes (marker[0].GetVisAttributes ());

  G4double lineWidth = GetLineWidth(pVA);
  glLineWidth(lineWidth);

  G4VMarker::FillStyle style = marker[0].GetFillStyle();

  G4bool filled = false;
  static G4bool hashedWarned = false;
  
  switch (style) {
  case G4VMarker::noFill: 
    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
    filled = false;
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
    filled = true;
    break;
    
  }



  MarkerSizeType sizeType;
  G4double size = GetMarkerSize(marker[0], sizeType);

  // Draw...
   if (sizeType == world) {  // Size specified in world coordinates.

     for (unsigned int a=0;a<marker.size();a++) {
       G4Point3D centre = marker[a].GetPosition();
       // A few useful quantities...
       DrawXYPolygon (shape, size, centre, pVA);
     }
   } else { // Size specified in screen (window) coordinates.
     // A few useful quantities...
     glPointSize (size);
     glBegin (GL_POINTS);
     for (unsigned int a=0;a<marker.size();a++) {
       G4Point3D centre = marker[a].GetPosition();
       glVertex3f(centre.x(),centre.y(),centre.z());
     }
     glEnd();
     //Antialiasing
     glEnable (GL_POINT_SMOOTH);
     //Transparency
     glEnable(GL_BLEND);
     glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

     // L. GARNIER 1 March 2009
     // Old method, we draw a bitmap instead of a GL_POINT. 
     // I remove it because it cost in term of computing performances
     // and gl2ps can't draw bitmaps

     //      glRasterPos3d(centre.x(),centre.y(),centre.z());
     //      const GLubyte* marker =
     //        G4OpenGLBitMapStore::GetBitMap(shape, size, filled);
     //      glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
     //      glBitmap(GLsizei(size), GLsizei(size), size/2., size/2., 0., 0., marker);
   }
}

void G4OpenGLSceneHandler::DrawXYPolygon
(G4OpenGLBitMapStore::Shape shape,
 G4double size,
 const G4Point3D& centre,
 const G4VisAttributes* pApplicableVisAtts)
{
  G4int nSides;
  G4double startPhi;
  if (shape == G4OpenGLBitMapStore::circle) {
    nSides = GetNoOfSides(pApplicableVisAtts);
    startPhi = 0.;
  } else {
    nSides = 4;
    startPhi = -pi / 4.;
  }

  const G4Vector3D& viewpointDirection =
    fpViewer -> GetViewParameters().GetViewpointDirection();
  const G4Vector3D& up = fpViewer->GetViewParameters().GetUpVector();
  const G4double dPhi = twopi / nSides;
  const G4double radius = size / 2.;
  G4Vector3D start = radius * (up.cross(viewpointDirection)).unit();
  G4double phi;
  G4int i;

  glBegin (GL_POLYGON);
  for (i = 0, phi = startPhi; i < nSides; i++, phi += dPhi) {
    G4Vector3D r = start; r.rotate(phi, viewpointDirection);
    G4Vector3D p = centre + r;
    glVertex3d (p.x(), p.y(), p.z());
  }
  glEnd ();
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

  // Loads G4Atts for picking...
  if (fpViewer->GetViewParameters().IsPicking()) {
    G4AttHolder* holder = new G4AttHolder;
    LoadAtts(polyhedron, holder);
    fPickMap[fPickName] = holder;
  }

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA =
    fpViewer -> GetApplicableVisAttributes (polyhedron.GetVisAttributes ());

  // Get view parameters that the user can force through the vis
  // attributes, thereby over-riding the current view parameter.
  G4ViewParameters::DrawingStyle drawing_style = GetDrawingStyle (pVA);

  //Get colour, etc...
  G4bool transparency_enabled = true;
  G4OpenGLViewer* pViewer = dynamic_cast<G4OpenGLViewer*>(fpViewer);
  if (pViewer) transparency_enabled = pViewer->transparency_enabled;
  const G4Colour& c = pVA->GetColour();
  GLfloat materialColour [4];
  materialColour [0] = c.GetRed ();
  materialColour [1] = c.GetGreen ();
  materialColour [2] = c.GetBlue ();
  if (transparency_enabled) {
    materialColour [3] = c.GetAlpha ();
  } else {
    materialColour [3] = 1.;
  }

  G4double lineWidth = GetLineWidth(pVA);
  glLineWidth(lineWidth);

  GLfloat clear_colour[4];
  glGetFloatv (GL_COLOR_CLEAR_VALUE, clear_colour);

  G4bool isAuxEdgeVisible = GetAuxEdgeVisible (pVA);

  G4bool clipping = pViewer->fVP.IsSection() || pViewer->fVP.IsCutaway();

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
    if (materialColour[3] < 1.) {
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
    glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
    break;
  case (G4ViewParameters::hsr):
    glEnable (GL_DEPTH_TEST);
    glDepthFunc (GL_LEQUAL);    
    if (materialColour[3] < 1.) {
      // Transparent...
      glDepthMask (0);  // Make depth buffer read-only.
      glDisable (GL_CULL_FACE);
      glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, materialColour);
    } else {
      // Opaque...
      glDepthMask (1);  // Make depth buffer writable (default).
      if (clipping) {
	glDisable (GL_CULL_FACE);
	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
      } else {
	glEnable (GL_CULL_FACE);
	glCullFace (GL_BACK);
	glPolygonMode (GL_FRONT, GL_FILL);
      }
      glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE, materialColour);
    }
    if (!fProcessing2D) glEnable (GL_LIGHTING);
    break;
  case (G4ViewParameters::wireframe):
  default:
    glEnable (GL_DEPTH_TEST);
    glDepthFunc (GL_LEQUAL);    //??? was GL_ALWAYS
    glDisable (GL_CULL_FACE);
    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
    glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
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
    G4int n;
    notLastFace = polyhedron.GetNextFacet(n, vertex, edgeFlag, normals);

    //Loop through the four edges of each G4Facet...
    G4int edgeCount = 0;
    for(edgeCount = 0; edgeCount < n; ++edgeCount) {
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
    if (n == 3) {
      edgeCount = 3;
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
    if (n > 4) {
      G4cerr <<
	"G4OpenGLSceneHandler::AddPrimitive(G4Polyhedron): WARNING";
      G4PhysicalVolumeModel* pPVModel =
	dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
      if (pPVModel) {
	G4VPhysicalVolume* pCurrentPV = pPVModel->GetCurrentPV();
	G4LogicalVolume* pCurrentLV = pPVModel->GetCurrentLV();
	G4cerr <<
	"\n  Volume " << pCurrentPV->GetName() <<
	", Solid " << pCurrentLV->GetSolid()->GetName() <<
	  " (" << pCurrentLV->GetSolid()->GetEntityType();
      }
      G4cerr<<
	"\n   G4Polyhedron facet with " << n << " edges" << G4endl;
    }

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
      if (materialColour[3] < 1.) {
	// Transparent...
	glDepthMask (0);  // Make depth buffer read-only.
	glDisable (GL_CULL_FACE);
	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
      } else {
	// Opaque...
	glDepthMask (1);  // Make depth buffer writable (default).
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
	if (materialColour[3] < 1.) {
	  // Transparent - don't paint...
	  goto end_of_drawing_through_stencil;
	}
	painting_colour = clear_colour;
      } else {  // drawing_style == G4ViewParameters::hlhsr
	painting_colour = materialColour;
      }
      if (materialColour[3] < 1.) {
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
      if (materialColour[3] < 1.) {
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
      glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
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
  glDepthMask (1);              // Revert to default for next primitive.
  glDisable (GL_LIGHTING);      // Revert to default for next primitive.
}

//Method for handling G4NURBS objects for drawing solids.
//Knots and Ctrl Pnts MUST be arrays of GLfloats.
void G4OpenGLSceneHandler::AddPrimitive (const G4NURBS& nurb) {

  // Loads G4Atts for picking...
  if (fpViewer->GetViewParameters().IsPicking()) {
    G4AttHolder* holder = new G4AttHolder;
    LoadAtts(nurb, holder);
    fPickMap[fPickName] = holder;
  }

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
    glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
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
