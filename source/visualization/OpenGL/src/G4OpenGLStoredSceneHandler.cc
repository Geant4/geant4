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
// $Id: G4OpenGLStoredSceneHandler.cc,v 1.46 2010-11-10 17:11:20 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// OpenGL stored scene - creates OpenGL display lists.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

// Included here - problems with HP compiler if not before other includes?
#include "G4NURBS.hh"

// Here follows a special for Mesa, the OpenGL emulator.  Does not affect
// other OpenGL's, as far as I'm aware.   John Allison 18/9/96.
#define CENTERLINE_CLPP  /* CenterLine C++ workaround: */
// Also seems to be required for HP's CC and AIX xlC, at least.

#include "G4OpenGLStoredSceneHandler.hh"

#include "G4PhysicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Polyhedron.hh"
#include "G4AttHolder.hh"
#include "G4OpenGLTransform3D.hh"
#include "G4OpenGLViewer.hh"

G4OpenGLStoredSceneHandler::PO::PO
(G4int id,
 const G4Transform3D& tr):
  fDisplayListId(id),
  fTransform(tr),
  fPickName(0)
{}

G4OpenGLStoredSceneHandler::TO::TO
(G4int id,
 const G4Transform3D& tr):
  fDisplayListId(id),
  fTransform(tr),
  fPickName(0),
  fStartTime(-DBL_MAX),
  fEndTime(DBL_MAX)
{}

G4OpenGLStoredSceneHandler::G4OpenGLStoredSceneHandler
(G4VGraphicsSystem& system,
 const G4String& name):
G4OpenGLSceneHandler (system, fSceneIdCount++, name),
fAddPrimitivePreambleNestingDepth (0),
fTopPODL (0)
{}

G4OpenGLStoredSceneHandler::~G4OpenGLStoredSceneHandler ()
{}

void G4OpenGLStoredSceneHandler::AddPrimitivePreamble(const G4Visible& visible)
{
  // Track nesting depth to avoid recursive calls, for example, from a
  // G4Polymarker that invokes a G4Circle...
  fAddPrimitivePreambleNestingDepth++;
  if (fAddPrimitivePreambleNestingDepth > 1) return;

  // Because of our need to control colour of transients (display by
  // time fading), display lists may only cover a single primitive.
  // So display list setup is here.

  if (fpViewer->GetViewParameters().IsPicking()) {
    fPickMap[++fPickName] = 0;
  }

  const G4Colour& c = GetColour (visible);

  if (fMemoryForDisplayLists) {
    fDisplayListId = glGenLists (1);
    if (glGetError() == GL_OUT_OF_MEMORY ||
	fDisplayListId > fDisplayListLimit) {
      G4cout <<
  "********************* WARNING! ********************"
  "\n*  Display list limit reached in OpenGL."
  "\n*  Continuing drawing WITHOUT STORING. Scene only partially refreshable."
  "\n*  Current limit: " << fDisplayListLimit <<
  ".  Change with \"/vis/ogl/set/displayListLimit\"."
  "\n***************************************************"
	     << G4endl;
      fMemoryForDisplayLists = false;
    }
  }
  if (fMemoryForDisplayLists) {
    if (fReadyForTransients) {
      TO to(fDisplayListId, *fpObjectTransformation);
      to.fPickName = fPickName;
      to.fColour = c;
      const G4VisAttributes* pVA =
	fpViewer->GetApplicableVisAttributes(visible.GetVisAttributes());
      to.fStartTime = pVA->GetStartTime();
      to.fEndTime = pVA->GetEndTime();
      fTOList.push_back(to);
      glDrawBuffer (GL_FRONT);
      glPushMatrix();
      G4OpenGLTransform3D oglt (*fpObjectTransformation);
      glMultMatrixd (oglt.GetGLMatrix ());
      glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
      glNewList (fDisplayListId, GL_COMPILE_AND_EXECUTE);
    }
    else {
      PO po(fDisplayListId, *fpObjectTransformation);
      po.fPickName = fPickName;
      fPOList.push_back(po);
      glNewList (fDisplayListId, GL_COMPILE);
      glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
    }
  } else {
    glDrawBuffer (GL_FRONT);
    glPushMatrix();
    G4OpenGLTransform3D oglt (*fpObjectTransformation);
    glMultMatrixd (oglt.GetGLMatrix ());
    glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
  }

  if (fProcessing2D) {
    // Push current 3D world matrices and load identity to define screen
    // coordinates...
    glMatrixMode (GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho (-1., 1., -1., 1., -G4OPENGL_FLT_BIG, G4OPENGL_FLT_BIG);
    glMatrixMode (GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    G4OpenGLTransform3D oglt (*fpObjectTransformation);
    glMultMatrixd (oglt.GetGLMatrix ());
    glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
  }
}

void G4OpenGLStoredSceneHandler::AddPrimitivePostamble()
{
  if (fProcessing2D) {
    // Pop current 3D world matrices back again...
    glMatrixMode (GL_PROJECTION);
    glPopMatrix();
    glMatrixMode (GL_MODELVIEW);
    glPopMatrix();
  }

  //  if ((glGetError() == GL_TABLE_TOO_LARGE) || (glGetError() == GL_OUT_OF_MEMORY)) {  // Could close?
  if (glGetError() == GL_OUT_OF_MEMORY) {  // Could close?
    G4cout <<
      "ERROR: G4OpenGLStoredSceneHandler::EndModeling: Failure to allocate"
      "  display List for fTopPODL - try OpenGL Immediated mode."
           << G4endl;
  }
  if (fMemoryForDisplayLists) {
    glEndList();
    if (glGetError() == GL_OUT_OF_MEMORY) {  // Could close?
      G4cout <<
        "ERROR: G4OpenGLStoredSceneHandler::EndModeling: Failure to allocate"
        "  display List for fTopPODL - try OpenGL Immediated mode."
             << G4endl;
    }
  }
  if (fReadyForTransients || !fMemoryForDisplayLists) {
    glPopMatrix();
    glFlush ();
    glDrawBuffer (GL_BACK);
  }
  fAddPrimitivePreambleNestingDepth--;
}

void G4OpenGLStoredSceneHandler::AddPrimitive (const G4Polyline& polyline)
{
  AddPrimitivePreamble(polyline);
  G4OpenGLSceneHandler::AddPrimitive(polyline);
  AddPrimitivePostamble();
}

void G4OpenGLStoredSceneHandler::AddPrimitive (const G4Polymarker& polymarker)
{
  AddPrimitivePreamble(polymarker);
  G4OpenGLSceneHandler::AddPrimitive(polymarker);
  AddPrimitivePostamble();
}

void G4OpenGLStoredSceneHandler::AddPrimitive (const G4Text& text)
{
  // Note: colour is still handled in
  // G4OpenGLSceneHandler::AddPrimitive(const G4Text&), so it still
  // gets into the display list
  AddPrimitivePreamble(text);
  G4OpenGLSceneHandler::AddPrimitive(text);
  AddPrimitivePostamble();
}

void G4OpenGLStoredSceneHandler::AddPrimitive (const G4Circle& circle)
{
  AddPrimitivePreamble(circle);
  G4OpenGLSceneHandler::AddPrimitive(circle);
  AddPrimitivePostamble();
}

void G4OpenGLStoredSceneHandler::AddPrimitive (const G4Square& square)
{
  AddPrimitivePreamble(square);
  G4OpenGLSceneHandler::AddPrimitive(square);
  AddPrimitivePostamble();
}

void G4OpenGLStoredSceneHandler::AddPrimitive (const G4Scale& scale)
{
  // Let base class split into primitives.
  G4OpenGLSceneHandler::AddPrimitive(scale);
}

void G4OpenGLStoredSceneHandler::AddPrimitive (const G4Polyhedron& polyhedron)
{
  // Note: colour is still handled in
  // G4OpenGLSceneHandler::AddPrimitive(const G4Polyhedron&), so it still
  // gets into the display list
  AddPrimitivePreamble(polyhedron);
  G4OpenGLSceneHandler::AddPrimitive(polyhedron);
  AddPrimitivePostamble();
}

void G4OpenGLStoredSceneHandler::AddPrimitive (const G4NURBS& nurbs)
{
  // Note: colour is still handled in
  // G4OpenGLSceneHandler::AddPrimitive(const G4NURBS&), so it still
  // gets into the display list
  AddPrimitivePreamble(nurbs);
  G4OpenGLSceneHandler::AddPrimitive(nurbs);
  AddPrimitivePostamble();
}

void G4OpenGLStoredSceneHandler::BeginPrimitives
(const G4Transform3D& objectTransformation)
{  
  G4OpenGLSceneHandler::BeginPrimitives (objectTransformation);

  // Display list setup moved to AddPrimitivePreamble.  See notes there.
}

void G4OpenGLStoredSceneHandler::EndPrimitives ()
{
  G4OpenGLSceneHandler::EndPrimitives ();
}

void G4OpenGLStoredSceneHandler::BeginPrimitives2D
(const G4Transform3D& objectTransformation)
{
  G4OpenGLSceneHandler::BeginPrimitives2D(objectTransformation);
}

void G4OpenGLStoredSceneHandler::EndPrimitives2D ()
{
  G4OpenGLSceneHandler::EndPrimitives2D ();
}

void G4OpenGLStoredSceneHandler::BeginModeling () {
  G4VSceneHandler::BeginModeling();
  ClearStore();  // ...and all that goes with it.
  /* Debug...
  fDisplayListId = glGenLists (1);
  G4cout << "OGL::fDisplayListId (start): " << fDisplayListId << G4endl;
  */
}

void G4OpenGLStoredSceneHandler::EndModeling () {
  // Make a List which calls the other lists.
  fTopPODL = glGenLists (1);
  if (glGetError() == GL_OUT_OF_MEMORY) {  // Could pre-allocate?
    G4cout <<
      "ERROR: G4OpenGLStoredSceneHandler::EndModeling: Failure to allocate"
      "  display List for fTopPODL - try OpenGL Immediated mode."
	   << G4endl;
  } else {
    glNewList (fTopPODL, GL_COMPILE_AND_EXECUTE); {
      for (size_t i = 0; i < fPOList.size (); i++) {
	glPushMatrix();
	G4OpenGLTransform3D oglt (fPOList[i].fTransform);
	glMultMatrixd (oglt.GetGLMatrix ());
	if (fpViewer->GetViewParameters().IsPicking())
	  glLoadName(fPOList[i].fPickName);
	glCallList (fPOList[i].fDisplayListId);
	glPopMatrix();
      }
    }
    glEndList ();
    if (glGetError() == GL_OUT_OF_MEMORY) {  // Could close?
      G4cout <<
        "ERROR: G4OpenGLStoredSceneHandler::EndModeling: Failure to allocate"
        "  display List for fTopPODL - try OpenGL Immediated mode."
             << G4endl;
    }
  }

  G4VSceneHandler::EndModeling ();

  /* Debug...
  fDisplayListId = glGenLists (1);
  G4cout << "OGL::fDisplayListId (end): " << fDisplayListId << G4endl;
  */
}

void G4OpenGLStoredSceneHandler::ClearStore () {

  //G4cout << "G4OpenGLStoredSceneHandler::ClearStore" << G4endl;

  G4VSceneHandler::ClearStore ();  // Sets need kernel visit, etc.

  // Delete OpenGL permanent display lists.
  for (size_t i = 0; i < fPOList.size (); i++)
    glDeleteLists (fPOList[i].fDisplayListId, 1);
  if (fTopPODL) glDeleteLists (fTopPODL, 1);
  fTopPODL = 0;

  // Clear other lists, dictionary, etc.
  fPOList.clear ();
  fSolidMap.clear ();
  ClearAndDestroyAtts();

  // ...and clear transient store...
  for (size_t i = 0; i < fTOList.size (); i++)
    glDeleteLists(fTOList[i].fDisplayListId, 1);
  fTOList.clear ();

  fMemoryForDisplayLists = true;
}

void G4OpenGLStoredSceneHandler::ClearTransientStore () {

  //G4cout << "G4OpenGLStoredSceneHandler::ClearTransientStore" << G4endl;

  G4VSceneHandler::ClearTransientStore ();

  // Delete OpenGL transient display lists and Transient Objects themselves.
  for (size_t i = 0; i < fTOList.size (); i++)
    glDeleteLists(fTOList[i].fDisplayListId, 1);
  fTOList.clear ();

  fMemoryForDisplayLists = true;

  // Make sure screen corresponds to graphical database...
  if (fpViewer) {
    fpViewer -> SetView ();
    fpViewer -> ClearView ();
    fpViewer -> DrawView ();
  }
}

void G4OpenGLStoredSceneHandler::RequestPrimitives (const G4VSolid& solid)
{
  if (fReadyForTransients) {
    // Always draw transient solids, e.g., hits represented as solids.
    // (As we have no control over the order of drawing of transient
    // objects, we cannot do anything about transparent ones, as
    // below, so always draw them.)
    G4VSceneHandler::RequestPrimitives (solid);
    return;
  }

  // For non-transient (run-duration) objects, ensure transparent
  // objects are drawn last.  The problem of
  // blending/transparency/alpha is quite a tricky one - see History
  // of opengl-V07-01-01/2/3.
  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA =
    fpViewer -> GetApplicableVisAttributes(fpVisAttribs);
  const G4Colour& c = pVA -> GetColour ();
  G4double opacity = c.GetAlpha ();

  if (!fSecondPass) {
    G4bool transparency_enabled = true;
    G4OpenGLViewer* pViewer = dynamic_cast<G4OpenGLViewer*>(fpViewer);
    if (pViewer) transparency_enabled = pViewer->transparency_enabled;
    if (transparency_enabled && opacity < 1.) {
      // On first pass, transparent objects are not drawn, but flag is set...
      fSecondPassRequested = true;
      return;
    }
  }

  // On second pass, opaque objects are not drwan...
  if (fSecondPass && opacity >= 1.) return;

  G4PhysicalVolumeModel* pPVModel =
    dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
 
  if (pPVModel) {
    // If part of the geometry hierarchy, i.e., from a
    // G4PhysicalVolumeModel, check if a display list already exists for
    // this solid, re-use it if possible.  We could be smarter, and
    // recognise repeated branches of the geometry hierarchy, for
    // example.  But this algorithm should be secure, I think...
    const G4VSolid* pSolid = &solid;
    EAxis axis = kRho;
    G4VPhysicalVolume* pCurrentPV = pPVModel->GetCurrentPV();
    if (pCurrentPV -> IsReplicated ()) {
      G4int nReplicas;
      G4double width;
      G4double offset;
      G4bool consuming;
      pCurrentPV->GetReplicationData(axis,nReplicas,width,offset,consuming);
    }
    // Provided it is not parametrised (because if so, the
    // solid's parameters might have been changed)...
    if (!(pCurrentPV -> IsParameterised ()) &&
	// Provided it is not replicated radially (because if so, the
	// solid's parameters will have been changed)...
	!(pCurrentPV -> IsReplicated () && axis == kRho) &&
	// ...and if the solid has already been rendered...
	(fSolidMap.find (pSolid) != fSolidMap.end ())) {
      fDisplayListId = fSolidMap [pSolid];
      PO po(fDisplayListId,*fpObjectTransformation);
      if (fpViewer->GetViewParameters().IsPicking()) {
	G4AttHolder* holder = new G4AttHolder;
	// Load G4Atts from G4VisAttributes, if any...
	const G4VisAttributes* va = pPVModel->GetCurrentLV()->GetVisAttributes();
	if (va) {
	  const std::map<G4String,G4AttDef>* vaDefs = va->GetAttDefs();
	  if (vaDefs) holder->AddAtts(va->CreateAttValues(), vaDefs);
	}
	// Load G4Atts from G4PhysicalVolumeModel...
	const std::map<G4String,G4AttDef>* defs = pPVModel->GetAttDefs();
	if (defs) holder->AddAtts(pPVModel->CreateCurrentAttValues(), defs);
	fPickMap[++fPickName] = holder;
	po.fPickName = fPickName;
      }
      fPOList.push_back(po);
    }
    else {
      G4VSceneHandler::RequestPrimitives (solid);
      fSolidMap [pSolid] = fDisplayListId;
    }
    return;
  }

  // Otherwise invoke base class method...
  G4VSceneHandler::RequestPrimitives (solid);
}

G4int G4OpenGLStoredSceneHandler::fSceneIdCount = 0;

G4int  G4OpenGLStoredSceneHandler::fDisplayListId = 0;
G4bool G4OpenGLStoredSceneHandler::fMemoryForDisplayLists = true;
G4int  G4OpenGLStoredSceneHandler::fDisplayListLimit = 50000;

#endif
