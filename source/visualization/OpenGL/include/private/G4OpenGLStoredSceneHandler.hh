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
// Andrew Walkden  10th February 1997
// G4OpenGLStoredSceneHandler - creates OpenGL Display lists.

#ifndef G4OPENGLSTOREDSCENEHANDLER_HH
#define G4OPENGLSTOREDSCENEHANDLER_HH

#include "globals.hh"
#include "G4OpenGLSceneHandler.hh"
#include "G4Text.hh"
#include <map>
#include <vector>

class G4OpenGLStored;

class G4OpenGLStoredSceneHandler: public G4OpenGLSceneHandler {

  friend class G4OpenGLStoredViewer;  // ..allows access to P/TODLs.

public:

  G4OpenGLStoredSceneHandler (G4VGraphicsSystem& system, const G4String& name = "");
  virtual ~G4OpenGLStoredSceneHandler ();
  void BeginPrimitives (const G4Transform3D& objectTransformation);
  void EndPrimitives ();
  void BeginPrimitives2D (const G4Transform3D& objectTransformation);
  void EndPrimitives2D ();
  void BeginModeling ();
  void EndModeling ();
  using G4VSceneHandler::AddPrimitive;
  void AddPrimitive (const G4Polyline&);
  void AddPrimitive (const G4Polymarker&);
  void AddPrimitive (const G4Circle&);
  void AddPrimitive (const G4Square&);
  void AddPrimitive (const G4Text&);
  void AddPrimitive (const G4Polyhedron&);
  void ClearStore ();
  void ClearTransientStore ();

protected:

  G4bool AddPrimitivePreamble(const G4VMarker& visible);
  G4bool AddPrimitivePreamble(const G4Polyline& visible);
  G4bool AddPrimitivePreamble(const G4Polyhedron& visible);
  // Return false if no further processing required.

  void AddPrimitivePostamble();

  // Two virtual functions for extra processing in a sub-class, for
  // example, to make a display tree.  They are to return true if the
  // visible object uses gl commands for drawing.  This is
  // predominantly true; a notable exception is Qt text.  In that
  // case, a display list does not need to be created; all relevant
  // information is assumed to be stored in the PO/TOList.
  virtual G4bool ExtraPOProcessing
  (const G4Visible&, size_t /*currentPOListIndex*/) {return true;}
  virtual G4bool ExtraTOProcessing
  (const G4Visible&, size_t /*currentTOListIndex*/) {return true;}

  static G4int fSceneIdCount;   // static counter for OpenGLStored scenes.

  // Display list management.  Static since there's only one OGL store.
  // Used to link a TODL or PODL to a display list.
  static G4int fDisplayListId;

  // Under some circumstances we need to prevent use of a display list.
  // For example, a transient display of the current time window during
  // a sequence of evolving time windows. This is set in
  // G4OpenGLStoredViewer::AddPrimitiveForASingleFrame and acted upon in
  // AddPrimitivePreambleInternal.
  G4bool fDoNotUseDisplayList;  // Avoid display list use if true

  // PODL = Persistent Object Display List.
  // This "top PODL" was made redundant when the PO list was
  // "unwrapped" 27th October 2011, but keep it for now in case we
  // need to wrap it again.
  GLint  fTopPODL;              // List which calls the other PODLs.

  // G4Text plus transform and 2/3D.
  struct G4TextPlus {
    G4TextPlus(const G4Text& text): fG4Text(text), fProcessing2D(false) {}
    G4Text fG4Text;
    G4bool fProcessing2D;
  };

  // PO = Persistent Object, i.e., run-durantion object, e.g., geometry.
  struct PO {
    PO();
    PO(const PO&);
    PO(G4int id, const G4Transform3D& tr = G4Transform3D());
    ~PO();
    PO& operator= (const PO&);
    G4int fDisplayListId;
    G4Transform3D fTransform;
    GLuint fPickName;
    G4Colour fColour;
    G4TextPlus* fpG4TextPlus;
    G4bool fMarkerOrPolyline;
  };
  std::vector<PO> fPOList; 
  
  // TO = Transient Object, e.g., trajectories.
  struct TO {
    TO();
    TO(const TO&);
    TO(G4int id, const G4Transform3D& tr = G4Transform3D());
    ~TO();
    TO& operator= (const TO&);
    G4int fDisplayListId;
    G4Transform3D fTransform;
    GLuint fPickName;
    G4double fStartTime, fEndTime;  // Time range (e.g., for trajectory steps).
    G4Colour fColour;
    G4TextPlus* fpG4TextPlus;
    G4bool fMarkerOrPolyline;
  };
  std::vector<TO> fTOList; 
  
  // Stop-gap solution of structure re-use.
  // A proper implementation would use geometry hierarchy.
  std::map <const G4VSolid*, G4int, std::less <const G4VSolid*> > fSolidMap;

private:
  bool AddPrimitivePreambleInternal(const G4Visible& visible, bool isMarker, bool isPolyline);

};

#endif
