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
// $Id: G4OpenGLStoredSceneHandler.hh,v 1.32 2010-11-10 17:10:49 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// G4OpenGLStoredSceneHandler - creates OpenGL Display lists.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLSTOREDSCENEHANDLER_HH
#define G4OPENGLSTOREDSCENEHANDLER_HH

#include "globals.hh"
#include "G4OpenGLSceneHandler.hh"
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
  void AddPrimitive (const G4Polyline&);
  void AddPrimitive (const G4Polymarker&);
  void AddPrimitive (const G4Circle&);
  void AddPrimitive (const G4Square&);
  void AddPrimitive (const G4Text&);
  void AddPrimitive (const G4Scale&);
  void AddPrimitive (const G4Polyhedron&);
  void AddPrimitive (const G4NURBS&);
  void ClearStore ();
  void ClearTransientStore ();

  static G4int GetDisplayListLimit() {return fDisplayListLimit;}
  static void SetDisplayListLimit(G4int lim) {fDisplayListLimit = lim;}

protected:

  void RequestPrimitives (const G4VSolid& solid);
  void AddPrimitivePreamble(const G4Visible& visible);
  void AddPrimitivePostamble();

  static G4int  fSceneIdCount;   // static counter for OpenGLStored scenes.
  // Display list management.  All static since there's only one OGL store.
  static G4int  fDisplayListId;  // Workspace.
  static G4bool fMemoryForDisplayLists;  // avoid memory overflow
  static G4int  fDisplayListLimit;       // avoid memory overflow
  G4int fAddPrimitivePreambleNestingDepth;
  
  // PODL = Persistent Object Display List.
  GLint  fTopPODL;                  // List which calls the other PODLs.
  // PO = Persistent Object, i.e., run-durantion object, e.g., geometry.
  struct PO {
    PO(G4int id, const G4Transform3D& tr = G4Transform3D());
    G4int fDisplayListId;
    G4Transform3D fTransform;
    GLuint fPickName;
  };
  std::vector<PO> fPOList; 
  
  // TO = Transient Object, e.g., trajectories.
  struct TO {
    TO(G4int id, const G4Transform3D& tr = G4Transform3D());
    G4int fDisplayListId;
    G4Transform3D fTransform;
    GLuint fPickName;
    G4double fStartTime, fEndTime;  // Time range (e.g., for trajectory steps).
    G4Colour fColour;
  };
  std::vector<TO> fTOList; 
  
  // Stop-gap solution of structure re-use.
  // A proper implementation would use geometry hierarchy.
  std::map <const G4VSolid*, G4int, std::less <const G4VSolid*> > fSolidMap;
};

#endif

#endif
