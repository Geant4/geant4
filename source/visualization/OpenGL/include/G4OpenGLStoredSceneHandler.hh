// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLStoredSceneHandler.hh,v 1.2 1999-05-10 14:03:45 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// G4OpenGLStoredSceneHandler - creates OpenGL Display lists.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLSTOREDSCENEHANDLER_HH
#define G4OPENGLSTOREDSCENEHANDLER_HH

#include "G4VSceneHandler.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLStoredViewer.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4OpenGLSceneHandler.hh"

class G4OpenGLStored;

class G4OpenGLStoredSceneHandler: public G4OpenGLSceneHandler {
    
public:
  typedef unsigned long G4VSolidPointer;
  G4OpenGLStoredSceneHandler (G4VGraphicsSystem& system, const G4String& name = "");
  virtual ~G4OpenGLStoredSceneHandler ();
  void BeginPrimitives (const G4Transform3D& objectTransformation);
  void EndPrimitives ();
  void BeginModeling ();
  void EndModeling ();
  static G4int GetSceneCount ();
  G4bool fMemoryForDisplayLists; // avoid memory overflow
  
private:
friend class G4OpenGLStoredViewer;
  // ..allows access to P/TODLs.
  void ClearStore ();
  void RequestPrimitives (const G4VSolid& solid);
  static G4int    fSceneIdCount;   // static counter for OpenGLStored scenes.
  static G4int    fSceneCount;     // No. of extanct scenes.
  G4int           fDisplayListId;  // Workspace.
  
  // PODL = Persistent Object Display List.
  GLint           fTopPODL;       // List which calls the other PODLs.
  RWTValOrderedVector<G4int> fPODLList; 
  RWTValOrderedVector<G4Transform3D> fPODLTransformList; 
  
  // TODL = Transient  Object Display List.
  RWTValOrderedVector<G4int> fTODLList; 
  RWTValOrderedVector<G4Transform3D> fTODLTransformList; 
  
  // Stop-gap solution of structure re-use.
  // A proper implementation would use geometry hierarchy.
  RWTValHashDictionary<G4VSolidPointer, G4int> fSolidDictionary;
};

inline G4int G4OpenGLStoredSceneHandler::GetSceneCount () {
  return fSceneCount;
}

#endif

#endif
