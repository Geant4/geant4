// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLStoredScene.hh,v 1.1 1999-01-07 16:14:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// G4OpenGLStoredScene - creates OpenGL Display lists.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLSTOREDSCENE_HH
#define G4OPENGLSTOREDSCENE_HH

#include "G4VScene.hh"
#include "G4OpenGLView.hh"
#include "G4OpenGLStoredView.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4OpenGLScene.hh"

class G4OpenGLStored;

class G4OpenGLStoredScene: public G4OpenGLScene {
    
public:
  typedef unsigned long G4VSolidPointer;
  G4OpenGLStoredScene (G4VGraphicsSystem& system, const G4String& name = "");
  ~G4OpenGLStoredScene ();
  void BeginPrimitives (const G4Transform3D& objectTransformation);
  void EndPrimitives ();
  void BeginModeling ();
  void EndModeling ();
  static G4int GetSceneCount ();
  G4bool fMemoryForDisplayLists; // avoid memory overflow
  
private:
friend class G4OpenGLStoredView;
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

inline G4int G4OpenGLStoredScene::GetSceneCount () {
  return fSceneCount;
}

#endif

#endif
