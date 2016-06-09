//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4OpenGLStoredSceneHandler.hh,v 1.14 2003/06/16 17:13:38 gunter Exp $
// GEANT4 tag $Name: geant4-06-00 $
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
#include <map>
#include <vector>

class G4OpenGLStored;

class G4OpenGLStoredSceneHandler: public G4OpenGLSceneHandler {
    
public:
  G4OpenGLStoredSceneHandler (G4VGraphicsSystem& system, const G4String& name = "");
  virtual ~G4OpenGLStoredSceneHandler ();
  void BeginPrimitives (const G4Transform3D& objectTransformation);
  void EndPrimitives ();
  void BeginModeling ();
  void EndModeling ();
  static G4int GetSceneCount ();
private:
  friend class G4OpenGLStoredViewer;
  // ..allows access to P/TODLs.
  void ClearStore ();
  void ClearTransientStore ();
  void RequestPrimitives (const G4VSolid& solid);
  static G4int     fSceneIdCount;   // static counter for OpenGLStored scenes.
  static G4int     fSceneCount;     // No. of extanct scene handlers.
  G4int            fDisplayListId;  // Workspace.
  G4bool  fMemoryForDisplayLists;   // avoid memory overflow
  G4ViewParameters fLastVP; // Keeps memory of last view parameters,
         // namely the ones which go with the graphical database
         // (display lists).  A view may interrogate this to decide if
         // it needs to visit the kernel, i.e., to remake the
         // graphical database.
  
  // PODL = Persistent Object Display List.
  GLint           fTopPODL;       // List which calls the other PODLs.
  std::vector<G4int> fPODLList; 
  std::vector<G4Transform3D> fPODLTransformList; 
  
  // TODL = Transient  Object Display List.
  std::vector<G4int> fTODLList; 
  std::vector<G4Transform3D> fTODLTransformList; 
  
  // Stop-gap solution of structure re-use.
  // A proper implementation would use geometry hierarchy.
  std::map <const G4VSolid*, G4int, std::less <const G4VSolid*> > fSolidMap;
};

inline G4int G4OpenGLStoredSceneHandler::GetSceneCount () {
  return fSceneCount;
}

#endif

#endif
