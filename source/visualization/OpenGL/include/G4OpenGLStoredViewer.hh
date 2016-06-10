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
// $Id: G4OpenGLStoredViewer.hh 91268 2015-06-29 07:04:24Z gcosmo $
//
//
// Andrew Walkden  7th February 1997
// Class G4OpenGLStoredViewer : Encapsulates the `storedness' of
//                              an OpenGL viewer, for inheritance by
//                              derived (X, Xm...) classes.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLSTOREDVIEWER_HH
#define G4OPENGLSTOREDVIEWER_HH

#include "G4OpenGLViewer.hh"

class G4OpenGLStoredSceneHandler;
class G4Colour;

class G4OpenGLStoredViewer: virtual public G4OpenGLViewer {
  
public:
  G4OpenGLStoredViewer (G4OpenGLStoredSceneHandler& scene);
  virtual ~G4OpenGLStoredViewer ();
  
protected:
  void KernelVisitDecision ();
  virtual G4bool CompareForKernelVisit(G4ViewParameters&);
  void DrawDisplayLists ();
  virtual void DisplayTimePOColourModification
  (G4Colour&, size_t /*currentPOListIndex*/) {}
  G4OpenGLStoredSceneHandler& fG4OpenGLStoredSceneHandler;
  G4ViewParameters fLastVP;  // Memory for making kernel visit decisions.
  
  // Two virtual functions to return sub-class selection.
  virtual G4bool POSelected(size_t) {return true;}
  virtual G4bool TOSelected(size_t) {return true;}
  
  G4bool fDepthTestEnable;
  G4Colour fOldDisplayListColor;
};

#endif

#endif
