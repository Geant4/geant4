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
// $Id: G4OpenGLViewer.hh,v 1.20 2007/05/08 11:04:11 allison Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// 
// Andrew Walkden  27th March 1996
// OpenGL viewer - opens window, hard copy, etc.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLVIEWER_HH
#define G4OPENGLVIEWER_HH

#include "G4VViewer.hh"
#include "G4OpenGL.hh"

class G4OpenGLSceneHandler;

// Base class for various OpenGLView classes.
class G4OpenGLViewer: virtual public G4VViewer {

  friend class G4OpenGLSceneHandler;
  friend class G4OpenGLImmediateSceneHandler;
  friend class G4OpenGLStoredSceneHandler;
  friend class G4OpenGLViewerMessenger;

public:
  void ClearView  ();

protected:
  G4OpenGLViewer (G4OpenGLSceneHandler& scene);
  virtual ~G4OpenGLViewer ();
  void SetView    ();
  void HaloingFirstPass ();
  void HaloingSecondPass ();
  void HLRFirstPass ();
  void HLRSecondPass ();
  void HLRThirdPass ();
  void InitializeGLView ();
  void Pick(GLdouble x, GLdouble y);
  virtual void CreateFontLists () {}
  virtual void print();
//////////////////////////////Vectored PostScript production functions///
  void printBuffer(GLint, GLfloat*);
  GLfloat* spewPrimitiveEPS (FILE*, GLfloat*);
  void spewSortedFeedback (FILE*, GLint, GLfloat*);
  void spewWireframeEPS (FILE*, GLint, GLfloat*, const char*);
  void print3DcolorVertex(GLint, GLint*, GLfloat*);
  G4float                           pointSize;
  char                              print_string[50];
  G4bool                            print_colour;
  G4bool                            vectored_ps;

  G4OpenGLSceneHandler& fOpenGLSceneHandler;
  G4Colour background;      //the OpenGL clear colour
  G4bool
    transparency_enabled,   //is alpha blending enabled?
    antialiasing_enabled,   //is antialiasing enabled?
    haloing_enabled;        //is haloing enabled for wireframe?
  G4double fStartTime, fEndTime;  // Time range (e.g., for trajectory steps).
  G4double fFadeFactor;  // 0: no fade; 1: maximum fade with time within range.
  G4bool fDisplayHeadTime;  // Display head time (fEndTime) in 2D text.
  G4double fDisplayHeadTimeX, fDisplayHeadTimeY;  // 2D screen coords.
  G4double fDisplayHeadTimeSize;  // Screen size.
  G4double fDisplayHeadTimeRed, fDisplayHeadTimeGreen, fDisplayHeadTimeBlue;
  G4bool fDisplayLightFront;// Display light front at head time originating at
  G4double fDisplayLightFrontX, fDisplayLightFrontY, fDisplayLightFrontZ,
    fDisplayLightFrontT;
  G4double fDisplayLightFrontRed, fDisplayLightFrontGreen, fDisplayLightFrontBlue;
};

typedef struct G4OpenGLViewerFeedback3Dcolor {
  GLfloat x;
  GLfloat y;
  GLfloat z;
  GLfloat red;
  GLfloat green;
  GLfloat blue;
  GLfloat alpha;
} Feedback3Dcolor;

#endif

#endif
