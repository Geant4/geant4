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
// $Id: G4OpenGLXViewer.hh 87695 2014-12-17 09:35:24Z gcosmo $
//
// 
// Andrew Walkden  7th February 1997
// G4OpenGLXViewer : Class to provide XWindows specific
//                   functionality for OpenGL in GEANT4

#ifdef G4VIS_BUILD_OPENGLX_DRIVER

#ifndef G4OPENGLXVIEWER_HH
#define G4OPENGLXVIEWER_HH

#include "G4OpenGLViewer.hh"
#include "globals.hh"

#include <X11/Intrinsic.h>

#include "G4OpenGL.hh"
#include <GL/glx.h>

class G4OpenGLSceneHandler;
class G4Text;

class G4OpenGLXViewer: virtual public G4OpenGLViewer {

  friend class G4OpenGLXViewerMessenger;
  friend class G4OpenGLXmViewer;

public:
  G4OpenGLXViewer (G4OpenGLSceneHandler& scene);
  virtual ~G4OpenGLXViewer ();
  void SetView ();
  void ShowView ();
#ifdef G4MULTITHREADED
  void SwitchToVisSubThread();
  void SwitchToMasterThread();
#endif
  void DrawText(const G4Text&);

protected:
  void GetXConnection ();
  void CreateGLXContext (XVisualInfo* vi);
  virtual void CreateMainWindow ();
  virtual void CreateFontLists ();

  static int snglBuf_RGBA[12];
  static int dblBuf_RGBA[13];

//////////////////////////////Pixmap (screen dump) production functions/////

  XWindowAttributes                 xwa;
  Display                           *dpy;
  static XVisualInfo                *vi_single_buffer;
  static XVisualInfo                *vi_double_buffer;
  XVisualInfo                       *vi_immediate,
                                    *vi_stored,
                                    *vi;
  Colormap                          cmap;
  XSetWindowAttributes              swa;
  GLXDrawable                       win;
  GLXContext                        cxMaster;
#ifdef G4MULTITHREADED
  GLXContext                        cxVisSubThread;
#endif
  XEvent                            event;
  G4int                             *attributeList,
                                    errorBase,
                                    eventBase,
                                    major,
                                    minor;
  XSizeHints                        *norm_hints;
  XWMHints                          *wm_hints;
  XClassHint                        *class_hints;
  Pixmap                            icon_pixmap;
  XSizeHints                        *size_hints;
  Atom                              Xatom;
  XTextProperty                     windowName,
                                    iconName;
  char                              charViewName [100];


private:
  G4OpenGLXViewer (const G4OpenGLXViewer&);
  G4OpenGLXViewer& operator = (const G4OpenGLXViewer&);
  GLXContext                        tmp_cx;
};

#endif

#endif
