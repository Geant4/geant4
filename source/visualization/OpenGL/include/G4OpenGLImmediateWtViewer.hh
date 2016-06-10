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
// $Id: G4OpenGLImmediateWtViewer.hh 85263 2014-10-27 08:58:31Z gcosmo $
//
// 
// Class G4OpenGLImmediateWtViewer : a class derived from
//   G4OpenGLWtViewer and G4OpenGLImmediateViewer.

#ifdef G4VIS_BUILD_OPENGLWT_DRIVER

#ifndef G4OPENGLIMMEDIATEWTVIEWER_HH
#define G4OPENGLIMMEDIATEWTVIEWER_HH

#include "G4OpenGLImmediateViewer.hh"
#include "G4OpenGLWtViewer.hh"
#include <Wt/WEvent>

#include "globals.hh"

class G4OpenGLImmediateSceneHandler;

class G4OpenGLImmediateWtViewer :
public G4OpenGLWtViewer, public G4OpenGLImmediateViewer, public Wt::WGLWidget {
   
public:
  G4OpenGLImmediateWtViewer (G4OpenGLImmediateSceneHandler& scene, Wt::WContainerWidget*, const G4String& name = "");
  ~G4OpenGLImmediateWtViewer ();
  void Initialise ();
  void resizeGL(int, int);
  void paintGL();
  void initializeGL ();
  void DrawView();
  void ShowView();
  //
  void popMatrix();
  void pushMatrix();
  void multMatrixd(const GLdouble*);
  void setMatrixUniforms();
  void loadIdentity();
  
  void ComputeView ();
  void FinishView();
  
private:
  //  void showEvent(QShowEvent event );
  void mousePressEvent(Wt::WMouseEvent *event);
  void mouseMoveEvent(Wt::WMouseEvent *event);
  void mouseDoubleClickEvent(Wt::WMouseEvent *event);
  void mouseReleaseEvent(Wt::WMouseEvent event);
  //  void WtContextMenuEvent(QContextMenuEvent e);
  void keyPressEvent (Wt::WKeyEvent *event);
  void paintEvent(Wt::WPaintDevice * event);

  void updateWWidget();
  
  // A client-side JavaScript matrix variable
  JavaScriptMatrix4x4 jsMatrix_;
  
  //  void ComputeView ();
  // implements G4VViewer::SetView() and ClearView()
  //  void SetView ();
  //  void ClearView ();

};

#endif

#endif
