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
// $Id$
//
// 
// Class G4OpenGLImmediateWtViewer : a class derived from
//   G4OpenGLWtViewer and G4OpenGLImmediateViewer.

#ifdef G4VIS_BUILD_OPENGLWT_DRIVER

#ifndef G4OPENGLIMMEDIATEWTVIEWER_HH
#define G4OPENGLIMMEDIATEWTVIEWER_HH

#include "G4OpenGLImmediateViewer.hh"
#include "G4OpenGLImmediateQtViewer.hh"
#include "G4OpenGLWtViewer.hh"
#include <Wt/WPaintedWidget>
#include <Wt/WEvent>

//#include <qgl.h> // include <qglwidget.h>

#include "globals.hh"

class G4OpenGLImmediateSceneHandler;
class QMouseEvent;

class G4OpenGLImmediateWtViewer : public G4VViewer,
  /*  public G4OpenGLImmediateViewer ,*/ public Wt::WPaintedWidget {
   
public:
  G4OpenGLImmediateWtViewer (G4OpenGLImmediateSceneHandler& scene, Wt::WContainerWidget *,
                const G4String& name = "");
  ~G4OpenGLImmediateWtViewer ();
  void DrawView();
  void ShowView();

private:
  G4OpenGLImmediateQtViewer * fQtViewer;
  //  void WtShowEvent(QShowEvent event );
  void WtWheelEvent(Wt::WMouseEvent event);
  void WtMousePressEvent(Wt::WMouseEvent event);
  void WtMouseMoveEvent(Wt::WMouseEvent event);
  void WtMouseDoubleClickEvent(Wt::WMouseEvent event);
  //  void WtMouseReleaseEvent(Wt::WMouseEvent event);
  //  void WtContextMenuEvent(QContextMenuEvent e);
  void WtKeyPressEvent (Wt::WKeyEvent event); 
  void paintEvent(Wt::WPaintDevice * event);

  QMouseEvent * ConvertWtMouseEventToQt(Wt::WMouseEvent event);
  QWheelEvent * ConvertWtWheelEventToQt(Wt::WMouseEvent event);
  QKeyEvent * ConvertWtKeyEventToQt(Wt::WKeyEvent event);

  // implements G4VViewer::SetView() and ClearView()
  void SetView ();
  void ClearView ();
};

#endif

#endif
