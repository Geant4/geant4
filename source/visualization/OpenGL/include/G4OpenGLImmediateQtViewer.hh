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
// $Id: G4OpenGLImmediateQtViewer.hh 83403 2014-08-21 15:07:30Z gcosmo $
//
// 
// Class G4OpenGLImmediateQtViewer : a class derived from
//   G4OpenGLQtViewer and G4OpenGLImmediateViewer.

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#ifndef G4OPENGLIMMEDIATEQTVIEWER_HH
#define G4OPENGLIMMEDIATEQTVIEWER_HH

#include "G4OpenGLImmediateViewer.hh"
#include "G4OpenGLQtViewer.hh"

#include <qgl.h> // include <qglwidget.h>

#include "globals.hh"

class G4OpenGLImmediateSceneHandler;

class G4OpenGLImmediateQtViewer:
  public G4OpenGLQtViewer, public G4OpenGLImmediateViewer, public QGLWidget {
  
public:
  G4OpenGLImmediateQtViewer (G4OpenGLImmediateSceneHandler& scene,
                const G4String& name = "");
  ~G4OpenGLImmediateQtViewer ();
  void Initialise ();
  void initializeGL ();
  void DrawView ();
  void resizeGL(int width,int height);
  void paintGL();
  void updateQWidget();
  void ShowView ();

protected:
  void showEvent(QShowEvent * event );
  void wheelEvent(QWheelEvent *event);
  void mousePressEvent(QMouseEvent *event);
  void mouseMoveEvent(QMouseEvent *event);
  void mouseDoubleClickEvent(QMouseEvent *event);
  void mouseReleaseEvent(QMouseEvent *event);
  void contextMenuEvent(QContextMenuEvent *e);
  void keyPressEvent (QKeyEvent * event); 
  void keyReleaseEvent (QKeyEvent * event);
  void paintEvent(QPaintEvent *event);
private:
  void ComputeView ();

};

#endif

#endif
