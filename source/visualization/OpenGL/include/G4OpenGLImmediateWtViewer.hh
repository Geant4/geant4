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
  void SetView();
  //
  void popMatrix();
  void pushMatrix();
  void multMatrixd(const GLdouble*);
  void setMatrixUniforms();
  void loadIdentity();
  
  void wtDrawArrays(GLenum mode,int first, G4int nPoints, std::vector<double> a_vertices);
  void enableClientState(int mode);
  void disableClientState(int mode);
  
  void ComputeView ();
  void drawScene ();
  void FinishView();
  
  Program shaderProgram_;

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
  void centerpoint(double &x, double &y, double &z);

  // Program and related variables
  AttribLocation vertexPositionAttribute_;
  AttribLocation vertexNormalAttribute_;
  UniformLocation pMatrixUniform_;
  UniformLocation cMatrixUniform_;
  UniformLocation mvMatrixUniform_;
  UniformLocation nMatrixUniform_;
  
  // A client-side JavaScript matrix variable
  JavaScriptMatrix4x4 jsMatrix_;
  
  // The so-called VBOs, Vertex Buffer Objects
  // This one contains both vertex (xyz) and normal (xyz) data
  Buffer objBuffer_;
  //  void ComputeView ();
  // implements G4VViewer::SetView() and ClearView()
  //  void SetView ();
  //  void ClearView ();

  // The shaders, in plain text format
  std::string vertexShader_;
  std::string fragmentShader_;
  
  
  // To avoid copying large constant data around, the data points are stored
  // in a global variable.
  std::vector<double> data;
  GLfloat* normals;
  std::vector <Buffer> VBO_Buffer;
#ifdef TEST_WT_EXAMPLE
  // Sets the shader source. Must be set before the widget is first rendered.
  void setShaders(const std::string &vertexShader,
                  const std::string &fragmentShader);
  void readObj(const std::string &fname,
               std::vector<double> &data);
  
  
protected:
  void drawCube();
private:
  
#endif


};

#endif

#endif
