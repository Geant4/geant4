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
// $Id: G4OpenGLViewer.hh 75567 2013-11-04 11:35:11Z gcosmo $
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
class G4OpenGL2PSAction;
class G4Text;
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
#include "G4OpenGLWtDrawer.hh"
#endif

// Base class for various OpenGLView classes.
class G4OpenGLViewer: virtual public G4VViewer {

  friend class G4OpenGLSceneHandler;
  friend class G4OpenGLImmediateSceneHandler;
  friend class G4OpenGLStoredSceneHandler;
  friend class G4OpenGLFileSceneHandler;
  friend class G4OpenGLViewerMessenger;

public:
  void ClearView  ();
//////////////////////////////Vectored PostScript production functions///
  void printEPS();

  // Special case for Wt, we want to have acces to the drawer
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
  inline G4OpenGLWtDrawer* getWtDrawer() {return fWtDrawer;}
  
  // Associate the Wt drawer to the OpenGLViewer and the OpenGLSceneHandler
  void setWtDrawer(G4OpenGLWtDrawer* drawer);
  G4OpenGLWtDrawer* fWtDrawer;
#endif

protected:
  G4OpenGLViewer (G4OpenGLSceneHandler& scene);
  virtual ~G4OpenGLViewer ();

private:
  G4OpenGLViewer(const G4OpenGLViewer&);
  G4OpenGLViewer& operator= (const G4OpenGLViewer&);

protected:
  void SetView    ();
  void ResetView ();

  virtual void DrawText(const G4Text&);
  void ChangePointSize(G4double size);
  void ChangeLineWidth(G4double width);
  void HaloingFirstPass ();
  void HaloingSecondPass ();
  void HLRFirstPass ();
  void HLRSecondPass ();
  void HLRThirdPass ();
  void InitializeGLView ();
  void ResizeGLView();
  void ResizeWindow(unsigned int, unsigned int);
  void Pick(GLdouble x, GLdouble y);
  virtual void CreateFontLists () {}
  void rotateScene (G4double dx, G4double dy);
  void rotateSceneToggle (G4double dx, G4double dy);
//////////////////////////////Vectored PostScript production functions///
  // print EPS file. Depending of fVectoredPs, it will print Vectored or not
  void setPrintSize(G4int,G4int);
  // set the new print size. 
  // -1 means 'print size' = 'window size'
  // Setting size greater than max OpenGL viewport size will set the size to
  // maximum
  void setPrintFilename(G4String name,G4bool inc);
  // set print filename. 
  // if inc, then the filename will be increment by one each time
  std::string getRealPrintFilename();
  unsigned int getWinWidth() const;
  unsigned int getWinHeight() const;
  G4bool sizeHasChanged();
  // return true if size has change since last redraw
  GLdouble getSceneNearWidth();
  GLdouble getSceneFarWidth();
  GLdouble getSceneDepth();
  G4bool isGl2psWriting();
  G4bool                            fPrintColour;
  G4bool                            fVectoredPs;

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
  G4OpenGL2PSAction* fGL2PSAction;

  G4double     fRot_sens;        // Rotation sensibility in degrees
  G4double     fPan_sens;        // Translation sensibility

public :
#ifdef G4OPENGL_VERSION_2
  // define the shaders for vertex and fragment in plain text format
  std::string vertexShader_;
  std::string fragmentShader_;


  // define the keyword shader to handle it in a better way for OpenGL and WebGL
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
  #define Shader Wt::WGLWidget::Shader
#else
  #define Shader GLuint
#endif
  
  // define some attributes and variables for OpenGL and WebGL
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
  Wt::WGLWidget::Program shaderProgram_;

  // Program and related variables
  Wt::WGLWidget::AttribLocation vertexPositionAttribute_;
  Wt::WGLWidget::AttribLocation vertexNormalAttribute_;
  Wt::WGLWidget::UniformLocation pMatrixUniform_;
  Wt::WGLWidget::UniformLocation cMatrixUniform_;
  Wt::WGLWidget::UniformLocation mvMatrixUniform_;
  Wt::WGLWidget::UniformLocation nMatrixUniform_;
  Wt::WGLWidget::UniformLocation tMatrixUniform_;
#else
  GLuint shaderProgram_;

  // Program and related variables
  GLuint vertexPositionAttribute_;
  GLuint vertexNormalAttribute_;
  GLuint pMatrixUniform_;
  GLuint cMatrixUniform_;
  GLuint mvMatrixUniform_;
  GLuint nMatrixUniform_;
  GLuint tMatrixUniform_;
#endif
  
#endif
  
private :
  static G4int                             fPrintSizeX;
  static G4int                             fPrintSizeY;
  static G4String                          fPrintFilename;
  static int                               fPrintFilenameIndex;
  unsigned int fWinSize_x, fWinSize_y;
  G4float                           fPointSize;
  G4bool fSizeHasChanged;
  int fGl2psDefaultLineWith;
  int fGl2psDefaultPointSize;

  // size of the OpenGL frame
  void rotateSceneThetaPhi(G4double dx, G4double dy);
  void rotateSceneInViewDirection (G4double dx, G4double dy);
  bool printGl2PS();
  G4int getRealPrintSizeX();
  G4int getRealPrintSizeY();
  GLubyte* grabPixels (int inColor,
		       unsigned int width,
		       unsigned int height);
  bool printNonVectoredEPS ();
  // print non vectored EPS files

  bool printVectoredEPS();
  // print vectored EPS files
};

#endif

#endif
