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
// $Id: G4OpenGLViewer.hh 103926 2017-05-03 13:43:27Z gcosmo $
//
// 
// Andrew Walkden  27th March 1996
// OpenGL viewer - opens window, hard copy, etc.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLVIEWER_HH
#define G4OPENGLVIEWER_HH

#include "G4VViewer.hh"
#include "G4OpenGL.hh"
#ifdef G4OPENGL_VERSION_2
#include "G4OpenGLVboDrawer.hh"
#endif

class G4OpenGLSceneHandler;
class G4OpenGL2PSAction;
class G4Text;

class G4OpenGLViewerPickMap {
  public :
  inline void setName(G4String n) {
    fName = n;
  }
  
  inline void setHitNumber(G4int n) {
    fHitNumber = n;
  }

  inline void setSubHitNumber(G4int n) {
    fSubHitNumber = n;
  }
  inline void setPickName(G4int n) {
    fPickName= n;
  }

  inline void addAttributes(G4String att) {
    fAttributes.push_back(att);
  }


  inline G4String getName() {
    return fName;
  }
  inline G4int getHitNumber() {
    return fHitNumber;
  }

  inline G4int getSubHitNumber() {
    return fSubHitNumber;
  }

  inline G4int getPickName() {
    return fPickName;
  }

  inline std::vector <G4String > getAttributes() {
    return fAttributes;
  }

  G4String print();
  
  private :
  G4String fName;
  G4int fHitNumber;
  G4int fSubHitNumber;
  G4int fPickName;
  std::vector <G4String > fAttributes;

};

// Base class for various OpenGLView classes.
class G4OpenGLViewer: virtual public G4VViewer {

  friend class G4OpenGLSceneHandler;
  friend class G4OpenGLImmediateSceneHandler;
  friend class G4OpenGLStoredSceneHandler;
  friend class G4OpenGLFileSceneHandler;
  friend class G4OpenGLViewerMessenger;

public:
  void ClearView  ();
  void ClearViewWithoutFlush ();
//////////////////////////////Vectored PostScript production functions///
  bool printEPS();
  virtual bool exportImage(std::string name="", int width=-1, int height=-1);

  bool setExportImageFormat(std::string format,bool quiet = false);
  // change the export image format according to thoses available for the current viewer

  // Special case for Wt, we want to have acces to the drawer
#ifdef G4OPENGL_VERSION_2
  inline G4OpenGLVboDrawer* getWtDrawer() {return fVboDrawer;}
  
  // Associate the Wt drawer to the OpenGLViewer and the OpenGLSceneHandler
  void setVboDrawer(G4OpenGLVboDrawer* drawer);
  G4OpenGLVboDrawer* fVboDrawer;

  inline bool isInitialized() {
    return fGlViewInitialized;
  }
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
  virtual G4String Pick(GLdouble x, GLdouble y);
  const std::vector < G4OpenGLViewerPickMap* > & GetPickDetails(GLdouble x, GLdouble y);
  virtual void CreateFontLists () {}
  void rotateScene (G4double dx, G4double dy);
  void rotateSceneToggle (G4double dx, G4double dy);
//////////////////////////////Vectored PostScript production functions///
  // print EPS file. Depending of fVectoredPs, it will print Vectored or not
  void setExportSize(G4int,G4int);
  // set the new print size. 
  // -1 means 'print size' = 'window size'
  // Setting size greater than max OpenGL viewport size will set the size to
  // maximum
  bool setExportFilename(G4String name,G4bool inc = true);
  // set export filename.
  // if inc, then the filename will be increment by one each time
  // try to guesss the correct format according to the extention

  std::string getRealPrintFilename();
  unsigned int getWinWidth() const;
  unsigned int getWinHeight() const;
  G4bool sizeHasChanged();
  // return true if size has change since last redraw
  GLdouble getSceneNearWidth();
  GLdouble getSceneFarWidth();
  GLdouble getSceneDepth();
  void addExportImageFormat(std::string format);
  // add a image format to the available export format list
  G4bool isGl2psWriting();
  G4bool isFramebufferReady();
  
  void g4GluPickMatrix(GLdouble x, GLdouble y, GLdouble width, GLdouble height,
                       GLint viewport[4]);
  // MESA implementation of gluPickMatrix
  
  void g4GluLookAt( GLdouble eyex, GLdouble eyey, GLdouble eyez,
                   GLdouble centerx, GLdouble centery, GLdouble
                   centerz,
                   GLdouble upx, GLdouble upy, GLdouble upz );
  // MESA implementation of gluLookAt
  void g4GlOrtho (GLdouble left, GLdouble right,
                  GLdouble bottom, GLdouble top,
                  GLdouble near, GLdouble far);
  // Redefinition on glOrtho to solve precision issues
  void g4GlFrustum (GLdouble left, GLdouble right,
                  GLdouble bottom, GLdouble top,
                  GLdouble near, GLdouble far);
  // Redefinition on glFrustum to solve precision issues

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
  unsigned int fWinSize_x;
  unsigned int fWinSize_y;
  std::vector < std::string > fExportImageFormatVector;
  std::string fDefaultExportImageFormat;
  std::string fExportImageFormat;
  int fExportFilenameIndex;
  G4int fPrintSizeX;
  G4int fPrintSizeY;


private :
  G4float fPointSize;
  G4String fExportFilename;
  G4String fDefaultExportFilename;
  G4bool fSizeHasChanged;
  int fGl2psDefaultLineWith;
  int fGl2psDefaultPointSize;
  bool fGlViewInitialized;
  
  // size of the OpenGL frame
  void rotateSceneThetaPhi(G4double dx, G4double dy);
  void rotateSceneInViewDirection (G4double dx, G4double dy);
  bool printGl2PS();
  G4int getRealExportWidth();
  G4int getRealExportHeight();
  GLubyte* grabPixels (int inColor,
		       unsigned int width,
		       unsigned int height);
  bool printNonVectoredEPS ();
  // print non vectored EPS files

  bool printVectoredEPS();
  // print vectored EPS files
  
  bool fIsGettingPickInfos;
  // Block SetView() during picking
  
#ifdef G4OPENGL_VERSION_2
public:
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
  inline Wt::WGLWidget::Program getShaderProgram() {
    return fShaderProgram;
  }
  inline Wt::WGLWidget::UniformLocation getShaderProjectionMatrix() {
    return fpMatrixUniform;
  }
  inline Wt::WGLWidget::UniformLocation getShaderTransformMatrix() {
    return ftMatrixUniform;
  }
#else
  inline GLuint getShaderProgram() {
    return fShaderProgram;
  }
  inline GLuint getShaderProjectionMatrix() {
    return fpMatrixUniform;
  }
  inline GLuint getShaderTransformMatrix() {
    return ftMatrixUniform;
  }
  inline GLuint getShaderViewModelMatrix() {
    return fmvMatrixUniform;
  }
#endif

protected :
  
  // define the keyword shader to handle it in a better way for OpenGL and WebGL
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
#define Shader Wt::WGLWidget::Shader
#else
#define Shader GLuint
#endif
  
  // define some attributes and variables for OpenGL and WebGL
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
  Wt::WGLWidget::Program fShaderProgram;
  
  // Program and related variables
  Wt::WGLWidget::AttribLocation fVertexPositionAttribute;
  Wt::WGLWidget::AttribLocation fVertexNormalAttribute;
  Wt::WGLWidget::UniformLocation fpMatrixUniform;
  Wt::WGLWidget::UniformLocation fcMatrixUniform;
  Wt::WGLWidget::UniformLocation fmvMatrixUniform;
  Wt::WGLWidget::UniformLocation fnMatrixUniform;
  Wt::WGLWidget::UniformLocation ftMatrixUniform;
#else
  GLuint fShaderProgram;
  
  // Program and related variables
  GLuint fVertexPositionAttribute;
  GLuint fVertexNormalAttribute;
  GLuint fpMatrixUniform;
  GLuint fcMatrixUniform;
  GLuint fmvMatrixUniform;
  GLuint fnMatrixUniform;
  GLuint ftMatrixUniform;
#endif
  
#endif
};

#endif

#endif
