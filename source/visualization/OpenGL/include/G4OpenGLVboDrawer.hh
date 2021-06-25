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
//
//
// G4OpenGLVboDrawer : Class to provide Wt and Qt specific
//                       functionality for OpenGL in GEANT4
// All references to Wt removed - 1/3/21 JA

#ifndef G4OpenGLVboDrawer_HH
#define G4OpenGLVboDrawer_HH

#include "G4OpenGL.hh"

#ifdef G4OPENGL_VERSION_2

// GL2PS
#define GL2PS_TEXT_B  4
#define GL2PS_TEXT_BL 5
#define GL2PS_TEXT_BR 6
#define GL2PS_EPS 1
#define GL2PS_PDF 2
#define GL2PS_PS 3
#define GL2PS_SVG 4


#define glEdgeFlag(a) fVboDrawer->empty()
#define glRenderMode(a) fVboDrawer->returnNULL()
#define glClipPlane(a,b) fVboDrawer->empty()
#define glGetIntegerv(a,b) fVboDrawer->empty()
#define glGetFloatv(a,b) fVboDrawer->empty()
#define glGetDoublev(a,b) fVboDrawer->empty()
#define glPassThrough fVboDrawer->empty()
#define glGetBooleanv fVboDrawer->empty()
#define glLoadName(a) fVboDrawer->empty()
#define glPushMatrix() fVboDrawer->empty()
#define glLoadIdentity() fVboDrawer->empty()
#define glPopMatrix() fVboDrawer->empty()
#define glCallList(a) fVboDrawer->empty()
#define glGenLists(a) fVboDrawer->returnNULL()
#define glVertex3d fVboDrawer->empty()
#define glBegin fVboDrawer->empty()
#define glEnd fVboDrawer->empty()
#define glNewList(a,b) fVboDrawer->empty()
#define glEndList() fVboDrawer->empty()
#define glPolygonMode(a,b) fVboDrawer->empty()
#define glDrawBuffer(a) fVboDrawer->empty()
#define glDeleteLists(a,b) fVboDrawer->empty()
#define glStencilFunc(a,b,c) fVboDrawer->empty()
#define glStencilOp(a,b,c) fVboDrawer->empty()
#define glColorMaterial(a,b) fVboDrawer->empty()
#define glLightfv(a,b,c) fVboDrawer->empty()
#define glScaled(a,b,c) fVboDrawer->empty()
#define gluLookAt fVboDrawer->empty()
#define gluPickMatrix fVboDrawer->empty()
#define glSelectBuffer(a,b) fVboDrawer->empty()
#define glInitNames() fVboDrawer->empty()
#define glPushNames(a) fVboDrawer->empty()
#define glPushName(a) fVboDrawer->empty()
#define glPixelStorei(a,b) fVboDrawer->empty()
#define glRasterPos3d(a,b,c) fVboDrawer->empty()
#define Geant4_gl2psTextOpt(a,b,c,d,e) fVboDrawer->empty()
#define glMaterialfv(a,b,c)  fVboDrawer->empty()
#define glCullFace(a) fVboDrawer->empty()
#define glReadBuffer(a) fVboDrawer->empty()
#define glReadPixels(a,b,c,d,e,f,g) fVboDrawer->empty()
#define glTranslatef(a,b,c) fVboDrawer->empty() // TO BE FIXED

// +--------------------------------+
// +        QT (OpenGL ES) case     +
// +--------------------------------+

class G4OpenGLImmediateQtViewer;

#define glOrtho fVboDrawer->vboGlOrtho
#define glFrustum fVboDrawer->vboGlFrustum
#define glMultMatrixf fVboDrawer->vboGlMultMatrixf
#define glMultMatrixd fVboDrawer->vboGlMultMatrixd
#define glMatrixMode fVboDrawer->vboGlMatrixMode
#define glPointSize fVboDrawer->vboGlPointSize
#define glColor3d fVboDrawer->vboGlColor3d
#define glColor4d fVboDrawer->vboGlColor4d
#define glColor4fv fVboDrawer->vboGlColor4fv


class G4OpenGLViewer;

class G4OpenGLVboDrawer {
public:
  G4OpenGLVboDrawer (G4OpenGLViewer*, std::string type);
  // Create a new OpenGL Drawer. Type could be one of the following :
  // OGL-ES, OGL-Stored, OGL-Immediate, OGL-VBO
  
  virtual ~G4OpenGLVboDrawer ();

  void vboGlMultMatrixf( const GLfloat *m );
  void vboGlMultMatrixd( const GLdouble *m );

  void vboGlFlush();
  void vboGlOrtho(GLdouble, GLdouble, GLdouble, GLdouble, GLdouble, GLdouble);
  void vboGlFrustum(GLdouble, GLdouble, GLdouble, GLdouble, GLdouble, GLdouble);
  void vboGlMatrixMode(GLenum);
  void vboGlPointSize(float size);
  inline void vboGlColor3d(int red,int green,int blue) {
    vboGlColor4d(red,green, blue, 1.0);
  }
  void vboGlColor4d(int red,int green,int blue,int alpha);
  void vboGlColor4fv(const GLfloat*);
  inline const char * getFragmentShaderSrc() {
    return fFragmentShaderSrc;
  }
  inline const char * getVertexShaderSrc() {
    return fVertexShaderSrc;
  }
  inline bool isVBO() {
    if ((fOGLType == "OGL-ES") || (fOGLType == "OGL-VBO")) {
      return true;
    }
    return false;
  }
  inline void empty() {}
  inline GLuint returnNULL() {
    return 0;
  }

private:
  const char *fFragmentShaderSrc;
  const char *fVertexShaderSrc;
  std::string fOGLType;
  GLenum fMatrixMode;
  
  G4OpenGLImmediateQtViewer* fVboViewer;
};

#endif // G4OPENGL_VERSION_2

#endif // G4OpenGLVboDrawer_HH

