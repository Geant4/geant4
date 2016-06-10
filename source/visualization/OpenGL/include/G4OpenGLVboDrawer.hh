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
// $Id: G4OpenGLVboDrawer.hh 74103 2014-06-23 07:52:38Z lgarnier $
//
//
// G4OpenGLVboDrawer : Class to provide Wt and Qt specific
//                       functionality for OpenGL in GEANT4

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
// +        WT (OpenGL ES) case     +
// +--------------------------------+

#ifdef G4VIS_BUILD_OPENGLWT_DRIVER

class G4OpenGLImmediateWtViewer;

// specific definition for WT :
// WARNING fVboDrawer should be the exact name of the object!

#define glGetError() Wt::WGLWidget::NONE

#define glOrtho fVboDrawer->vboGlOrtho
#define glFrustum fVboDrawer->vboGlFrustum
#define glViewport fVboDrawer->vboGlViewport
#define glEnable fVboDrawer->vboGlEnable
#define glDisable fVboDrawer->vboGlDisable
#define glBlendFunc fVboDrawer->vboGlBlendFunc
#define glClear fVboDrawer->vboGlClear
#define glClearColor fVboDrawer->vboGlClearColor
#define glClearDepth fVboDrawer->vboGlClearDepth
#define glDepthFunc fVboDrawer->vboGlDepthFunc
#define glDepthMask fVboDrawer->vboGlDepthMask
#define glFlush fVboDrawer->vboGlFlush
#define glColorMask fVboDrawer->vboGlColorMask
#define glLineWidth fVboDrawer->vboGlLineWidth
#define glUniformMatrix4 fVboDrawer->vboGlUniformMatrix4
#define glDrawArrays fVboDrawer->vboGlDrawArrays
#define glCreateBuffer fVboDrawer->vboGlCreateBuffer
#define glVertexPointer fVboDrawer->vboGlVertexPointer
#define glBindBuffer fVboDrawer->vboGlBindBuffer
#define glDeleteBuffer fVboDrawer->vboGlDeleteBuffer
#define glBufferDatafv  fVboDrawer->vboGlBufferDatafv
#define glBufferDataiv  fVboDrawer->vboGlBufferDataiv
#define glGetAttribLocation fVboDrawer->vboGlGetAttribLocation
#define glEnableVertexAttribArray fVboDrawer->vboGlEnableVertexAttribArray
#define glDisableVertexAttribArray fVboDrawer->vboGlDisableVertexAttribArray
#define glShaderSource fVboDrawer->vboGlShaderSource
#define glCompileShader fVboDrawer->vboGlCompileShader
#define glCreateShader fVboDrawer->vboGlCreateShader
#define glCreateProgram fVboDrawer->vboGlCreateProgram
#define glAttachShader fVboDrawer->vboGlAttachShader
#define glLinkProgram fVboDrawer->vboGlLinkProgram
#define glUseProgram fVboDrawer->vboGlUseProgram
#define glDrawElements fVboDrawer->vboGlDrawElements
#define glVertexAttribPointer fVboDrawer->vboGlVertexAttribPointer
#define glGetUniformLocation fVboDrawer->vboGlGetUniformLocation
#define glPointSize fVboDrawer->vboGlPointSize
#define glColor3d fVboDrawer->vboGlColor3d
#define glColor4d fVboDrawer->vboGlColor4d
#define glColor4fv fVboDrawer->vboGlColor4fv
#define glMultMatrixd fVboDrawer->vboGlMultMatrixd
#define glMultMatrixf fVboDrawer->vboGlMultMatrixf
#define glGetUniformLocation fVboDrawer->vboGlGetUniformLocation
#define glGetAttribLocation fVboDrawer->vboGlGetAttribLocation
#define glMatrixMode fVboDrawer->vboGlMatrixMode


// Only used in fvboDrawer->VboDrawer to be compatible between Wt and Qt framework
#define glUniform1f fVboViewer->uniform1f
#define glUniform4fv fVboViewer->uniform4fv
#define glUniformMatrix4dv fVboDrawer->vboGlUniformMatrix4;
#define glUniformMatrix4fv fVboDrawer->vboGlUniformMatrix4fv;



#define GL_VIEWPORT Wt::WGLWidget::VIEWPORT
#define GL_RGBA Wt::WGLWidget::RGBA
#define GL_ONE_MINUS_SRC_ALPHA Wt::WGLWidget::ONE_MINUS_SRC_ALPHA
#define GL_BLEND Wt::WGLWidget::BLEND
#define GL_SRC_ALPHA Wt::WGLWidget::SRC_ALPHA
#define GL_LEQUAL Wt::WGLWidget::LEQUAL
#define GL_FALSE false
#define GL_LESS Wt::WGLWidget::LESS
#define GL_SELECT Wt::WGLWidget::SELECT
#define GL_TRUE true
#define GL_RGB Wt::WGLWidget::RGB
#define GL_CURRENT_RASTER_POSITION_VALID Wt::WGLWidget::CURRENT_RASTER_POSITION_VALID
#define GL_ONE Wt::WGLWidget::ONE
#define GL_ZERO Wt::WGLWidget::ZERO
#define GL_COLOR_INDEX Wt::WGLWidget::COLOR_INDEX
#define GL_LINE_TOKEN Wt::WGLWidget::LINE_TOKEN
#define GL_LINE_RESET_TOKEN Wt::WGLWidget::LINE_RESET_TOKEN
#define GL_POLYGON_TOKEN Wt::WGLWidget::POLYGON_TOKEN
#define GL_FEEDBACK Wt::WGLWidget::FEEDBACK
#define GL_COLOR_CLEAR_VALUE Wt::WGLWidget::COLOR_CLEAR_VALUE
#define GL_BITMAP_TOKEN Wt::WGLWidget::BITMAP_TOKEN
#define GL_DRAW_PIXEL_TOKEN Wt::WGLWidget::DRAW_PIXEL_TOKEN
#define GL_COPY_PIXEL_TOKEN Wt::WGLWidget::COPY_PIXEL_TOKEN
#define GL_PASS_THROUGH_TOKEN Wt::WGLWidget::PASS_THROUGH_TOKEN
#define GL_3D_COLOR Wt::WGLWidget::3D_COLOR
#define GL_DEPTH_TEST Wt::WGLWidget::DEPTH_TEST
#define GL_FRONT Wt::WGLWidget::FRONT
#define GL_BACK Wt::WGLWidget::BACK
#define GL_FRONT_AND_BACK Wt::WGLWidget::FRONT_AND_BACK
#define GL_OUT_OF_MEMORY Wt::WGLWidget::OUT_OF_MEMORY
#define GL_LINE_STRIP Wt::WGLWidget::LINE_STRIP
#define GL_QUADS Wt::WGLWidget::QUADS
#define GL_LINE_LOOP Wt::WGLWidget::LINE_LOOP
#define GL_LINES Wt::WGLWidget::LINES
#define GL_POINTS Wt::WGLWidget::POINTS
#define GL_TRIANGLES Wt::WGLWidget::TRIANGLES
#define GL_TRIANGLE_STRIP Wt::WGLWidget::TRIANGLE_STRIP
#define GL_TRIANGLE_FAN Wt::WGLWidget::TRIANGLE_FAN
#define GL_FLOAT Wt::WGLWidget::FLOAT
#define GL_STENCIL_TEST Wt::WGLWidget::STENCIL_TEST
#define GL_ALWAYS Wt::WGLWidget::ALWAYS
#define GL_INVERT Wt::WGLWidget::INVERT
#define GL_COMPILE_AND_EXECUTE Wt::WGLWidget::COMPILE_AND_EXECUTE
#define GL_COMPILE Wt::WGLWidget::COMPILE
#define GL_COLOR_BUFFER_BIT Wt::WGLWidget::COLOR_BUFFER_BIT
#define GL_DEPTH_BUFFER_BIT Wt::WGLWidget::DEPTH_BUFFER_BIT
#define GL_STENCIL_BUFFER_BIT Wt::WGLWidget::STENCIL_BUFFER_BIT
#define GL_UNSIGNED_BYTE Wt::WGLWidget::UNSIGNED_BYTE
#define GL_ARRAY_BUFFER Wt::WGLWidget::ARRAY_BUFFER
#define GL_ELEMENT_ARRAY_BUFFER Wt::WGLWidget::ELEMENT_ARRAY_BUFFER
#define GL_RENDER Wt::WGLWidget::RENDER
#define GL_LUMINANCE Wt::WGLWidget::LUMINANCE
#define GL_STATIC_DRAW Wt::WGLWidget::STATIC_DRAW
#define GL_FRAGMENT_SHADER Wt::WGLWidget::FRAGMENT_SHADER
#define GL_VERTEX_SHADER Wt::WGLWidget::VERTEX_SHADER
#define GL_UNSIGNED_INT Wt::WGLWidget::UNSIGNED_INT
#define GL_UNSIGNED_SHORT Wt::WGLWidget::UNSIGNED_SHORT
#define GL_CULL_FACE Wt::WGLWidget::CULL_FACE
#define GL_MAX_VIEWPORT_DIMS Wt::WGLWidget::MAX_VIEWPORT_DIMS
#define GL_PROJECTION Wt::WGLWidget::FRAGMENT_SHADER  // Not the good value, but should be ok, work together with GL_MODELVIEW
#define GL_MODELVIEW Wt::WGLWidget::VERTEX_SHADER // Not the good value, but should be ok, work together with GL_PROJECTION

// to be implemented
#define GL_LINE 0
#define GL_FILL 0
#define GL_PROJECTION_MATRIX 0
#define GL_UNPACK_SWAP_BYTES 0
#define GL_UNPACK_LSB_FIRST 0
#define GL_UNPACK_SKIP_ROWS 0
#define GL_UNPACK_LOW_LENGHT 0
#define GL_UNPACK_SKIP_PIXELS 0
#define GL_UNPACK_ALIGNMENT 0
#define GL_UNPACK_ROW_LENGTH 0
#define GL_CLIP_PLANE0 Wt::WGLWidget::NONE
#define GL_CLIP_PLANE1 Wt::WGLWidget::NONE
#define GL_CLIP_PLANE2 Wt::WGLWidget::NONE
#define GL_CLIP_PLANE3 Wt::WGLWidget::NONE
#define GL_CLIP_PLANE4 Wt::WGLWidget::NONE
#define GL_COLOR_MATERIAL Wt::WGLWidget::NONE
#define GL_AMBIENT_AND_DIFFUSE Wt::WGLWidget::NONE
#define GL_POLYGON 0
#define GL_LIGHTING Wt::WGLWidget::NONE
#define GL_POINT_SMOOTH Wt::WGLWidget::NONE
#define GL_LINE_SMOOTH Wt::WGLWidget::NONE
#define GL_POLYGON_SMOOTH Wt::WGLWidget::NONE
#define GL_LIGHT0  Wt::WGLWidget::NONE
#define GL_AMBIENT Wt::WGLWidget::NONE
#define GL_DIFFUSE Wt::WGLWidget::NONE
#define GL_POSITION Wt::WGLWidget::NONE

#define GLenum Wt::WGLWidget::GLenum
#define GLchar char
typedef unsigned char GLboolean;
typedef unsigned int GLbitfield;
typedef void GLvoid;
typedef char GLbyte;
typedef short GLshort;
typedef int GLint;
typedef unsigned char GLubyte;
typedef unsigned short GLushort;
typedef unsigned int GLuint;
typedef int GLsizei;
typedef float GLfloat;
typedef float GLclampf;
typedef double GLdouble;
typedef double GLclampd;


#else

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


#endif // G4VIS_BUILD_OPENGLQT_DRIVER


class G4OpenGLViewer;

class G4OpenGLVboDrawer {
public:
  G4OpenGLVboDrawer (G4OpenGLViewer*, std::string type);
  // Create a new OpenGL Drawer. Type could be one of the following :
  // OGL-ES, OGL-Stored, OGL-Immediate, OGL-VBO
  
  virtual ~G4OpenGLVboDrawer ();

// WT specific
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
  void vboGlClear(Wt::WFlags< GLenum > mask);
  void vboGlUniformMatrix4(const Wt::WGLWidget::UniformLocation &location, const Wt::WMatrix4x4 &mat);
  void vboGlUniformMatrix4(const Wt::WGLWidget::UniformLocation &location, const double* matrix);
  void vboGlUniformMatrix4fv(const Wt::WGLWidget::UniformLocation &location, const float* matrix);
  void vboGlUniformMatrix4(const Wt::WGLWidget::UniformLocation &location, const Wt::WGLWidget::JavaScriptMatrix4x4 &mat);
  Wt::WGLWidget::Buffer vboGlCreateBuffer();
  void vboGlBindBuffer(GLenum target, Wt::WGLWidget::Buffer buffer);
  void vboGlDeleteBuffer(Wt::WGLWidget::Buffer buffer);
  void vboGlVertexAttribPointer(Wt::WGLWidget::AttribLocation location, int size, GLenum type, bool normalized, unsigned stride, unsigned offset);
  void vboGlShaderSource(Wt::WGLWidget::Shader shader, GLsizei , const GLchar **src, const GLint *);
  void vboGlCompileShader(Wt::WGLWidget::Shader shader);
  Wt::WGLWidget::Shader vboGlCreateShader(GLenum shader);
  Wt::WGLWidget::Program vboGlCreateProgram();
  void vboGlAttachShader(Wt::WGLWidget::Program program, Wt::WGLWidget::Shader shader);
  void vboGlLinkProgram(Wt::WGLWidget::Program program);
  void vboGlUseProgram(Wt::WGLWidget::Program program);
  void vboGlEnableVertexAttribArray(Wt::WGLWidget::AttribLocation pointer);
  void vboGlDisableVertexAttribArray(Wt::WGLWidget::AttribLocation pointer);
  Wt::WGLWidget::UniformLocation vboGlGetUniformLocation(Wt::WGLWidget::Program programm,const std::string &src);
  Wt::WGLWidget::AttribLocation vboGlGetAttribLocation(Wt::WGLWidget::Program shader,const std::string &src);

  void vboGlClearColor (double r, double g, double b, double a);
  void vboGlClearDepth(double depth);
  void vboGlViewport(int x, int y, unsigned width, unsigned height);
  void vboGlEnable(GLenum cap);
  void vboGlDisable(GLenum cap);
  void vboGlBlendFunc (GLenum sfactor, GLenum dfactor);
  void vboGlDepthFunc (GLenum func);
  void vboGlDepthMask(bool flag);
  void vboGlColorMask (bool red, bool green, bool blue, bool alpha);
  void vboGlLineWidth(double width);
  void vboGlDrawArrays(GLenum mode, int first, unsigned count);
  void vboGlBufferDatafv(GLenum target, const std::vector<double>::iterator begin, const std::vector<double>::iterator end, GLenum usage);
  void vboGlBufferDataiv(GLenum target, const std::vector<unsigned short>::iterator begin, const std::vector<unsigned short>::iterator end, GLenum usage, GLenum type);
  void vboGlDrawElements(GLenum mode, unsigned count, GLenum type, unsigned offset);
  void vboGlMultMatrixf( const GLfloat *m );
  void vboGlMultMatrixd( const GLdouble *m );
#else
  void vboGlMultMatrixf( const GLfloat *m );
  void vboGlMultMatrixd( const GLdouble *m );
#endif
  
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
  
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
  G4OpenGLImmediateWtViewer* fVboViewer;
#else
  G4OpenGLImmediateQtViewer* fVboViewer;
#endif // G4VIS_BUILD_OPENGLWT_DRIVER
};

#endif // G4OPENGL_VERSION_2

#endif // G4OpenGLVboDrawer_HH

