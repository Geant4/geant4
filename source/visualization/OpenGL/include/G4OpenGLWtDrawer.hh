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
// $Id: G4OpenGLWtViewer.hh 74103 2013-09-23 07:52:38Z lgarnier $
//
//
// G4OpenGLWtViewer : Class to provide WindowsNT specific
//                       functionality for OpenGL in GEANT4

#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
#ifndef G4OPENGLWTDRAWER_HH
#define G4OPENGLWTDRAWER_HH

#  include <Wt/WGLWidget>

// specific definition for WT :
// WARNING fWtDrawer should be the exact name of the object!

#define glViewport(a, b, width, height) (fWtDrawer->wglViewport(a, b, width, height))
#define glEnable fWtDrawer->wglEnable
#define glDisable fWtDrawer->wglDisable
#define glBlendFunc fWtDrawer->wglBlendFunc
#define glClear(a) (fWtDrawer->wglClear(a))
#define glClearColor fWtDrawer->wglClearColor
#define glClearDepth fWtDrawer->wglClearDepth
#define glDepthFunc fWtDrawer->wglDepthFunc
#define glDepthMask fWtDrawer->wglDepthMask
#define glFlush() fWtDrawer->wglFlush()
#define glColorMask fWtDrawer->wglColorMask
#define glLineWidth fWtDrawer->wglLineWidth
#define glUniformMatrix4 fWtDrawer->wglUniformMatrix4
#define glDrawArrays fWtDrawer->wglDrawArrays
#define glCreateBuffer fWtDrawer->wglCreateBuffer
#define glVertexPointer fWtDrawer->wglVertexPointer
#define glBindBuffer fWtDrawer->wglBindBuffer
#define glDeleteBuffer fWtDrawer->wglDeleteBuffer
#define glBufferDatafv  fWtDrawer->wglBufferDatafv
#define glBufferDataiv  fWtDrawer->wglBufferDataiv
#define glGetAttribLocation fWtDrawer->wglGetAttribLocation
#define glEnableVertexAttribArray fWtDrawer->wglEnableVertexAttribArray
#define glDisableVertexAttribArray fWtDrawer->wglDisableVertexAttribArray
#define glShaderSource fWtDrawer->wglShaderSource
#define glCompileShader fWtDrawer->wglCompileShader
#define glCreateShader fWtDrawer->wglCreateShader
#define glCreateProgram fWtDrawer->wglCreateProgram
#define glAttachShader fWtDrawer->wglAttachShader
#define glLinkProgram fWtDrawer->wglLinkProgram
#define glUseProgram fWtDrawer->wglUseProgram
#define glDrawElements fWtDrawer->wglDrawElements
#define glVertexAttribPointer fWtDrawer->wglVertexAttribPointer
#define glGetUniformLocation fWtDrawer->wglGetUniformLocation
#define glPointSize fWtDrawer->wglPointSize
#define glColor4d fWtDrawer->wglColor4d
#define glColor4fv fWtDrawer->wglColor4fv
#define glMultMatrixd fWtDrawer->wglMultMatrixd
#define glGetUniformLocation fWtDrawer->wglGetUniformLocation
#define glGetAttribLocation fWtDrawer->wglGetAttribLocation

#define glEdgeFlag(a) (printf("?"))// FIXME : Small debug message for: printf("glEdgeFlag Not implemented\n"))
#define glRenderMode(a) (printf("glRenderMode Not implemented\n"))
#define glClipPlane(a, b) (printf("glClipPane Not implemented\n"))
#define glGetIntegerv(a, b) (printf("glGetIntegerv Not implemented %d %d\n",a,*b))
#define glGetFloatv(a, b) (printf("glGetFloatv Not implemented\n"))
#define glGetDoublev(a,b) (printf("glDoublev Not implemented\n"))
#define glPassThrough(a) (printf("pathTrough Not implemented\n"))
#define glGetBooleanv(a,b) (printf("glGetBooleanv Not implemented\n"))
#define glLoadName(a) (printf("glLoadName Not implemented\n"))
#define glColor3d(a,b,c) (printf("glColor3d Not implemented\n"))
#define glMatrixMode(a) (printf("glMatrixMode Not implemented\n"))
#define glPushMatrix() (printf("glPushMatrix Not implemented\n"))
#define glLoadIdentity() (printf("glLoadIdentity Not implemented\n"))
#define glOrtho(a,b,c,d,e,f) (printf("glOrtho Not implemented %f %f %f %f %f %f\n",a,b,c,d,e,f))
#define glPopMatrix() (printf("glPopMatrix Not implemented\n"))
#define glCallList(a) (printf("glCallList Not implemented\n"))
#define glGenLists(a) (printf("glGenLists Not implemented\n"))
#define glGetError() (printf("glGetError Not implemented\n"))
#define glVertex3d(a,b,c) (printf("glVertex3d Not implemented\n"))
#define glBegin (printf("glBegin Not implemented\n"))
#define glEnd() (printf("glEnd Not implemented\n"))
#define glNewList(a,b) (printf("glNewList Not implemented\n"))
#define glEndList() (printf("glEndList Not implemented\n"))
#define glPolygonMode(a,b) (printf("glPolygonMode Not implemented\n"))
#define glDrawBuffer(a) (printf("glDrawBuffer Not implemented\n"))
#define glDeleteLists(a,b) (printf("glDeleteLists Not implemented\n"))
#define glStencilFunc(a,b,c) (printf("glStencilFunc Not implemented\n"))
#define glStencilOp(a,b,c) (printf("glStencilOp Not implemented\n"))
#define glColorMaterial(a,b) (printf("glColorMaterial Not implemented\n"))
#define glLightfv(a,b,c) (printf("glLightfv Not implemented %u %u %f\n",a,b,*c))
#define glScaled(a,b,c) (printf("glScaled Not implemented\n"))
#define glFrustum(a,b,c,d,e,f) (printf("glFrustum Not implemented\n"))
#define gluLookAt(a,b,c,d,e,f,g,h,i) (printf("gluLookAt Not implemented %f %f %f %f %f %f %f %f %f\n",a,b,c,d,e,f,g,h,i))
#define gluPickMatrix(a,b,c,d,e) (printf("gluPickMatrix Not implemented %f %f %f %f %d\n",a,b,c,d,*e))
#define glSelectBuffer(a,b) (printf("glSelectBuffer Not implemented\n"))
#define glInitNames() (printf("glInitNames Not implemented\n"))
#define glPushNames(a) (printf("glPushNames Not implemented\n"))
#define glPushName(a) (printf("glPushName Not implemented\n"))
#define glPixelStorei(a,b) (printf("glPixelStorei Not implemented\n"))
#define glRasterPos3d(a,b,c) (printf("glRasterPos3d Not implemented\n"))
#define Geant4_gl2psTextOpt(a,b,c,d,e) (printf("gl2psTextOpt Not implemented\n"))
#define glMaterialfv(a,b,c)  (printf("glMaterialfv Not implemented\n"))
#define glCullFace(a) (printf("glCullFace Not implemented\n"))
#define glReadBuffer(a) (printf("glReadBuffer Not implemented\n"))
#define glReadPixels(a,b,c,d,e,f,g) (printf("glReadPixels Not implemented\n"))

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
#define GL_UNSIGNED_BYTE UNSIGNED_BYTE
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
#define GL_PROJECTION 0
#define GL_MODELVIEW 0
#define GL_CLIP_PLANE0 Wt::WGLWidget::ZERO
#define GL_CLIP_PLANE1 Wt::WGLWidget::ZERO
#define GL_CLIP_PLANE2 Wt::WGLWidget::ZERO
#define GL_CLIP_PLANE3 Wt::WGLWidget::ZERO
#define GL_CLIP_PLANE4 Wt::WGLWidget::ZERO
#define GL_COLOR_MATERIAL Wt::WGLWidget::ZERO
#define GL_AMBIENT_AND_DIFFUSE Wt::WGLWidget::ZERO
#define GL_POLYGON 0
#define GL_LIGHTING Wt::WGLWidget::ZERO
#define GL_POINT_SMOOTH Wt::WGLWidget::ZERO
#define GL_LINE_SMOOTH Wt::WGLWidget::ZERO
#define GL_POLYGON_SMOOTH Wt::WGLWidget::ZERO
#define GL_LIGHT0  Wt::WGLWidget::ZERO
#define GL_AMBIENT Wt::WGLWidget::ZERO
#define GL_DIFFUSE Wt::WGLWidget::ZERO
#define GL_POSITION Wt::WGLWidget::ZERO

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

// GL2PS
#define GL2PS_TEXT_B  4
#define GL2PS_TEXT_BL 5
#define GL2PS_TEXT_BR 6

class G4OpenGLViewer;
class G4OpenGLImmediateWtViewer;

class G4OpenGLWtDrawer {
public:
  G4OpenGLWtDrawer (G4OpenGLViewer*);
  virtual ~G4OpenGLWtDrawer ();
  void wglClearColor (double r, double g, double b, double a);
  void wglClearDepth(double depth);
  void wglClear(Wt::WFlags< GLenum > mask);
  void wglFlush();
  void wglViewport(int x, int y, unsigned width, unsigned height);
  void wglEnable(GLenum cap);
  void wglDisable(GLenum cap);
  void wglBlendFunc (GLenum sfactor, GLenum dfactor);
  void wglDepthFunc (GLenum func);
  void wglDepthMask(bool flag);
  void wglColorMask (bool red, bool green, bool blue, bool alpha);
  void wglLineWidth(double width);
  void wglUniformMatrix4(const Wt::WGLWidget::UniformLocation &location, const Wt::WMatrix4x4 &mat);
  void wglUniformMatrix4(const Wt::WGLWidget::UniformLocation &location, const Wt::WGLWidget::JavaScriptMatrix4x4 &mat);
  void wglDrawArrays(GLenum mode, int first, unsigned count);
  Wt::WGLWidget::Buffer wglCreateBuffer();
  void wglBindBuffer(GLenum target, Wt::WGLWidget::Buffer buffer);
  void wglDeleteBuffer(Wt::WGLWidget::Buffer buffer);
  void wglBufferDatafv(GLenum target, const std::vector<double>::iterator begin, const std::vector<double>::iterator end, GLenum usage);
  void wglBufferDataiv(GLenum target, const std::vector<unsigned short>::iterator begin, const std::vector<unsigned short>::iterator end, GLenum usage);
  void wglDrawElements(GLenum mode, unsigned count, GLenum type, unsigned offset);
  void wglVertexAttribPointer(Wt::WGLWidget::AttribLocation location, int size, GLenum type, bool normalized, unsigned stride, unsigned offset);
  void wglPointSize(float size);
  void wglColor4d(int red,int green,int blue,int alpha);
  void wglColor4fv(const GLfloat*);
  void wglMultMatrixd( const GLdouble *m );
  void wglShaderSource(Wt::WGLWidget::Shader shader, GLsizei , const GLchar **src, const GLint *);
  void wglCompileShader(Wt::WGLWidget::Shader shader);
  Wt::WGLWidget::Shader wglCreateShader(GLenum shader);
  Wt::WGLWidget::Program wglCreateProgram();
  void wglAttachShader(Wt::WGLWidget::Program program, Wt::WGLWidget::Shader shader);
  void wglLinkProgram(Wt::WGLWidget::Program program);
  void wglUseProgram(Wt::WGLWidget::Program program);
  void wglEnableVertexAttribArray(Wt::WGLWidget::AttribLocation pointer);
  void wglDisableVertexAttribArray(Wt::WGLWidget::AttribLocation pointer);
  Wt::WGLWidget::UniformLocation wglGetUniformLocation(Wt::WGLWidget::Program programm,const std::string &src);
  Wt::WGLWidget::AttribLocation wglGetAttribLocation(Wt::WGLWidget::Shader shader,const std::string &src);

  private:
  G4OpenGLImmediateWtViewer* fWtViewer;
  
};

#endif

#endif
