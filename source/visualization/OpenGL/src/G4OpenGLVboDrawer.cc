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
// $Id: G4OpenGLVboDrawer.cc 74103 2013-09-23 07:52:38Z lgarnier $
//
//
// G4OpenGLWtViewer : Class to provide Vertex Buffer Object (VBO) specific
//                     functionality for OpenGL > 2.0 in GEANT4
//

#include "G4OpenGLViewer.hh"
#ifdef G4OPENGL_VERSION_2

#include "G4OpenGLVboDrawer.hh"

#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
#include "G4OpenGLImmediateWtViewer.hh"
#else
#include "G4OpenGLImmediateQtViewer.hh"
#endif


//////////////////////////////////////////////////////////////////////////////
G4OpenGLVboDrawer::G4OpenGLVboDrawer (G4OpenGLViewer* viewer,
                                      std::string type
                                      ):
fVboViewer(NULL),
fOGLType(type)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
  G4OpenGLImmediateWtViewer* v = dynamic_cast<G4OpenGLImmediateWtViewer*>(viewer);
#else
  G4OpenGLImmediateQtViewer* v = dynamic_cast<G4OpenGLImmediateQtViewer*>(viewer);
#endif
  if (v) {
    fVboViewer = v;
  }
  
  fFragmentShaderSrc =
  "#ifdef GL_ES\n"
  "precision highp float;\n"
  "#endif\n"
  "\n"
  "varying vec3 vLightWeighting;\n"
  "uniform vec4 uPointColor; // Point Color\n"
  "\n"
  "void main(void) {\n"
  "  vec4 matColor = uPointColor;\n"
  "  gl_FragColor = vec4(matColor.rgb, matColor.a);\n"
  "}\n";
  
  
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
  fVertexShaderSrc =
  "attribute vec3 aVertexPosition;\n"
  "attribute vec3 aVertexNormal;\n"
  "\n"
  "uniform mat4 uMVMatrix; // [M]odel[V]iew matrix\n"
  "uniform mat4 uCMatrix;  // Client-side manipulated [C]amera matrix\n"
  "uniform mat4 uPMatrix;  // Perspective [P]rojection matrix\n"
  "uniform mat4 uNMatrix;  // [N]ormal transformation\n"
  "// uNMatrix is the transpose of the inverse of uCMatrix * uMVMatrix\n"
  "uniform mat4 uTMatrix;  // [T]ransformation  matrix\n"
  "uniform float uPointSize;  // Point size\n"
  "\n"
  "varying vec3 vLightWeighting;\n"
  "\n"
  "void main(void) {\n"
  "  // Calculate the position of this vertex\n"
  "  gl_Position = uPMatrix * uCMatrix * uMVMatrix * uTMatrix * vec4(aVertexPosition, 1.0);\n"
  "\n"
  "  // Phong shading\n"
  "  vec3 transformedNormal = normalize((uNMatrix * vec4(normalize(aVertexNormal), 0)).xyz);\n"
  "  vec3 lightingDirection = normalize(vec3(1, 1, 1));\n"
  "  float directionalLightWeighting = max(dot(transformedNormal, lightingDirection), 0.0);\n"
  "  vec3 uAmbientLightColor = vec3(0.2, 0.2, 0.2);\n"
  "  vec3 uDirectionalColor = vec3(0.8, 0.8, 0.8);\n"
  "  gl_PointSize = uPointSize;\n"
  "  vLightWeighting = uAmbientLightColor + uDirectionalColor * directionalLightWeighting;\n"
  "}\n";
  
#else
  
  
  fVertexShaderSrc =
  "attribute highp vec4 aVertexPosition;\n"
  "attribute vec3 aVertexNormal;\n"
  "uniform highp mat4 uCMatrix;\n"
  "uniform highp mat4 uPMatrix;  // Perspective [P]rojection matrix\n"
  "uniform highp mat4 uMVMatrix; // [M]odel[V]iew matrix\n"
  "uniform highp mat4 uTMatrix;  // [T]ransformation  matrix\n"
  "uniform float uPointSize;  // Point size\n"
  "void main(void)\n"
  "{\n"
  "   gl_Position = uPMatrix * uCMatrix * uMVMatrix * uTMatrix * aVertexPosition;\n"
  "  // Phong shading\n"
  //                                               "  vec3 transformedNormal = normalize((uNMatrix * vec4(normalize(aVertexNormal), 0)).xyz);\n"
  "  vec3 lightingDirection = normalize(vec3(1, 1, 1));\n"
  //                                               "  float directionalLightWeighting = max(dot(transformedNormal, lightingDirection), 0.0);\n"
  //                                               "  vec3 uAmbientLightColor = vec3(0.2, 0.2, 0.2);\n"
  //                                               "  vec3 uDirectionalColor = vec3(0.8, 0.8, 0.8);\n"
  "  gl_PointSize = uPointSize;\n"
  //                                               "  vLightWeighting = uAmbientLightColor + uDirectionalColor * directionalLightWeighting;\n"
  "}";
#endif
  
  
}

//////////////////////////////////////////////////////////////////////////////
G4OpenGLVboDrawer::~G4OpenGLVboDrawer (
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
}
// +--------------------------------+
// +        WT (OpenGL ES) case     +
// +--------------------------------+

#ifdef G4VIS_BUILD_OPENGLWT_DRIVER

void G4OpenGLVboDrawer:: vboGlClear(Wt::WFlags< GLenum > mask) {
  if (fVboViewer->isInitialized()) {
    fVboViewer->clear(mask);
  }
}


void G4OpenGLVboDrawer:: vboGlUniformMatrix4(const Wt::WGLWidget::UniformLocation &location, const Wt::WMatrix4x4 &m) {
  if (fVboViewer) {
    vboGlUseProgram(fVboViewer->getShaderProgram());
    fVboViewer->uniformMatrix4(location, m);
  }
}


void G4OpenGLVboDrawer:: vboGlUniformMatrix4(const Wt::WGLWidget::UniformLocation &location, const double* matrix) {
  if (fVboViewer) {
    Wt::WMatrix4x4 mat(
                       (double) matrix[0], (double) matrix[4], (double) matrix[8], (double) matrix[12],
                       (double) matrix[1], (double) matrix[5], (double) matrix[9], (double) matrix[13],
                       (double) matrix[2], (double) matrix[6], (double) matrix[10],(double) matrix[14],
                       (double) matrix[3],(double) matrix[7],(double) matrix[11],(double) matrix[15]);
    
    fVboViewer->uniformMatrix4(location, mat);
  }
}

void G4OpenGLVboDrawer:: vboGlUniformMatrix4fv(const Wt::WGLWidget::UniformLocation &location, const float* matrix) {
  if (fVboViewer) {
    Wt::WMatrix4x4 mat(
                       (double) matrix[0], (double) matrix[4], (double) matrix[8], (double) matrix[12],
                       (double) matrix[1], (double) matrix[5], (double) matrix[9], (double) matrix[13],
                       (double) matrix[2], (double) matrix[6], (double) matrix[10],(double) matrix[14],
                       (double) matrix[3],(double) matrix[7],(double) matrix[11],(double) matrix[15]);
    
    fVboViewer->uniformMatrix4(location, mat);
  }
}

void G4OpenGLVboDrawer:: vboGlUniformMatrix4(const Wt::WGLWidget::UniformLocation &location, const Wt::WGLWidget::JavaScriptMatrix4x4 &m) {
  if (fVboViewer->isInitialized()) {
    fVboViewer->uniformMatrix4(location, m);
  }
}

Wt::WGLWidget::Buffer G4OpenGLVboDrawer:: vboGlCreateBuffer(){
  if (fVboViewer) {
    return fVboViewer->createBuffer();
  }
  return Wt::WGLWidget::Buffer();
}

void G4OpenGLVboDrawer:: vboGlBindBuffer(GLenum target, Wt::WGLWidget::Buffer buffer){
  if (fVboViewer) {
    fVboViewer->bindBuffer(target,buffer);
  }
}

void G4OpenGLVboDrawer:: vboGlDeleteBuffer(Wt::WGLWidget::Buffer buffer){
  if (fVboViewer == NULL) return;
  fVboViewer->deleteBuffer(buffer);
}

void G4OpenGLVboDrawer:: vboGlVertexAttribPointer(Wt::WGLWidget::AttribLocation location, int size, GLenum type, bool normalized, unsigned stride, unsigned offset){
  if (fVboViewer->isInitialized()) {
    fVboViewer->vertexAttribPointer(location, size,type, normalized, stride, offset);
  }
}

Wt::WGLWidget::Program G4OpenGLVboDrawer:: vboGlCreateProgram(){
  if (fVboViewer) {
    return fVboViewer->createProgram();
  }
  return Wt::WGLWidget::Program();
}

void G4OpenGLVboDrawer:: vboGlAttachShader(Wt::WGLWidget::Program program, Shader shader){
  if (fVboViewer) {
    fVboViewer->attachShader(program,shader);
  }
}

void G4OpenGLVboDrawer:: vboGlLinkProgram(Wt::WGLWidget::Program program){
  if (fVboViewer) {
    fVboViewer->linkProgram(program);
  }
}

void G4OpenGLVboDrawer:: vboGlUseProgram(Wt::WGLWidget::Program program){
  if (fVboViewer) {
    fVboViewer->useProgram(program);
  }
}

void G4OpenGLVboDrawer:: vboGlEnableVertexAttribArray(Wt::WGLWidget::AttribLocation pointer){
  if (fVboViewer) {
    fVboViewer->enableVertexAttribArray(pointer);
  }
}

void G4OpenGLVboDrawer:: vboGlDisableVertexAttribArray(Wt::WGLWidget::AttribLocation pointer){
  if (fVboViewer) {
    fVboViewer->disableVertexAttribArray(pointer);
  }
}

Wt::WGLWidget::UniformLocation G4OpenGLVboDrawer:: vboGlGetUniformLocation(Wt::WGLWidget::Program programm,const std::string &src){
  if (fVboViewer) {
    return fVboViewer->getUniformLocation(programm, src);
  }
  return Wt::WGLWidget::UniformLocation();
}

Wt::WGLWidget::AttribLocation G4OpenGLVboDrawer:: vboGlGetAttribLocation(Wt::WGLWidget::Program shader,const std::string &src){
  if (fVboViewer) {
    return fVboViewer->getAttribLocation(shader, src);
  }
  return Wt::WGLWidget::AttribLocation();
}





void G4OpenGLVboDrawer::vboGlClearColor (double r, double g, double b, double a) {
  
  if (fVboViewer->isInitialized() ) {
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
    fVboViewer->clearColor(r,g,b,a);
#else
#endif
  }
}

void G4OpenGLVboDrawer::vboGlClearDepth(double depth) {
  if (fVboViewer->isInitialized()) {
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
    fVboViewer->clearDepth(depth);
#else
    glClearDepth(depth);
#endif
  }
}


void G4OpenGLVboDrawer::vboGlFlush() {
  if (fVboViewer->isInitialized()) {
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
    fVboViewer->flush();
#else
#endif
  }
}

void G4OpenGLVboDrawer:: vboGlViewport(int x, int y, unsigned width, unsigned height) {
  if (fVboViewer->isInitialized()) {
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
    fVboViewer->viewport(x,y,width,height);
#else
#endif
  }
}

void G4OpenGLVboDrawer:: vboGlEnable(GLenum cap) {
  if (fVboViewer->isInitialized()) {
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
    if (cap != Wt::WGLWidget::NONE) {
      fVboViewer->enable(cap);
    }
#else
#endif
  }
}

void G4OpenGLVboDrawer:: vboGlDisable(GLenum cap) {
  if (fVboViewer->isInitialized()) {
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
    if (cap != Wt::WGLWidget::NONE) {
      fVboViewer->disable(cap);
    }
#else
#endif
  }
}

void G4OpenGLVboDrawer:: vboGlBlendFunc (GLenum sfactor, GLenum dfactor) {
  if (fVboViewer->isInitialized()) {
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
    fVboViewer->blendFunc(sfactor, dfactor);
#else
#endif
    
  }
}

void G4OpenGLVboDrawer:: vboGlDepthFunc (GLenum func) {
  if (fVboViewer->isInitialized()) {
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
    fVboViewer->depthFunc(func);
#else
#endif
  }
}

void G4OpenGLVboDrawer:: vboGlDepthMask(bool flag) {
  if (fVboViewer->isInitialized()) {
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
    fVboViewer->depthMask(flag);
#else
#endif
  }
}

void G4OpenGLVboDrawer:: vboGlColorMask (bool red, bool green, bool blue, bool alpha) {
  if (fVboViewer->isInitialized()) {
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
    fVboViewer->colorMask(red,green,blue,alpha);
#else
#endif
  }
}

void G4OpenGLVboDrawer:: vboGlLineWidth(double width) {
  if (fVboViewer->isInitialized()) {
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
    fVboViewer->lineWidth(width);
#else
#endif
  }
}


void G4OpenGLVboDrawer:: vboGlDrawArrays(GLenum mode, int first, unsigned count){
  if (fVboViewer->isInitialized()) {
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
    fVboViewer->drawArrays(mode,first, count);
#else
#endif
  }
}

void G4OpenGLVboDrawer:: vboGlDrawElements(GLenum mode, unsigned count, GLenum type, unsigned offset){
  if (fVboViewer->isInitialized()) {
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
    fVboViewer->drawElements(mode,count,type,offset);
#else
#endif
  }
}

void G4OpenGLVboDrawer:: vboGlBufferDatafv(GLenum target, const std::vector<double>::iterator begin, const std::vector<double>::iterator end, GLenum usage){
  if (fVboViewer) {
    if (fVboViewer->isInitialized()) {
      fVboViewer->bufferDatafv(target,begin,end,usage);
    }
  }
}

void G4OpenGLVboDrawer:: vboGlBufferDataiv(GLenum target, const std::vector<unsigned short>::iterator begin, const std::vector<unsigned short>::iterator end, GLenum usage, GLenum type){
  if (fVboViewer) {
    if (fVboViewer->isInitialized()) {
      fVboViewer->bufferDataiv(target,begin,end,usage,type);
    }
  }
}


void G4OpenGLVboDrawer:: vboGlMultMatrixd(const GLdouble *matrix){
  if (fVboViewer) {
    if (fVboViewer->isInitialized()) {
      Wt::WMatrix4x4 mat(
                         (double) matrix[0], (double) matrix[4], (double) matrix[8], (double) matrix[12],
                         (double) matrix[1], (double) matrix[5], (double) matrix[9], (double) matrix[13],
                         (double) matrix[2], (double) matrix[6], (double) matrix[10],(double) matrix[14],
                         (double) matrix[3],(double) matrix[7],(double) matrix[11],(double) matrix[15]);
      
      // FIXME !
      fVboViewer->uniformMatrix4(fVboViewer->getShaderTransformMatrix(),mat);
      if (fMatrixMode == GL_MODELVIEW) {
        //      fVboViewer->uniformMatrix4(fVboViewer->getShaderTransformMatrix(),mat);
      } else {
        G4cerr << "glMultMatrixd could only be used in GL_MODELVIEW mode" << G4endl;
      }
    }
  }
}

void G4OpenGLVboDrawer:: vboGlMultMatrixf(const GLfloat *matrix){
  if (fVboViewer) {
    if (fVboViewer->isInitialized()) {
      Wt::WMatrix4x4 mat(
                         matrix[0], matrix[4], matrix[8], matrix[12],
                         matrix[1], matrix[5], matrix[9], matrix[13],
                         matrix[2], matrix[6], matrix[10], matrix[14],
                         matrix[3], matrix[7], matrix[11], matrix[15]);
      
      if (fMatrixMode == GL_MODELVIEW) {
        fVboViewer->uniformMatrix4(fVboViewer->getShaderTransformMatrix(),mat);
      } else {
        G4cerr << "glMultMatrixf could only be used in GL_MODELVIEW mode" << G4endl;
      }
    }
  }
}

void G4OpenGLVboDrawer:: vboGlShaderSource(Shader shader, GLsizei , const GLchar **src, const GLint *){
  if (fVboViewer) {
    std::string s = *src;
    fVboViewer->shaderSource(shader, s);
  }
}

void G4OpenGLVboDrawer:: vboGlCompileShader(Shader shader){
  if (fVboViewer) {
    fVboViewer->compileShader(shader);
  }
}

Shader G4OpenGLVboDrawer:: vboGlCreateShader(GLenum shader){
  if (fVboViewer) {
    return fVboViewer->createShader(shader);
  }
  return Shader();
}

#else

// +--------------------------------+
// +        QT (OpenGL ES) case     +
// +--------------------------------+

void G4OpenGLVboDrawer:: vboGlMultMatrixf(const GLfloat *matrix){
  if (fVboViewer) {
    if (fVboViewer->isInitialized()) {
      // FIXME
      // glUniformMatrix4fv(12, 1, 0, 0x7fff5fbf5d00)
      //  Error: GL_INVALID_OPERATION
      
      if (fMatrixMode == GL_MODELVIEW) {
        glUniformMatrix4fv(fVboViewer->getShaderTransformMatrix(),1,0,matrix);
      } else {
        G4cerr << "glMultMatrixf could only be used in GL_MODELVIEW mode" << G4endl;
      }
    }
  }
}


void G4OpenGLVboDrawer:: vboGlMultMatrixd(const GLdouble *matrix){
  if (fVboViewer) {
    if (fVboViewer->isInitialized()) {
      // FIXME !
      //    if (fMatrixMode == GL_MODELVIEW) {
      //      printf("G4OpenGLVboDrawer:: vboGlMultMatrixd %d %d\n",fVboViewer->getShaderTransformMatrix(), matrix);
      //!! TEST !!
      float mat[16] = {
        matrix[0],matrix[1],matrix[2],matrix[3],
        matrix[4],matrix[5],matrix[6],matrix[7],
        matrix[8],matrix[9],matrix[10],matrix[11],
        matrix[12],matrix[13],matrix[14],matrix[15]
      };

      glUniformMatrix4fv(fVboViewer->getShaderTransformMatrix(),1,0,mat);
      GLenum e = glGetError();
      printf("GL error : %d",e);
      //    } else {
      //      G4cerr << "glMultMatrixd could only be used in GL_MODELVIEW mode" << G4endl;
      //    }
    }
  }
}



#endif
// +--------------------------------+
// +        All case     +
// +--------------------------------+



void G4OpenGLVboDrawer::vboGlOrtho(GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble zNear, GLdouble zFar) {
  if (fVboViewer) {
    if (fVboViewer->isInitialized()) {
      printf("glOrtho implemented --- %f %f %f %f %f %f \n",left, right, bottom, top, zNear, zFar);
      float a = 2.0f / (right - left);
      float b = 2.0f / (top - bottom);
      float c = -2.0f / (zFar - zNear);
      
      float tx = - (right + left)/(right - left);
      float ty = - (top + bottom)/(top - bottom);
      float tz = - (zFar + zNear)/(zFar - zNear);
      
      float ortho[16] = {
        a, 0, 0, 0,
        0, b, 0, 0,
        0, 0, c, 0,
        tx, ty, tz, 1
      };
      // FIXME :
      //  glUniformMatrix4fv(0, 1, 0, 0x7fff5fbf5d00)
      // Error: GL_INVALID_OPERATION
      
      if (fMatrixMode == GL_PROJECTION) {
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
        vboGlUniformMatrix4fv(fVboViewer->getShaderProjectionMatrix(), ortho);
#else
        glUniformMatrix4fv(fVboViewer->getShaderProjectionMatrix(), 1, 0, ortho);
#endif
      } else {
        G4cerr << "glFrustum could only be used in GL_PROJECTION mode" << G4endl;
      }
    }
  }
}


void G4OpenGLVboDrawer::vboGlFrustum(GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble zNear, GLdouble zFar) {
  if (fVboViewer) {
    if (fVboViewer->isInitialized()) {
      float deltaX = right - left;
      float deltaY = top - bottom;
      float deltaZ = zFar - zNear;
      
      float a = 2.0f * zNear / deltaX;
      float b = 2.0f * zNear / deltaY;
      float c = (right + left) / deltaX;
      float d = (top + bottom) / deltaY;
      float e = -(zFar + zNear) / (zFar - zNear);
      float f = -2.0f * zFar * zNear / deltaZ;
      
      float proj[16] = {
        a, 0, 0, 0,
        0, b, 0, 0,
        c, d, e, -1.0f,
        0, 0, f, 0
      };
      
      if (fMatrixMode == GL_PROJECTION) {
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
        vboGlUniformMatrix4fv(fVboViewer->getShaderProjectionMatrix(), proj);
#else
        glUniformMatrix4fv(fVboViewer->getShaderProjectionMatrix(), 1, 0, proj);
#endif
      } else {
        G4cerr << "glFrustrum could only be used in GL_PROJECTION mode" << G4endl;
      }
    }
  }
}


void G4OpenGLVboDrawer::vboGlMatrixMode(GLenum a) {
  if (fVboViewer) {
    if (fVboViewer->isInitialized()) {
      printf("G4OpenGLVboDrawer::vboGlMatrixMode CHANGED :%d \n",a);
      fMatrixMode = a;
    }
  }
}


void G4OpenGLVboDrawer::vboGlColor4d(int red,int green,int blue,int alpha) {
  if (fVboViewer) {
    if (fVboViewer->isInitialized()) {
      //    double color [] = { red, green, blue, alpha };
      // FIXME : REMOVE /2 , used to render transparents for testing purpose
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
      double color [] = { red, green, blue, 0.7 };
      glUniform4fv (fVboViewer->getUniformLocation(fVboViewer->getShaderProgram(), "uPointColor"),color);
#else
      alpha = 0.7;
      glUniform4f (glGetUniformLocation(fVboViewer->getShaderProgram(), "uPointColor"),red, green, blue, alpha);
#endif
    }
  }
}

void G4OpenGLVboDrawer:: vboGlColor4fv(const GLfloat* data) {
  if (fVboViewer) {
    if (fVboViewer->isInitialized()) {
      double color [] = { (data[0]), (data[1]), (data[2]), 0.7};
      // FIXME : REMOVE /2 , used to render transparents for testing purpose
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
      glUniform4fv (fVboViewer->getUniformLocation(fVboViewer->getShaderProgram(), "uPointColor"),color);
#else
      glUniform4f (glGetUniformLocation(fVboViewer->getShaderProgram(), "uPointColor"),color[0],color[1],color[2], color[3]);
#endif
    }
  }
}

void G4OpenGLVboDrawer:: vboGlPointSize(float size) {
  if (fVboViewer) {
    if (fVboViewer->isInitialized()) {
#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
      glUniform1f(fVboViewer->getUniformLocation(fVboViewer->getShaderProgram(), "uPointSize"),size);
#else
      glUniform1f (glGetUniformLocation(fVboViewer->getShaderProgram(), "uPointSize"),size);
#endif
    }
  }
}

#endif

