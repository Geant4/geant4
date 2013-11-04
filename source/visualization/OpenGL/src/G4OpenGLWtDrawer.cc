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
// $Id: G4OpenGLWtDrawer.cc 74103 2013-09-23 07:52:38Z lgarnier $
//
//
// G4OpenGLWtViewer : Class to provide Wt specific
//                     functionality for OpenGL in GEANT4
//
// 27/06/2003 : G.Barrand : implementation (at last !).

#ifdef G4VIS_BUILD_OPENGLWT_DRIVER

#include "G4OpenGLWtDrawer.hh"
#include "G4OpenGLImmediateWtViewer.hh"

//////////////////////////////////////////////////////////////////////////////
G4OpenGLWtDrawer::G4OpenGLWtDrawer (G4OpenGLViewer* viewer
                                    ):
fWtViewer(NULL)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  G4OpenGLImmediateWtViewer* v = dynamic_cast<G4OpenGLImmediateWtViewer*>(viewer);
  if (v) {
    fWtViewer = v;
  }
}

//////////////////////////////////////////////////////////////////////////////
G4OpenGLWtDrawer::~G4OpenGLWtDrawer (
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
}

void G4OpenGLWtDrawer::wglClearColor (double r, double g, double b, double a) {
  
  if (fWtViewer) {
    fWtViewer->clearColor(r,g,b,a);
  }
}

void G4OpenGLWtDrawer::wglClearDepth(double depth) {
  if (fWtViewer) {
    fWtViewer->clearDepth(depth);
  }
}

void G4OpenGLWtDrawer::wglClear(Wt::WFlags< GLenum > mask) {
  if (fWtViewer) {
    fWtViewer->clear(mask);
  }
}

void G4OpenGLWtDrawer::wglFlush() {
  if (fWtViewer) {
    fWtViewer->flush();
  }
}

void G4OpenGLWtDrawer::wglViewport(int x, int y, unsigned width, unsigned height) {
  if (fWtViewer) {
    fWtViewer->viewport(x,y,width,height);
  }
}

void G4OpenGLWtDrawer::wglEnable(GLenum cap) {
  if (fWtViewer) {
    fWtViewer->enable(cap);
  }
}

void G4OpenGLWtDrawer::wglDisable(GLenum cap) {
  if (fWtViewer) {
    fWtViewer->disable(cap);
  }
}

void G4OpenGLWtDrawer::wglBlendFunc (GLenum sfactor, GLenum dfactor) {
  if (fWtViewer) {
    fWtViewer->blendFunc(sfactor, dfactor);
  }
}

void G4OpenGLWtDrawer::wglDepthFunc (GLenum func) {
  if (fWtViewer) {
    fWtViewer->depthFunc(func);
  }
}

void G4OpenGLWtDrawer::wglDepthMask(bool flag) {
  if (fWtViewer) {
    fWtViewer->depthMask(flag);
  }
}

void G4OpenGLWtDrawer::wglColorMask (bool red, bool green, bool blue, bool alpha) {
  if (fWtViewer) {
    fWtViewer->colorMask(red,green,blue,alpha);
  }
}

void G4OpenGLWtDrawer::wglLineWidth(double width) {
  if (fWtViewer) {
    fWtViewer->lineWidth(width);
  }
}

void G4OpenGLWtDrawer::wglUniformMatrix4(const Wt::WGLWidget::UniformLocation &location, const Wt::WMatrix4x4 &m) {
  if (fWtViewer) {
    fWtViewer->uniformMatrix4(location, m);
  }
}

void G4OpenGLWtDrawer::wglUniformMatrix4(const Wt::WGLWidget::UniformLocation &location, const Wt::WGLWidget::JavaScriptMatrix4x4 &m) {
  if (fWtViewer) {
    fWtViewer->uniformMatrix4(location, m);
  }
}

void G4OpenGLWtDrawer::wglDrawArrays(GLenum mode, int first, unsigned count){
  if (fWtViewer) {
    fWtViewer->drawArrays(mode,first, count);
  }
}

void G4OpenGLWtDrawer::wglDrawElements(GLenum mode, unsigned count, GLenum type, unsigned offset){
  if (fWtViewer) {
    fWtViewer->drawElements(mode,count,type,offset);
  }
}

Wt::WGLWidget::Buffer G4OpenGLWtDrawer::wglCreateBuffer(){
  if (fWtViewer) {
    return fWtViewer->createBuffer();
  }
  return NULL;
}

void G4OpenGLWtDrawer::wglBindBuffer(GLenum target, Wt::WGLWidget::Buffer buffer){
  if (fWtViewer) {
    fWtViewer->bindBuffer(target,buffer);
  }
}

void G4OpenGLWtDrawer::wglDeleteBuffer(Wt::WGLWidget::Buffer buffer){
  if (fWtViewer == NULL) return;
  fWtViewer->deleteBuffer(buffer);
}

void G4OpenGLWtDrawer::wglBufferDatafv(GLenum target, const std::vector<double>::iterator begin, const std::vector<double>::iterator end, GLenum usage){
  if (fWtViewer) {
    fWtViewer->bufferDatafv(target,begin,end,usage);
  }
}

void G4OpenGLWtDrawer::wglBufferDataiv(GLenum target, const std::vector<unsigned short>::iterator begin, const std::vector<unsigned short>::iterator end, GLenum usage){
  if (fWtViewer) {
    fWtViewer->bufferDataiv(target,begin,end,usage, GL_UNSIGNED_SHORT);
  }
}

void G4OpenGLWtDrawer::wglVertexAttribPointer(Wt::WGLWidget::AttribLocation location, int size, GLenum type, bool normalized, unsigned stride, unsigned offset){
  if (fWtViewer) {
    fWtViewer->vertexAttribPointer(location, size,type, normalized, stride, offset);
  }
}

void G4OpenGLWtDrawer::wglMultMatrixd(const GLdouble *matrix){
  if (fWtViewer) {
    Wt::WMatrix4x4 mat(
                       (double) matrix[0], (double) matrix[4], (double) matrix[8], (double) matrix[12],
                       (double) matrix[1], (double) matrix[5], (double) matrix[9], (double) matrix[13],
                       (double) matrix[2], (double) matrix[6], (double) matrix[10],(double) matrix[14],
                       (double) matrix[3],(double) matrix[7],(double) matrix[11],(double) matrix[15]);
    
    fWtViewer->uniformMatrix4(fWtViewer->getUniformLocation(fWtViewer->shaderProgram_, "uTMatrix"),mat);
  }
}

void G4OpenGLWtDrawer::wglPointSize(float size) {
  if (fWtViewer) {
    fWtViewer->uniform1f (fWtViewer->getUniformLocation(fWtViewer->shaderProgram_, "uPointSize"),size);
  }
}

void G4OpenGLWtDrawer::wglColor4d(int red,int green,int blue,int alpha) {
  if (fWtViewer) {
    //    double color [] = { red, green, blue, alpha };
    double color [] = { red, green, blue, 0.7 };
    // FIXME : REMOVE /2 , used to render transparents for testing purpose
    fWtViewer->uniform4fv (fWtViewer->getUniformLocation(fWtViewer->shaderProgram_, "uPointColor"),color);
  }
}

void G4OpenGLWtDrawer::wglColor4fv(const GLfloat* data) {
  if (fWtViewer) {
    //    double color [] = { (data[0]), (data[1]), (data[2]), (data[3])};
    double color [] = { (data[0]), (data[1]), (data[2]), 0.7};
    // FIXME : REMOVE /2 , used to render transparents for testing purpose
    fWtViewer->uniform4fv (fWtViewer->getUniformLocation(fWtViewer->shaderProgram_, "uPointColor"),color);
  }
}

void G4OpenGLWtDrawer::wglShaderSource(Shader shader, GLsizei , const GLchar **src, const GLint *){
  if (fWtViewer) {
    std::string s = *src;
    fWtViewer->shaderSource(shader, s);
  }
}

void G4OpenGLWtDrawer::wglCompileShader(Shader shader){
  if (fWtViewer) {
    fWtViewer->compileShader(shader);
  }
}

Shader G4OpenGLWtDrawer::wglCreateShader(GLenum shader){
  if (fWtViewer) {
    return fWtViewer->createShader(shader);
  }
  return NULL;
}

Wt::WGLWidget::Program G4OpenGLWtDrawer::wglCreateProgram(){
  if (fWtViewer) {
    return fWtViewer->createProgram();
  }
  return NULL;
}

void G4OpenGLWtDrawer::wglAttachShader(Wt::WGLWidget::Program program, Shader shader){
  if (fWtViewer) {
    fWtViewer->attachShader(program,shader);
  }
}

void G4OpenGLWtDrawer::wglLinkProgram(Wt::WGLWidget::Program program){
  if (fWtViewer) {
    fWtViewer->linkProgram(program);
  }
}

void G4OpenGLWtDrawer::wglUseProgram(Wt::WGLWidget::Program program){
  if (fWtViewer) {
    fWtViewer->useProgram(program);
  }
}

void G4OpenGLWtDrawer::wglEnableVertexAttribArray(Wt::WGLWidget::AttribLocation pointer){
  if (fWtViewer) {
    fWtViewer->enableVertexAttribArray(pointer);
  }
}

void G4OpenGLWtDrawer::wglDisableVertexAttribArray(Wt::WGLWidget::AttribLocation pointer){
  if (fWtViewer) {
    fWtViewer->disableVertexAttribArray(pointer);
  }
}

Wt::WGLWidget::UniformLocation G4OpenGLWtDrawer::wglGetUniformLocation(Wt::WGLWidget::Program programm,const std::string &src){
  if (fWtViewer) {
    return fWtViewer->getUniformLocation(programm, src);
  }
  return NULL;
}

Wt::WGLWidget::AttribLocation G4OpenGLWtDrawer::wglGetAttribLocation(Shader shader,const std::string &src){
  if (fWtViewer) {
    return fWtViewer->getAttribLocation(shader, src);
  }
  return NULL;
}


#endif
