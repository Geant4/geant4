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
// $Id: G4OpenGLViewer.cc 103926 2017-05-03 13:43:27Z gcosmo $
//
// 
// Andrew Walkden  27th March 1996
// OpenGL view - opens window, hard copy, etc.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLSceneHandler.hh"
#include "G4OpenGLTransform3D.hh"
#include "G4OpenGL2PSAction.hh"

#include "G4Scene.hh"
#include "G4VisExtent.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"
#include "G4Plane3D.hh"
#include "G4AttHolder.hh"
#include "G4AttCheck.hh"
#include "G4Text.hh"

#ifdef G4OPENGL_VERSION_2
// We need to have a Wt gl drawer because we will draw inside the WtGL component (ImmediateWtViewer)
#include "G4OpenGLVboDrawer.hh"
#endif

// GL2PS
#include "Geant4_gl2ps.h"

#include <sstream>
#include <string>
#include <iomanip>

G4OpenGLViewer::G4OpenGLViewer (G4OpenGLSceneHandler& scene):
G4VViewer (scene, -1),
#ifdef G4OPENGL_VERSION_2
fVboDrawer(NULL),
#endif
fPrintColour (true),
fVectoredPs (true),
fOpenGLSceneHandler(scene),
background (G4Colour(0.,0.,0.)),
transparency_enabled (true),
antialiasing_enabled (false),
haloing_enabled (false),
fStartTime(-DBL_MAX),
fEndTime(DBL_MAX),
fFadeFactor(0.),
fDisplayHeadTime(false),
fDisplayHeadTimeX(-0.9),
fDisplayHeadTimeY(-0.9),
fDisplayHeadTimeSize(24.),
fDisplayHeadTimeRed(0.),
fDisplayHeadTimeGreen(1.),
fDisplayHeadTimeBlue(1.),
fDisplayLightFront(false),
fDisplayLightFrontX(0.),
fDisplayLightFrontY(0.),
fDisplayLightFrontZ(0.),
fDisplayLightFrontT(0.),
fDisplayLightFrontRed(0.),
fDisplayLightFrontGreen(1.),
fDisplayLightFrontBlue(0.),
fRot_sens(1.),
fPan_sens(0.01),
fWinSize_x(0),
fWinSize_y(0),
fDefaultExportImageFormat("pdf"),
fExportImageFormat("pdf"),
fExportFilenameIndex(0),
fPrintSizeX(-1),
fPrintSizeY(-1),
fPointSize (0),
fDefaultExportFilename("G4OpenGL"),
fSizeHasChanged(0),
fGl2psDefaultLineWith(1),
fGl2psDefaultPointSize(2),
fGlViewInitialized(false),
fIsGettingPickInfos(false)
#ifdef G4OPENGL_VERSION_2
,fShaderProgram(0)
,fVertexPositionAttribute(0)
,fVertexNormalAttribute(0)
,fpMatrixUniform(0)
,fcMatrixUniform(0)
,fmvMatrixUniform(0)
,fnMatrixUniform(0)
#endif
{
  // Make changes to view parameters for OpenGL...
  fVP.SetAutoRefresh(true);
  fDefaultVP.SetAutoRefresh(true);

  fGL2PSAction = new G4OpenGL2PSAction();

  // add supported export image format
  addExportImageFormat("eps");
  addExportImageFormat("ps");
  addExportImageFormat("pdf");
  addExportImageFormat("svg");

  // Change the default name
  fExportFilename += fDefaultExportFilename + "_" + GetShortName().data();
  
  //  glClearColor (0.0, 0.0, 0.0, 0.0);
  //  glClearDepth (1.0);
  //  glDisable (GL_BLEND);
  //  glDisable (GL_LINE_SMOOTH);
  //  glDisable (GL_POLYGON_SMOOTH);

}

G4OpenGLViewer::~G4OpenGLViewer ()
{
  delete fGL2PSAction;
}

void G4OpenGLViewer::InitializeGLView () 
{
#ifdef G4OPENGL_VERSION_2
  if (fVboDrawer) {
    
    // First, load a simple shader
    fShaderProgram = glCreateProgram();
    Shader vertexShader = glCreateShader(GL_VERTEX_SHADER);
    const char * vSrc = fVboDrawer->getVertexShaderSrc();
    glShaderSource(vertexShader, 1, &vSrc, NULL);
    glCompileShader(vertexShader);
    glAttachShader(fShaderProgram, vertexShader);
    
    Shader fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    const char * fSrc = fVboDrawer->getFragmentShaderSrc();
    glShaderSource(fragmentShader, 1, &fSrc, NULL);
    glCompileShader(fragmentShader);
    
    glAttachShader(fShaderProgram, fragmentShader);
    glLinkProgram(fShaderProgram);
    glUseProgram(fShaderProgram);
    
    //   UniformLocation uColor = getUniformLocation(fShaderProgram, "uColor");
    //   uniform4fv(uColor, [0.0, 0.3, 0.0, 1.0]);
    
    // Extract the references to the attributes from the shader.
    
    fVertexPositionAttribute =
    glGetAttribLocation(fShaderProgram, "aVertexPosition");
    
    
    glEnableVertexAttribArray(fVertexPositionAttribute);
    
    // Extract the references the uniforms from the shader
    fpMatrixUniform  = glGetUniformLocation(fShaderProgram, "uPMatrix");
    fcMatrixUniform  = glGetUniformLocation(fShaderProgram, "uCMatrix");
    fmvMatrixUniform = glGetUniformLocation(fShaderProgram, "uMVMatrix");
    fnMatrixUniform  = glGetUniformLocation(fShaderProgram, "uNMatrix");
    ftMatrixUniform  = glGetUniformLocation(fShaderProgram, "uTMatrix");
    
    /*    glUniformMatrix4fv(fcMatrixUniform, 1, 0, identity);
     glUniformMatrix4fv(fpMatrixUniform, 1, 0, identity);
     glUniformMatrix4fv(ftMatrixUniform, 1, 0, identity);
     glUniformMatrix4fv(fmvMatrixUniform, 1, 0, identity);
     */
    // We have to set that in order to avoid calls on opengl commands before all is ready
    fGlViewInitialized = true;
  }
#endif
  
  if (fWinSize_x == 0) {
    fWinSize_x = fVP.GetWindowSizeHintX();
  }
  if (fWinSize_y == 0) {
    fWinSize_y = fVP.GetWindowSizeHintY();
  }

  glClearColor (0.0, 0.0, 0.0, 0.0);
  glClearDepth (1.0);
#ifndef G4OPENGL_VERSION_2
  glDisable (GL_LINE_SMOOTH);
  glDisable (GL_POLYGON_SMOOTH);
#endif

// clear the buffers and window?
  ClearView ();
  FinishView ();
  
  glDepthFunc (GL_LEQUAL);
  glDepthMask (GL_TRUE);
  
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
}

void G4OpenGLViewer::ClearView () {
  ClearViewWithoutFlush();

  if(!isFramebufferReady()) {
    return;
  }

  glFlush();
}


void G4OpenGLViewer::ClearViewWithoutFlush () {
  // Ready for clear ?
  // See : http://lists.apple.com/archives/mac-opengl/2012/Jul/msg00038.html
  if(!isFramebufferReady()) {
    return;
  }
  
  glClearColor (background.GetRed(),
                background.GetGreen(),
                background.GetBlue(),
                1.);
  glClearDepth (1.0);
  //Below line does not compile with Mesa includes.
  //glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glClear (GL_COLOR_BUFFER_BIT);
  glClear (GL_DEPTH_BUFFER_BIT);
  glClear (GL_STENCIL_BUFFER_BIT);
}


void G4OpenGLViewer::ResizeWindow(unsigned int aWidth, unsigned int aHeight) {
  if ((fWinSize_x != aWidth) || (fWinSize_y != aHeight)) {
    fWinSize_x = aWidth;
    fWinSize_y = aHeight;
    fSizeHasChanged = true;
  } else {
    fSizeHasChanged = false;
  }
}

/**
 * Set the viewport of the scene
 * MAXIMUM SIZE is :
 * GLint dims[2];
 * glGetIntegerv(GL_MAX_VIEWPORT_DIMS, dims);
 */
void G4OpenGLViewer::ResizeGLView()
{
  // Check size
  GLint dims[2];
  dims[0] = 0;
  dims[1] = 0;

  glGetIntegerv(GL_MAX_VIEWPORT_DIMS, dims);

  if ((dims[0] !=0 ) && (dims[1] !=0)) {

    if (fWinSize_x > (unsigned)dims[0]) {
      G4cerr << "Try to resize view greater than max X viewport dimension. Desired size "<<fWinSize_x <<" is resize to "<<  dims[0] << G4endl;
      fWinSize_x = dims[0];
    }
    if (fWinSize_y > (unsigned)dims[1]) {
      G4cerr << "Try to resize view greater than max Y viewport dimension. Desired size "<<fWinSize_y <<" is resize to "<<  dims[1] << G4endl;
      fWinSize_y = dims[1];
    }
  }
    
  glViewport(0, 0, fWinSize_x,fWinSize_y);   


}


void G4OpenGLViewer::SetView () {
  // if getting pick infos, should not resize the view.
  if (fIsGettingPickInfos) return;
  
  if (!fSceneHandler.GetScene()) {
    return;
  }
  // Calculates view representation based on extent of object being
  // viewed and (initial) viewpoint.  (Note: it can change later due
  // to user interaction via visualization system's GUI.)
  
  // Lighting.
  GLfloat lightPosition [4];
  lightPosition [0] = fVP.GetActualLightpointDirection().x();
  lightPosition [1] = fVP.GetActualLightpointDirection().y();
  lightPosition [2] = fVP.GetActualLightpointDirection().z();
  lightPosition [3] = 0.;
  // Light position is "true" light direction, so must come after gluLookAt.
  GLfloat ambient [] = { 0.2f, 0.2f, 0.2f, 1.f};
  GLfloat diffuse [] = { 0.8f, 0.8f, 0.8f, 1.f};
  glEnable (GL_LIGHT0);
  glLightfv (GL_LIGHT0, GL_AMBIENT, ambient);
  glLightfv (GL_LIGHT0, GL_DIFFUSE, diffuse);
  
  G4double ratioX = 1;
  G4double ratioY = 1;
  if (fWinSize_y > fWinSize_x) {
    ratioX = ((G4double)fWinSize_y) / ((G4double)fWinSize_x);
  }
  if (fWinSize_x > fWinSize_y) {
    ratioY = ((G4double)fWinSize_x) / ((G4double)fWinSize_y);
  }
  
  // Get radius of scene, etc.
  // Note that this procedure properly takes into account zoom, dolly and pan.
  const G4Point3D targetPoint
    = fSceneHandler.GetScene()->GetStandardTargetPoint()
    + fVP.GetCurrentTargetPoint ();
  G4double radius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
  if(radius<=0.) radius = 1.;
  const G4double cameraDistance = fVP.GetCameraDistance (radius);
  const G4Point3D cameraPosition =
    targetPoint + cameraDistance * fVP.GetViewpointDirection().unit();
  const GLdouble pnear  = fVP.GetNearDistance (cameraDistance, radius);
  const GLdouble pfar   = fVP.GetFarDistance  (cameraDistance, pnear, radius);
  const GLdouble right  = fVP.GetFrontHalfHeight (pnear, radius) * ratioY;
  const GLdouble left   = -right;
  const GLdouble top    = fVP.GetFrontHalfHeight (pnear, radius) * ratioX;
  const GLdouble bottom = -top;
  
  // FIXME
  ResizeGLView();
  //SHOULD SetWindowsSizeHint()...

  glMatrixMode (GL_PROJECTION); // set up Frustum.
  glLoadIdentity();

  const G4Vector3D scaleFactor = fVP.GetScaleFactor();
  glScaled(scaleFactor.x(),scaleFactor.y(),scaleFactor.z());
  
  if (fVP.GetFieldHalfAngle() == 0.) {
    g4GlOrtho (left, right, bottom, top, pnear, pfar);
  }
  else {
    g4GlFrustum (left, right, bottom, top, pnear, pfar);
  }  

  glMatrixMode (GL_MODELVIEW); // apply further transformations to scene.
  glLoadIdentity();
  
  const G4Normal3D& upVector = fVP.GetUpVector ();  
  G4Point3D gltarget;
  if (cameraDistance > 1.e-6 * radius) {
    gltarget = targetPoint;
  }
  else {
    gltarget = targetPoint - radius * fVP.GetViewpointDirection().unit();
  }

  const G4Point3D& pCamera = cameraPosition;  // An alias for brevity.

  g4GluLookAt (pCamera.x(),  pCamera.y(),  pCamera.z(),       // Viewpoint.
             gltarget.x(), gltarget.y(), gltarget.z(),      // Target point.
             upVector.x(), upVector.y(), upVector.z());     // Up vector.
  // Light position is "true" light direction, so must come after gluLookAt.
  glLightfv (GL_LIGHT0, GL_POSITION, lightPosition);

  // The idea is to use back-to-back clipping planes.  This can cut an object
  // down to just a few pixels, which can make it difficult to see.  So, for
  // now, comment this out and use the generic (Boolean) method, via
  // G4VSolid* G4OpenGLSceneHandler::CreateSectionSolid ()
  // { return G4VSceneHandler::CreateSectionSolid(); }
//  if (fVP.IsSection () ) {  // pair of back to back clip planes.
//    const G4Plane3D& sp = fVP.GetSectionPlane ();
//    double sArray[4];
//    sArray[0] = sp.a();
//    sArray[1] = sp.b();
//    sArray[2] = sp.c();
//    sArray[3] = sp.d() + radius * 1.e-05;
//    glClipPlane (GL_CLIP_PLANE0, sArray);
//    glEnable (GL_CLIP_PLANE0);
//    sArray[0] = -sp.a();
//    sArray[1] = -sp.b();
//    sArray[2] = -sp.c();
//    sArray[3] = -sp.d() + radius * 1.e-05;
//    glClipPlane (GL_CLIP_PLANE1, sArray);
//    glEnable (GL_CLIP_PLANE1);
//  } else {
//    glDisable (GL_CLIP_PLANE0);
//    glDisable (GL_CLIP_PLANE1);
//  }

  // What we call intersection of cutaways is easy in OpenGL.  You
  // just keep cutting.  Unions are more tricky - you have to have
  // multiple passes and this is handled in
  // G4OpenGLImmediate/StoredViewer::ProcessView.
  const G4Planes& cutaways = fVP.GetCutawayPlanes();
  size_t nPlanes = cutaways.size();
  if (fVP.IsCutaway() &&
      fVP.GetCutawayMode() == G4ViewParameters::cutawayIntersection &&
      nPlanes > 0) {
    double a[4];
    a[0] = cutaways[0].a();
    a[1] = cutaways[0].b();
    a[2] = cutaways[0].c();
    a[3] = cutaways[0].d();
    glClipPlane (GL_CLIP_PLANE2, a);
    glEnable (GL_CLIP_PLANE2);
    if (nPlanes > 1) {
      a[0] = cutaways[1].a();
      a[1] = cutaways[1].b();
      a[2] = cutaways[1].c();
      a[3] = cutaways[1].d();
      glClipPlane (GL_CLIP_PLANE3, a);
      glEnable (GL_CLIP_PLANE3);
    }
    if (nPlanes > 2) {
      a[0] = cutaways[2].a();
      a[1] = cutaways[2].b();
      a[2] = cutaways[2].c();
      a[3] = cutaways[2].d();
      glClipPlane (GL_CLIP_PLANE4, a);
      glEnable (GL_CLIP_PLANE4);
    }
  } else {
    glDisable (GL_CLIP_PLANE2);
    glDisable (GL_CLIP_PLANE3);
    glDisable (GL_CLIP_PLANE4);
  }

  // Background.
  background = fVP.GetBackgroundColour ();

}



void G4OpenGLViewer::ResetView () {
  G4VViewer::ResetView();
  fRot_sens = 1;
  fPan_sens = 0.01;
}


void G4OpenGLViewer::HaloingFirstPass () {
  
  //To perform haloing, first Draw all information to the depth buffer
  //alone, using a chunky line width, and then Draw all info again, to
  //the colour buffer, setting a thinner line width an the depth testing 
  //function to less than or equal, so if two lines cross, the one 
  //passing behind the other will not pass the depth test, and so not
  //get rendered either side of the infront line for a short distance.

  //First, disable writing to the colo(u)r buffer...
  glColorMask (GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

  //Now enable writing to the depth buffer...
  glDepthMask (GL_TRUE);
  glDepthFunc (GL_LESS);
  glClearDepth (1.0);

  //Finally, set the line width to something wide...
  ChangeLineWidth(3.0);

}

void G4OpenGLViewer::HaloingSecondPass () {

  //And finally, turn the colour buffer back on with a sesible line width...
  glColorMask (GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
  glDepthFunc (GL_LEQUAL);
  ChangeLineWidth(1.0);

}

G4String G4OpenGLViewer::Pick(GLdouble x, GLdouble y)
{
  const std::vector < G4OpenGLViewerPickMap* > & pickMap = GetPickDetails(x,y);
  G4String txt = "";
  if (pickMap.size() == 0) {
//        txt += "No hits recorded.";;
  } else {
    for (unsigned int a=0; a < pickMap.size(); a++) {
      if (pickMap[a]->getAttributes().size() > 0) {
        txt += pickMap[a]->print();
      }
    }
  }
  return txt;
}

const std::vector < G4OpenGLViewerPickMap* > & G4OpenGLViewer::GetPickDetails(GLdouble x, GLdouble y)
{
  static std::vector < G4OpenGLViewerPickMap* > pickMapVector;
  for (auto pickMap: pickMapVector) {
    delete pickMap;
  }
  pickMapVector.clear();
  
  const G4int BUFSIZE = 512;
  GLuint selectBuffer[BUFSIZE];
  glSelectBuffer(BUFSIZE, selectBuffer);
  glRenderMode(GL_SELECT);
  glInitNames();
  glPushName(0);
  glMatrixMode(GL_PROJECTION);
  G4double currentProjectionMatrix[16];
  glGetDoublev(GL_PROJECTION_MATRIX, currentProjectionMatrix);
  glPushMatrix();
  glLoadIdentity();
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
/*  G4cout
  << "viewport, x,y: "
  << viewport[0] << ',' << viewport[1] << ',' << viewport[2] << ',' << viewport[3]
  << ", " << x << ',' << y
  << G4endl;
*/
  fIsGettingPickInfos = true;
  // Define 5x5 pixel pick area
  g4GluPickMatrix(x, viewport[3] - y, 5., 5., viewport);
  glMultMatrixd(currentProjectionMatrix);
  glMatrixMode(GL_MODELVIEW);
  DrawView();
  GLint hits = glRenderMode(GL_RENDER);
  fIsGettingPickInfos = false;
  if (hits < 0) {
    G4cout << "Too many hits.  Zoom in to reduce overlaps." << G4endl;
    goto reset_and_return;
  }
  if (hits > 0) {
    GLuint* p = selectBuffer;
    for (GLint i = 0; i < hits; ++i) {
      G4OpenGLViewerPickMap* pickMap = new G4OpenGLViewerPickMap();
      GLuint nnames = *p++;
      // This bit of debug code or...
      //GLuint zmin = *p++;
      //GLuint zmax = *p++;
      //G4cout << "Hit " << i << ": " << nnames << " names"
      //       << "\nzmin: " << zmin << ", zmax: " << zmax << G4endl;
      // ...just increment the pointer
      p++;
      p++;
      for (GLuint j = 0; j < nnames; ++j) {
        GLuint name = *p++;
        pickMap->setHitNumber(i);
        pickMap->setSubHitNumber(j);
        pickMap->setPickName(name);
        std::map<GLuint, G4AttHolder*>::iterator iter =
        fOpenGLSceneHandler.fPickMap.find(name);
        if (iter != fOpenGLSceneHandler.fPickMap.end()) {
          G4AttHolder* attHolder = iter->second;
          if(attHolder && attHolder->GetAttDefs().size()) {
            for (size_t iAtt = 0; iAtt < attHolder->GetAttDefs().size(); ++iAtt) {
              std::ostringstream oss;
              oss << G4AttCheck(attHolder->GetAttValues()[iAtt],
                                attHolder->GetAttDefs()[iAtt]);
              pickMap->addAttributes(oss.str());
            }
          }
        }
      }
      pickMapVector.push_back(pickMap);
    }
  }

reset_and_return:
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  return pickMapVector;
}

GLubyte* G4OpenGLViewer::grabPixels
(int inColor, unsigned int width, unsigned int height) {
  
  GLubyte* buffer;
  GLint swapbytes, lsbfirst, rowlength;
  GLint skiprows, skippixels, alignment;
  GLenum format;
  int size;

  if (inColor) {
    format = GL_RGB;
    size = width*height*3;
  } else {
    format = GL_LUMINANCE;
    size = width*height*1;
  }

  buffer = new GLubyte[size];
  if (buffer == NULL)
    return NULL;

  glGetIntegerv (GL_UNPACK_SWAP_BYTES, &swapbytes);
  glGetIntegerv (GL_UNPACK_LSB_FIRST, &lsbfirst);
  glGetIntegerv (GL_UNPACK_ROW_LENGTH, &rowlength);

  glGetIntegerv (GL_UNPACK_SKIP_ROWS, &skiprows);
  glGetIntegerv (GL_UNPACK_SKIP_PIXELS, &skippixels);
  glGetIntegerv (GL_UNPACK_ALIGNMENT, &alignment);

  glPixelStorei (GL_UNPACK_SWAP_BYTES, GL_FALSE);
  glPixelStorei (GL_UNPACK_LSB_FIRST, GL_FALSE);
  glPixelStorei (GL_UNPACK_ROW_LENGTH, 0);

  glPixelStorei (GL_UNPACK_SKIP_ROWS, 0);
  glPixelStorei (GL_UNPACK_SKIP_PIXELS, 0);
  glPixelStorei (GL_UNPACK_ALIGNMENT, 1);

  glReadBuffer(GL_FRONT);
  glReadPixels (0, 0, (GLsizei)width, (GLsizei)height, format, GL_UNSIGNED_BYTE, (GLvoid*) buffer);

  glPixelStorei (GL_UNPACK_SWAP_BYTES, swapbytes);
  glPixelStorei (GL_UNPACK_LSB_FIRST, lsbfirst);
  glPixelStorei (GL_UNPACK_ROW_LENGTH, rowlength);
  
  glPixelStorei (GL_UNPACK_SKIP_ROWS, skiprows);
  glPixelStorei (GL_UNPACK_SKIP_PIXELS, skippixels);
  glPixelStorei (GL_UNPACK_ALIGNMENT, alignment);
  
  return buffer;
}

bool G4OpenGLViewer::printEPS() {
  bool res;

  // Change the LC_NUMERIC value in order to have "." separtor and not ","
  // This case is only useful for French, Canadien...
  size_t len = strlen(setlocale(LC_NUMERIC,NULL));
  char* oldLocale = (char*)(malloc(len+1));
  if(oldLocale!=NULL) strncpy(oldLocale,setlocale(LC_NUMERIC,NULL),len);
  setlocale(LC_NUMERIC,"C");

  if (((fExportImageFormat == "eps") || (fExportImageFormat == "ps")) && (!fVectoredPs)) {
    res = printNonVectoredEPS();
  } else {
    res = printVectoredEPS();
  }

  // restore the local
  if (oldLocale) {
    setlocale(LC_NUMERIC,oldLocale);
    free(oldLocale);
  }

  if (res == false) {
    G4cerr << "Error saving file... " << getRealPrintFilename().c_str() << G4endl;
  } else {
    G4cout << "File " << getRealPrintFilename().c_str() << " size: " << getRealExportWidth() << "x" << getRealExportHeight() << " has been saved " << G4endl;

    // increment index if necessary
    if ( fExportFilenameIndex != -1) {
      fExportFilenameIndex++;
    }
  }

  return res;
}

bool G4OpenGLViewer::printVectoredEPS() {
  return printGl2PS();
}

bool G4OpenGLViewer::printNonVectoredEPS () {

  int width = getRealExportWidth();
  int height = getRealExportHeight();

  FILE* fp;
  GLubyte* pixels;
  GLubyte* curpix;
  int components, pos, i;

  pixels = grabPixels (fPrintColour, width, height);

  if (pixels == NULL) {
      G4cerr << "Failed to get pixels from OpenGl viewport" << G4endl;
    return false;
  }
  if (fPrintColour) {
    components = 3;
  } else {
    components = 1;
  }
  std::string name = getRealPrintFilename();
  fp = fopen (name.c_str(), "w");
  if (fp == NULL) {
    G4cerr << "Can't open filename " << name.c_str() << G4endl;
    return false;
  }
  
  fprintf (fp, "%%!PS-Adobe-2.0 EPSF-1.2\n");
  fprintf (fp, "%%%%Title: %s\n", name.c_str());
  fprintf (fp, "%%%%Creator: OpenGL pixmap render output\n");
  fprintf (fp, "%%%%BoundingBox: 0 0 %d %d\n", width, height);
  fprintf (fp, "%%%%EndComments\n");
  fprintf (fp, "gsave\n");
  fprintf (fp, "/bwproc {\n");
  fprintf (fp, "    rgbproc\n");
  fprintf (fp, "    dup length 3 idiv string 0 3 0 \n");
  fprintf (fp, "    5 -1 roll {\n");
  fprintf (fp, "    add 2 1 roll 1 sub dup 0 eq\n");
  fprintf (fp, "    { pop 3 idiv 3 -1 roll dup 4 -1 roll dup\n");
  fprintf (fp, "       3 1 roll 5 -1 roll } put 1 add 3 0 \n");
  fprintf (fp, "    { 2 1 roll } ifelse\n");
  fprintf (fp, "    }forall\n");
  fprintf (fp, "    pop pop pop\n");
  fprintf (fp, "} def\n");
  fprintf (fp, "systemdict /colorimage known not {\n");
  fprintf (fp, "   /colorimage {\n");
  fprintf (fp, "       pop\n");
  fprintf (fp, "       pop\n");
  fprintf (fp, "       /rgbproc exch def\n");
  fprintf (fp, "       { bwproc } image\n");
  fprintf (fp, "   }  def\n");
  fprintf (fp, "} if\n");
  fprintf (fp, "/picstr %d string def\n", width * components);
  fprintf (fp, "%d %d scale\n", width, height);
  fprintf (fp, "%d %d %d\n", width, height, 8);
  fprintf (fp, "[%d 0 0 %d 0 0]\n", width, height);
  fprintf (fp, "{currentfile picstr readhexstring pop}\n");
  fprintf (fp, "false %d\n", components);
  fprintf (fp, "colorimage\n");
  
  curpix = (GLubyte*) pixels;
  pos = 0;
  for (i = width*height*components; i>0; i--) {
    fprintf (fp, "%02hx ", (unsigned short)(*(curpix++)));
    if (++pos >= 32) {
      fprintf (fp, "\n");
      pos = 0; 
    }
  }
  if (pos)
    fprintf (fp, "\n");

  fprintf (fp, "grestore\n");
  fprintf (fp, "showpage\n");
  delete [] pixels;
  fclose (fp);

  // Reset for next time (useful is size change)
  //  fPrintSizeX = -1;
  //  fPrintSizeY = -1;

  return true;
}

/** Return if gl2ps is currently writing
 */
bool G4OpenGLViewer::isGl2psWriting() {

  if (!fGL2PSAction) return false;
  if (fGL2PSAction->fileWritingEnabled()) {
    return true;
  }
  return false;
}


G4bool G4OpenGLViewer::isFramebufferReady() {
  bool check = false;
#ifdef G4VIS_BUILD_OPENGLQT_DRIVER
  check = true;
#endif
#ifdef G4VIS_BUILD_OPENGLX_DRIVER
  check = false;
#endif
#ifdef G4VIS_BUILD_OPENGLXM_DRIVER
  check = false;
#endif
#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER
  check = false;
#endif

#if GL_ARB_framebuffer_object
  if (check) {
//    if ( glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_UNDEFINED) {
//      return false;
//    }
  }
#endif
    return true;
}


/* Draw Gl2Ps text if needed
 */
void G4OpenGLViewer::DrawText(const G4Text& g4text)
{
  // gl2ps or GL window ?
  if (isGl2psWriting()) {

    G4VSceneHandler::MarkerSizeType sizeType;
    G4double size = fSceneHandler.GetMarkerSize(g4text,sizeType);
    G4Point3D position = g4text.GetPosition();

    G4String textString = g4text.GetText();
    
    glRasterPos3d(position.x(),position.y(),position.z());
    GLint align = GL2PS_TEXT_B;

    switch (g4text.GetLayout()) {
    case G4Text::left: align = GL2PS_TEXT_BL; break;
    case G4Text::centre: align = GL2PS_TEXT_B; break;
    case G4Text::right: align = GL2PS_TEXT_BR;
    }
    
    gl2psTextOpt(textString.c_str(),"Times-Roman",GLshort(size),align,0);

  } else {

    static G4int callCount = 0;
    ++callCount;
    //if (callCount <= 10 || callCount%100 == 0) {
    if (callCount <= 1) {
      G4cout <<
	"G4OpenGLViewer::DrawText: Not implemented for \""
	     << fName <<
	"\"\n  Called with "
	     << g4text
	     << G4endl;
    }
  }
}

/** Change PointSize on gl2ps if needed
 */
void G4OpenGLViewer::ChangePointSize(G4double size) {

  if (isGl2psWriting()) {
    fGL2PSAction->setPointSize(int(size));
  } else {
    glPointSize (size);
  }
}


/** Change LineSize on gl2ps if needed
 */
void G4OpenGLViewer::ChangeLineWidth(G4double width) {

  if (isGl2psWriting()) {
    fGL2PSAction->setLineWidth(int(width));
  } else {
    glLineWidth (width);
  }
}

/**
 Export image with the given name with width and height
 Several cases :
 If name is "", filename will have the default value
 If name is "toto.png", set the name to "toto" and the format to "png". No incremented suffix is added.
 If name is "toto", set the name to "toto" and the format to default (or current format if specify).
   Will also add an incremented suffix at the end of the file
*/
bool G4OpenGLViewer::exportImage(std::string name, int width, int height) {

  if (! setExportFilename(name)) {
    return false;
  }

  if ((width != -1) && (height != -1)) {
    setExportSize(width, height);
  }

  if (fExportImageFormat == "eps") {
    fGL2PSAction->setExportImageFormat(GL2PS_EPS);
  } else if (fExportImageFormat == "ps") {
    fGL2PSAction->setExportImageFormat(GL2PS_PS);
  } else if (fExportImageFormat == "svg") {
    fGL2PSAction->setExportImageFormat(GL2PS_SVG);
  } else if (fExportImageFormat == "pdf") {
    fGL2PSAction->setExportImageFormat(GL2PS_PDF);
  } else {
    setExportImageFormat(fExportImageFormat,true); // will display a message if this format is not correct for the current viewer
    return false;
  }
  return printEPS();
}


bool G4OpenGLViewer::printGl2PS() {

  int width = getRealExportWidth();
  int height = getRealExportHeight();
  bool res = true;

  // no need to redraw at each new primitive for printgl2PS
  G4OpenGLSceneHandler& oglSceneHandler = dynamic_cast<G4OpenGLSceneHandler&>(fSceneHandler);
  G4OpenGLSceneHandler::FlushAction originalFlushAction = oglSceneHandler.GetFlushAction();
  oglSceneHandler.SetFlushAction(G4OpenGLSceneHandler::never);

  if (!fGL2PSAction) return false;

  fGL2PSAction->setFileName(getRealPrintFilename().c_str());
  // try to resize
  int X = fWinSize_x;
  int Y = fWinSize_y;

  fWinSize_x = width;
  fWinSize_y = height;
  // Laurent G. 16/03/10 : Not the good way to do. 
  // We should draw in a new offscreen context instead of
  // resizing and drawing in current window...
  // This should be solve when we will do an offscreen method
  // to render OpenGL
  // See : 
  // http://developer.apple.com/Mac/library/documentation/GraphicsImaging/Conceptual/OpenGL-MacProgGuide/opengl_offscreen/opengl_offscreen.html
  // http://www.songho.ca/opengl/gl_fbo.html

   ResizeGLView();
   bool extendBuffer = true;
   bool endWriteAction = false;
   bool beginWriteAction = true;
   bool filePointerOk = true;
   while ((extendBuffer) && (! endWriteAction) && (filePointerOk)) {
     
     beginWriteAction = fGL2PSAction->enableFileWriting();
     // 3 cases :
     // - true
     // - false && ! fGL2PSAction->fileWritingEnabled() => bad file name
     // - false && fGL2PSAction->fileWritingEnabled() => buffer size problem ?

     filePointerOk = fGL2PSAction->fileWritingEnabled();
       
     if (beginWriteAction) {
       
       // Set the viewport
       // By default, we choose the line width (trajectories...)
       fGL2PSAction->setLineWidth(fGl2psDefaultLineWith);
       // By default, we choose the point size (markers...)
       fGL2PSAction->setPointSize(fGl2psDefaultPointSize);
       
       DrawView ();
       endWriteAction = fGL2PSAction->disableFileWriting();
     }
     if (filePointerOk) {
       if ((! endWriteAction) || (! beginWriteAction)) {
         extendBuffer = fGL2PSAction->extendBufferSize();
       }
     }
   }
   fGL2PSAction->resetBufferSizeParameters();

   if (!extendBuffer ) {
     G4cerr << "ERROR: gl2ps buffer size is not big enough to print this geometry. Try to extend it. No output produced"<< G4endl;
     res = false;
   }
   if (!beginWriteAction ) {
     G4cerr << "ERROR: saving file "<<getRealPrintFilename().c_str()<<". Check read/write access. No output produced" << G4endl;
     res = false;
   }
   if (!endWriteAction ) {
     G4cerr << "gl2ps error. No output produced" << G4endl;
     res = false;
   }
  fWinSize_x = X;
  fWinSize_y = Y;

  oglSceneHandler.SetFlushAction(originalFlushAction);

  // Reset for next time (useful is size change)
  //  fPrintSizeX = 0;
  //  fPrintSizeY = 0;

  return res;
}

unsigned int G4OpenGLViewer::getWinWidth() const{
  return fWinSize_x;
}

unsigned int G4OpenGLViewer::getWinHeight() const{
  return fWinSize_y;
}

G4bool G4OpenGLViewer::sizeHasChanged() {
  return fSizeHasChanged;
}

G4int G4OpenGLViewer::getRealExportWidth() {
  if (fPrintSizeX == -1) {
    return fWinSize_x;
  }
  GLint dims[2];
  glGetIntegerv(GL_MAX_VIEWPORT_DIMS, dims);

  // L.Garnier 01-2010: Some problems with mac 10.6
  if ((dims[0] !=0 ) && (dims[1] !=0)) {
    if (fPrintSizeX > dims[0]){
      return dims[0];
    }
  }
  if (fPrintSizeX < -1){
    return 0;
  }
  return fPrintSizeX;
}

G4int G4OpenGLViewer::getRealExportHeight() {
  if (fPrintSizeY == -1) {
    return fWinSize_y;
  }
  GLint dims[2];
  glGetIntegerv(GL_MAX_VIEWPORT_DIMS, dims);

  // L.Garnier 01-2010: Some problems with mac 10.6
  if ((dims[0] !=0 ) && (dims[1] !=0)) {
    if (fPrintSizeY > dims[1]){
      return dims[1];
    }
  }
  if (fPrintSizeY < -1){
    return 0;
  }
  return fPrintSizeY;
}

void G4OpenGLViewer::setExportSize(G4int X, G4int Y) {
  fPrintSizeX = X;
  fPrintSizeY = Y;
}

/**
 If name is "" or "!", filename and extension will have the default value.
 If name is "toto.png", set the name to "toto" and the format to "png". No incremented suffix is added.
 If name is "toto", set the name to "toto" and the format to default (or current format if specify).
 If name is the same as previous, do not reset incremented suffix.
*/
bool G4OpenGLViewer::setExportFilename(G4String name,G4bool inc) {
  if (name == "!") {
    name = "";
  }

  if (inc) {
    if ((name != "") && (fExportFilename != name)) {
      fExportFilenameIndex=0;
    }
  } else {
    fExportFilenameIndex=-1;
  }

  if (name.size() == 0) {
    name = getRealPrintFilename().c_str();
  } else {
    // guess format by extention
    std::string extension = name.substr(name.find_last_of(".") + 1);
    // no format
    if (name.size() != extension.size()) {
      if (! setExportImageFormat(extension, false)) {
        return false;
      }
    }
    // get the name
    fExportFilename = name.substr(0,name.find_last_of("."));
  }
  return true;
}

std::string G4OpenGLViewer::getRealPrintFilename() {
  std::string temp = fExportFilename;
  if (fExportFilenameIndex != -1) {
    temp += std::string("_");
    std::ostringstream os;
    os << std::setw(4) << std::setfill('0') << fExportFilenameIndex;
    std::string nb_str = os.str();
    temp += nb_str;
  }
  temp += "."+fExportImageFormat;
  return temp;
}

GLdouble G4OpenGLViewer::getSceneNearWidth()
{
  if (!fSceneHandler.GetScene()) {
    return 0;
  }
  const G4Point3D targetPoint
    = fSceneHandler.GetScene()->GetStandardTargetPoint()
    + fVP.GetCurrentTargetPoint ();
  G4double radius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
  if(radius<=0.) radius = 1.;
  const G4double cameraDistance = fVP.GetCameraDistance (radius);
  const GLdouble pnear   = fVP.GetNearDistance (cameraDistance, radius);
  return 2 * fVP.GetFrontHalfHeight (pnear, radius);
}

GLdouble G4OpenGLViewer::getSceneFarWidth()
{
  if (!fSceneHandler.GetScene()) {
    return 0;
  }
  const G4Point3D targetPoint
    = fSceneHandler.GetScene()->GetStandardTargetPoint()
    + fVP.GetCurrentTargetPoint ();
  G4double radius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
  if(radius<=0.) radius = 1.;
  const G4double cameraDistance = fVP.GetCameraDistance (radius);
  const GLdouble pnear   = fVP.GetNearDistance (cameraDistance, radius);
  const GLdouble pfar    = fVP.GetFarDistance  (cameraDistance, pnear, radius);
  return 2 * fVP.GetFrontHalfHeight (pfar, radius);
}


GLdouble G4OpenGLViewer::getSceneDepth()
{
  if (!fSceneHandler.GetScene()) {
    return 0;
  }
  const G4Point3D targetPoint
    = fSceneHandler.GetScene()->GetStandardTargetPoint()
    + fVP.GetCurrentTargetPoint ();
  G4double radius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
  if(radius<=0.) radius = 1.;
  const G4double cameraDistance = fVP.GetCameraDistance (radius);
  const GLdouble pnear   = fVP.GetNearDistance (cameraDistance, radius);
  return fVP.GetFarDistance  (cameraDistance, pnear, radius)- pnear;
}



void G4OpenGLViewer::rotateScene(G4double dx, G4double dy)
{
  if (fVP.GetRotationStyle() == G4ViewParameters::freeRotation) {
    rotateSceneInViewDirection(dx,dy);
  } else {
    if( dx != 0) {
      rotateSceneThetaPhi(dx,0);
    }
    if( dy != 0) {
      rotateSceneThetaPhi(0,dy);
    }
  }
}


void G4OpenGLViewer::rotateSceneToggle(G4double dx, G4double dy)
{
  if (fVP.GetRotationStyle() != G4ViewParameters::freeRotation) {
    rotateSceneInViewDirection(dx,dy);
  } else {
    if( dx != 0) {
      rotateSceneThetaPhi(dx,0);
    }
    if( dy != 0) {
      rotateSceneThetaPhi(0,dy);
    }
  }
}

void G4OpenGLViewer::rotateSceneThetaPhi(G4double dx, G4double dy)
{
  if (!fSceneHandler.GetScene()) {
    return;
  }

  G4Vector3D vp;
  G4Vector3D up;
  
  G4Vector3D xprime;
  G4Vector3D yprime;
  G4Vector3D zprime;
  
  G4double delta_alpha;
  G4double delta_theta;
  
  G4Vector3D new_vp;
  G4Vector3D new_up;
  
  G4double cosalpha;
  G4double sinalpha;
  
  G4Vector3D a1;
  G4Vector3D a2;
  G4Vector3D delta;
  G4Vector3D viewPoint;

    
  //phi spin stuff here
  
  vp = fVP.GetViewpointDirection ().unit ();
  up = fVP.GetUpVector ().unit ();
  
  yprime = (up.cross(vp)).unit();
  zprime = (vp.cross(yprime)).unit();
  
  if (fVP.GetLightsMoveWithCamera()) {
    delta_alpha = dy * fRot_sens;
    delta_theta = -dx * fRot_sens;
  } else {
    delta_alpha = -dy * fRot_sens;
    delta_theta = dx * fRot_sens;
  }    
  
  delta_alpha *= deg;
  delta_theta *= deg;
  
  new_vp = std::cos(delta_alpha) * vp + std::sin(delta_alpha) * zprime;
  
  // to avoid z rotation flipping
  // to allow more than 360∞ rotation

  if (fVP.GetLightsMoveWithCamera()) {
    new_up = (new_vp.cross(yprime)).unit();
    if (new_vp.z()*vp.z() <0) {
      new_up.set(new_up.x(),-new_up.y(),new_up.z());
    }
  } else {
    new_up = up;
    if (new_vp.z()*vp.z() <0) {
      new_up.set(new_up.x(),-new_up.y(),new_up.z());
    }
  }
  fVP.SetUpVector(new_up);
  ////////////////
  // Rotates by fixed azimuthal angle delta_theta.
  
  cosalpha = new_up.dot (new_vp.unit());
  sinalpha = std::sqrt (1. - std::pow (cosalpha, 2));
  yprime = (new_up.cross (new_vp.unit())).unit ();
  xprime = yprime.cross (new_up);
  // Projection of vp on plane perpendicular to up...
  a1 = sinalpha * xprime;
  // Required new projection...
  a2 = sinalpha * (std::cos (delta_theta) * xprime + std::sin (delta_theta) * yprime);
  // Required Increment vector...
  delta = a2 - a1;
  // So new viewpoint is...
  viewPoint = new_vp.unit() + delta;
  
  fVP.SetViewAndLights (viewPoint);
}


void G4OpenGLViewer::rotateSceneInViewDirection(G4double dx, G4double dy)
{
  if (!fSceneHandler.GetScene()) {
    return;
  }

  G4Vector3D vp;
  G4Vector3D up;
  
  G4Vector3D xprime;
  G4Vector3D yprime;
  G4Vector3D zprime;
  
  G4Vector3D new_vp;
  G4Vector3D new_up;
  
  G4Vector3D a1;
  G4Vector3D a2;
  G4Vector3D delta;
  G4Vector3D viewPoint;

  dx = dx/100;
  dy = dy/100;

  //phi spin stuff here
  
  vp = fVP.GetViewpointDirection ().unit();
  up = fVP.GetUpVector ().unit();

  G4Vector3D zPrimeVector = G4Vector3D(up.y()*vp.z()-up.z()*vp.y(),
                             up.z()*vp.x()-up.x()*vp.z(),
                             up.x()*vp.y()-up.y()*vp.x());

  viewPoint = vp/fRot_sens + (zPrimeVector*dx - up*dy) ;
  new_up = G4Vector3D(viewPoint.y()*zPrimeVector.z()-viewPoint.z()*zPrimeVector.y(),
                       viewPoint.z()*zPrimeVector.x()-viewPoint.x()*zPrimeVector.z(),
                       viewPoint.x()*zPrimeVector.y()-viewPoint.y()*zPrimeVector.x());

  G4Vector3D new_upUnit = new_up.unit();
  
  

   fVP.SetUpVector(new_upUnit);
   fVP.SetViewAndLights (viewPoint);
}


void G4OpenGLViewer::addExportImageFormat(std::string format) {
  fExportImageFormatVector.push_back(format);
}

bool G4OpenGLViewer::setExportImageFormat(std::string format, bool quiet) {
  bool found = false;
  std::string list;
  for (unsigned int a=0; a<fExportImageFormatVector.size(); a++) {
    list +=fExportImageFormatVector.at(a) + " ";

    if (fExportImageFormatVector.at(a) == format) {
      if (! quiet) {
        G4cout << " Changing export format to \"" << format << "\"" << G4endl;
      }
      if (format != fExportImageFormat) {
        fExportFilenameIndex = 0;
        fExportImageFormat = format;
      }
      return true;
    }
  }
  if (! found) {
    if (format.size() == 0) {
      G4cout << " Current formats availables are : " << list << G4endl;
    } else {
      G4cerr << " Format \"" << format << "\" is not available for the selected viewer. Current formats availables are : " << list << G4endl;
    }
  }
  return false;
}


// From MESA implementation :
// http://www.techques.com/question/1-8660454/gluPickMatrix-code-from-Mesa

void G4OpenGLViewer::g4GluPickMatrix(GLdouble x, GLdouble y, GLdouble width, GLdouble height,
                     GLint viewport[4])
  {
    GLdouble mat[16];
    GLdouble sx, sy;
    GLdouble tx, ty;
    
    sx = viewport[2] / width;
    sy = viewport[3] / height;
    tx = (viewport[2] + 2.0 * (viewport[0] - x)) / width;
    ty = (viewport[3] + 2.0 * (viewport[1] - y)) / height;
    
#define M(row, col) mat[col*4+row]
    M(0, 0) = sx;
    M(0, 1) = 0.0;
    M(0, 2) = 0.0;
    M(0, 3) = tx;
    M(1, 0) = 0.0;
    M(1, 1) = sy;
    M(1, 2) = 0.0;
    M(1, 3) = ty;
    M(2, 0) = 0.0;
    M(2, 1) = 0.0;
    M(2, 2) = 1.0;
    M(2, 3) = 0.0;
    M(3, 0) = 0.0;
    M(3, 1) = 0.0;
    M(3, 2) = 0.0;
    M(3, 3) = 1.0;
#undef M
    
    glMultMatrixd(mat);
}





// From MESA implementation :
// https://github.com/jlamarche/iOS-OpenGLES-Stuff/blob/master/Wavefront%20OBJ%20Loader/Classes/gluLookAt.m
// or http://www.daniweb.com/software-development/game-development/threads/308901/lookat-matrix-source-code

void G4OpenGLViewer::g4GluLookAt( GLdouble eyex, GLdouble eyey, GLdouble eyez,
                        GLdouble centerx, GLdouble centery, GLdouble
                        centerz,
                        GLdouble upx, GLdouble upy, GLdouble upz )
{
	GLdouble mat[16];
	GLdouble x[3], y[3], z[3];
	GLdouble mag;
  
	/* Make rotation matrix */
  
	/* Z vector */
	z[0] = eyex - centerx;
	z[1] = eyey - centery;
	z[2] = eyez - centerz;
	mag = std::sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2]);
	if (mag) {			/* mpichler, 19950515 */
		z[0] /= mag;
		z[1] /= mag;
		z[2] /= mag;
	}
  
	/* Y vector */
	y[0] = upx;
	y[1] = upy;
	y[2] = upz;
  
	/* X vector = Y cross Z */
	x[0] = y[1] * z[2] - y[2] * z[1];
	x[1] = -y[0] * z[2] + y[2] * z[0];
	x[2] = y[0] * z[1] - y[1] * z[0];
  
	/* Recompute Y = Z cross X */
	y[0] = z[1] * x[2] - z[2] * x[1];
	y[1] = -z[0] * x[2] + z[2] * x[0];
	y[2] = z[0] * x[1] - z[1] * x[0];
  
	/* mpichler, 19950515 */
	/* cross product gives area of parallelogram, which is < 1.0 for
	 * non-perpendicular unit-length vectors; so normalize x, y here
	 */
  
	mag = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
	if (mag) {
		x[0] /= mag;
		x[1] /= mag;
		x[2] /= mag;
	}
  
	mag = std::sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
	if (mag) {
		y[0] /= mag;
		y[1] /= mag;
		y[2] /= mag;
	}
  
#define M(row,col)  mat[col*4+row]
	M(0, 0) = x[0];
	M(0, 1) = x[1];
	M(0, 2) = x[2];
	M(0, 3) = 0.0;
	M(1, 0) = y[0];
	M(1, 1) = y[1];
	M(1, 2) = y[2];
	M(1, 3) = 0.0;
	M(2, 0) = z[0];
	M(2, 1) = z[1];
	M(2, 2) = z[2];
	M(2, 3) = 0.0;
	M(3, 0) = 0.0;
	M(3, 1) = 0.0;
	M(3, 2) = 0.0;
	M(3, 3) = 1.0;
#undef M
	glMultMatrixd(mat);
  
	/* Translate Eye to Origin */
	glTranslated(-eyex, -eyey, -eyez);
}

void G4OpenGLViewer::g4GlOrtho (GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble zNear, GLdouble zFar) {
  //  glOrtho (left, right, bottom, top, near, far);
  
  GLdouble a = 2.0 / (right - left);
  GLdouble b = 2.0 / (top - bottom);
  GLdouble c = -2.0 / (zFar - zNear);
  
  GLdouble tx = - (right + left)/(right - left);
  GLdouble ty = - (top + bottom)/(top - bottom);
  GLdouble tz = - (zFar + zNear)/(zFar - zNear);
  
  GLdouble ortho[16] = {
    a, 0, 0, 0,
    0, b, 0, 0,
    0, 0, c, 0,
    tx, ty, tz, 1
  };
  glMultMatrixd(ortho);
  
}


void G4OpenGLViewer::g4GlFrustum (GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble zNear, GLdouble zFar) {
  //  glFrustum (left, right, bottom, top, near, far);
  
  GLdouble deltaX = right - left;
  GLdouble deltaY = top - bottom;
  GLdouble deltaZ = zFar - zNear;
  
  GLdouble a = 2.0f * zNear / deltaX;
  GLdouble b = 2.0f * zNear / deltaY;
  GLdouble c = (right + left) / deltaX;
  GLdouble d = (top + bottom) / deltaY;
  GLdouble e = -(zFar + zNear) / (zFar - zNear);
  GLdouble f = -2.0f * zFar * zNear / deltaZ;
  
  GLdouble proj[16] = {
    a, 0, 0, 0,
    0, b, 0, 0,
    c, d, e, -1.0f,
    0, 0, f, 0
  };
  
  glMultMatrixd(proj);
  
}


#ifdef G4OPENGL_VERSION_2

// Associate the VBO drawer to the OpenGLViewer and the OpenGLSceneHandler
void G4OpenGLViewer::setVboDrawer(G4OpenGLVboDrawer* drawer) {
  fVboDrawer = drawer;
  try {
    G4OpenGLSceneHandler& sh = dynamic_cast<G4OpenGLSceneHandler&>(fSceneHandler);
    sh.setVboDrawer(fVboDrawer);
  } catch(std::bad_cast exp) { }
}

#endif


G4String G4OpenGLViewerPickMap::print() {
  std::ostringstream txt;

  txt << fName;

  txt << "Hit: " << fHitNumber << ", Sub-hit: " << fSubHitNumber << ", PickName: " << fPickName << "\n";
  
  for (unsigned int a=0; a<fAttributes.size(); a++) {
    txt << fAttributes[a] << "\n";
  }
  return txt.str();
}

#endif
