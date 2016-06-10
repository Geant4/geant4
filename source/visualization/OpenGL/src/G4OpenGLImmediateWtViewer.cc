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
// $Id: G4OpenGLImmediateWtViewer.cc 85263 2014-10-27 08:58:31Z gcosmo $
//
//
// Class G4OpenGLImmediateWtViewer : a class derived from G4OpenGLWtViewer and
//                                G4OpenGLImmediateViewer.

#ifdef G4VIS_BUILD_OPENGLWT_DRIVER

#include "G4OpenGLImmediateWtViewer.hh"
#include "G4OpenGLImmediateSceneHandler.hh"

#include "G4ios.hh"
#define G4DEBUG_VIS_OGL 1

G4OpenGLImmediateWtViewer::G4OpenGLImmediateWtViewer
(G4OpenGLImmediateSceneHandler& sceneHandler,
  Wt::WContainerWidget* aParent,
  const G4String&  name):
  G4VViewer (sceneHandler, sceneHandler.IncrementViewCount (), name),
  G4OpenGLViewer (sceneHandler),
  G4OpenGLWtViewer (sceneHandler),
  G4OpenGLImmediateViewer (sceneHandler),
  Wt::WGLWidget(aParent)

{
// Create a new drawer
  // register the WtDrawer to the OpenGLViewer
  setVboDrawer(new G4OpenGLVboDrawer(this,"OGL-ES"));

  // Add the GL Widget to its parent
  aParent->addWidget(this);

  fHasToRepaint = false;
  fIsRepainting = false;

#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateWtViewer INIT\n");
#endif
  
  
  if (fViewId < 0) return;  // In case error in base class instantiation.
}

G4OpenGLImmediateWtViewer::~G4OpenGLImmediateWtViewer() {
}

void G4OpenGLImmediateWtViewer::Initialise() {
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateWtViewer::Initialise \n");
#endif
  fReadyToPaint = false;
  CreateMainWindow (this,Wt::WString(fName));
  CreateFontLists ();

  fReadyToPaint = true;
}

void G4OpenGLImmediateWtViewer::initializeGL () {

  InitializeGLView ();

  // If a double buffer context has been forced upon us, ignore the
  // back buffer for this OpenGLImmediate view.
  //  glDrawBuffer (GL_FRONT); // FIXME : Ne marche pas avec cette ligne, mais affiche le run correctement...

  if (fSceneHandler.GetScene() == 0) {
    fHasToRepaint =false;
  } else {
    fHasToRepaint =true;
  }

  // In order to know where to look at, calculate the centerpoint of the
  // scene
  double cx, cy, cz;
  cx = cy = cz = 0.;
  
  // Transform the world so that we look at the centerpoint of the scene
  Wt::WMatrix4x4 worldTransform;
  worldTransform.lookAt(
                        cx, cy, cz + 10, // camera position
                        cx, cy, cz,      // looking at
                        0, 1, 0);        // 'up' vector
  
  // We want to be able to change the camera position client-side. In
  // order to do so, the world transformation matrix must be stored in
  // a matrix that can be manipulated from JavaScript.
  jsMatrix_ = createJavaScriptMatrix4();
  setJavaScriptMatrix4(jsMatrix_, worldTransform);
  
  // This installs a client-side mouse handler that modifies the
  // world transformation matrix. Like WMatrix4x4::lookAt, this works
  // by specifying a center point and an up direction; mouse movements
  // will allow the camera to be moved around the center point.
  setClientSideLookAtHandler(jsMatrix_, // the name of the JS matrix
                             cx, cy, cz,                       // the center point
                             0, 1, 0,                          // the up direction
                             0.005, 0.005);                    // 'speed' factors
  // Alternative: this installs a client-side mouse handler that allows
  // to 'walk' around: go forward, backward, turn left, turn right, ...
  //setClientSideWalkHandler(jsMatrix_, 0.05, 0.005);

  
  // Set the clear color to a transparant background
  glClearColor(0, 0, 0, 0);
  
  // Reset Z-buffer, enable Z-buffering
  glClearDepth(1);
  glEnable(DEPTH_TEST);
  glDepthFunc(LEQUAL);
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLWtViewer initializeGL END\n");
#endif
  
  // Set the component visible
  show() ;
  
}

void  G4OpenGLImmediateWtViewer::DrawView() {
  updateWWidget();
}


void G4OpenGLImmediateWtViewer::ComputeView () {
  
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLWtViewer::ComputeView %d %d   VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV\n",getWinWidth(), getWinHeight());
#endif
  
  // If a double buffer context has been forced upon us, ignore the
  // back buffer for this OpenGL view.
  //  glDrawBuffer (GL_FRONT);

  G4ViewParameters::DrawingStyle dstyle = GetViewParameters().GetDrawingStyle();

  if(dstyle!=G4ViewParameters::hlr &&
     haloing_enabled) {

    HaloingFirstPass ();
    NeedKernelVisit ();
    ProcessView ();
    FinishView();
#ifdef G4DEBUG_VIS_OGL
    printf("G4OpenGLWtViewer::ComputeView First ProcessView ok\n");
#endif
    HaloingSecondPass ();

  }

  NeedKernelVisit ();  // Always need to visit G4 kernel.
  ProcessView ();

/* FIXME
 if (isRecording()) {
    savePPMToTemp();
  }
 */
  
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLWtViewer::ComputeView %d %d ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n",getWinWidth(), getWinHeight());
#endif
  fHasToRepaint = true;
}

/**
 - Lors du resize de la fenetre, on doit non pas redessiner le detecteur, mais aussi les evenements
 */
void G4OpenGLImmediateWtViewer::resizeGL(
 int width
,int height)
{
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateWtViewer resizeGL %d %d\n",width,height);
#endif
  G4OpenGLWtViewer::resizeGL(width,height);

  // Set the viewport size.
  glViewport(0, 0, width, height);

  // Set projection matrix to some fixed values
/*  Wt::WMatrix4x4 proj;
  proj.perspective(45, ((double)width)/height, 1, 40);
  glUniformMatrix4(fpMatrixUniform, proj);
*/
  SetView();
  //  updateWWidget();
}


void G4OpenGLImmediateWtViewer::paintGL() {

  if (fIsRepainting) {
    //    return ;
  }
  if (!fReadyToPaint) {
    fReadyToPaint= true;
    return;
  }
  if ((getWinWidth() == 0) && (getWinHeight() == 0)) {
      return;
  }

  // DO NOT RESIZE IF SIZE HAS NOT CHANGE
  if ( !fHasToRepaint) {
    double sw = 0;
    double sh = 0;
    //    if (!isMaximized() && !isFullScreen()) {
    sw = width().value();
    sh = height().value();
    //     } else {
    //       sw = frameGeometry().width();
    //       sh = frameGeometry().height();
    //     }
    if ((getWinWidth() == sw) &&(getWinHeight() == sh)) {
      return;

    } else if ((sw == 0) && (sh == 0)) { // NOT A TOP LEVEL WIDGET
      if (((getWinWidth() == width().value())) &&(getWinHeight() == height().value())) {
        return;
      }
    }
  }
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateWtViewer::paintGL VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV ready %d\n",fReadyToPaint);
#endif
  
  SetView();

  ClearView(); //ok, put the background correct
  
  // SHOULD BE REMOVED

  // Configure the shader: set the uniforms
  // Uniforms are 'configurable constants' in a shader: they are
  // identical for every point that has to be drawn.
  // Set the camera transformation to the value of a client-side JS matrix
  glUniformMatrix4(fcMatrixUniform, jsMatrix_);
  // Often, a model matrix is used to move the model around. We're happy
  // with the location of the model, so we leave it as the unit matrix
  Wt::WMatrix4x4 modelMatrix;
  glUniformMatrix4(fmvMatrixUniform, modelMatrix);
  // The next one is a bit complicated. In desktop OpenGL, a shader
  // has the gl_NormalMatrix matrix available in the shader language,
  // a matrix that is used to transform normals to e.g. implement proper
  // Phong shading (google will help you to find a detailed explanation
  // of why you need it). It is the transposed inverse of the model view
  // matrix. Unfortunately, this matrix is not available in WebGL, so if
  // you want to do phong shading, you must calculate it yourself.
  // Wt provides methods to calculate the transposed inverse of a matrix,
  // when client-side JS matrices are involved. Here, we inverse-transpose
  // the product of the client-side camera matrix and the model matrix.
  glUniformMatrix4(fnMatrixUniform, (jsMatrix_ * modelMatrix).inverted().transposed());

  // Create a new Buffer
  Buffer objBuffer_2 = glCreateBuffer(); //glGenBuffers(1,&objBuffer_2)
  
  // Bind this buffer
  glBindBuffer(GL_ARRAY_BUFFER, objBuffer_2);
  // SHOULD BE REMOVED END


  ComputeView();

  fHasToRepaint = false; // could be set to false by ComputeView

#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateQtViewer::paintGL ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ready %d\n\n\n",fReadyToPaint);
#endif
  fIsRepainting = false;

}


void G4OpenGLImmediateWtViewer::mousePressEvent(Wt::WMouseEvent *event)
{
  // boutons et position
  G4MousePressEvent(event);
}

void G4OpenGLImmediateWtViewer::keyPressEvent (Wt::WKeyEvent *event)
{
  G4keyPressEvent(event);
}

/**
 void G4OpenGLImmediateWtViewer::showEvent (QShowEvent *)
 {
 fHasToRepaint = true;
 }
 */


/**
 * This function was build in order to make a zoom on double clic event.
 * It was think to build a rubberband on the zoom area, but never work fine
 */
void G4OpenGLImmediateWtViewer::mouseDoubleClickEvent(Wt::WMouseEvent *)
{
  G4MouseDoubleClickEvent();
}


void G4OpenGLImmediateWtViewer::mouseReleaseEvent(Wt::WMouseEvent )
{
  G4MouseReleaseEvent();
}


void G4OpenGLImmediateWtViewer::mouseMoveEvent(Wt::WMouseEvent *event)
{
  G4MouseMoveEvent(event);
}


/**
 void G4OpenGLImmediateWtViewer::contextMenuEvent(QContextMenuEvent *e)
 {
 G4manageContextMenuEvent(e);
 }
 */

void G4OpenGLImmediateWtViewer::paintEvent(Wt::WPaintDevice * /* painter */) {
  if ( fHasToRepaint) {
    updateGL();
  }
}



void G4OpenGLImmediateWtViewer::FinishView()
{
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLWtViewer::FinishView() \n");
#endif
  flush ();
  
  // L. Garnier 10/2009 : Not necessary and cause problems on mac OS X 10.6
  //  fWindow->swapBuffers ();
}



void G4OpenGLImmediateWtViewer::popMatrix() {
}

void G4OpenGLImmediateWtViewer::pushMatrix() {
}

void G4OpenGLImmediateWtViewer::multMatrixd(const GLdouble* /* m */) {
  //  mMatrix = mMatrix * m;
}

void G4OpenGLImmediateWtViewer::loadIdentity() {
  mMatrix.setToIdentity ();
}


void G4OpenGLImmediateWtViewer::setMatrixUniforms() {
  /*
   UniformLocation  pUniform = getUniformLocation(shaderProgram, "uPMatrix");
   uniformMatrix4fv(pUniform, false, new Float32Array(perspectiveMatrix.flatten()));
   
   UniformLocation  mvUniform = getUniformLocation(shaderProgram, "uMVMatrix");
   uniformMatrix4fv(mvUniform, false, new Float32Array(mvMatrix.flatten()));
   */
}


void G4OpenGLImmediateWtViewer::ShowView (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateWtViewer ShowView\n");
#endif
  repaintSlot();
  //  activateWindow();
}


void G4OpenGLImmediateWtViewer::updateWWidget() {
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateWtViewer updateWWidget\n");
#endif
  fHasToRepaint= true;
  //  updateGL();
  repaintGL(PAINT_GL | RESIZE_GL);
  //  paintGL() ;
  fHasToRepaint= false;
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateWtViewer updateWWidget END\n");
#endif
}


#endif
