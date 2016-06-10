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
// $Id: G4OpenGLImmediateWtViewer.cc 75567 2013-11-04 11:35:11Z gcosmo $
//
//
// Class G4OpenGLImmediateWtViewer : a class derived from G4OpenGLWtViewer and
//                                G4OpenGLImmediateViewer.

#ifdef G4VIS_BUILD_OPENGLWT_DRIVER

#include "G4OpenGLImmediateWtViewer.hh"
#include "G4OpenGLImmediateSceneHandler.hh"

#include "G4ios.hh"

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
  fWtDrawer = new G4OpenGLWtDrawer(this);
  
  // register the WtDrawer to the OpenGLViewer
  setWtDrawer(fWtDrawer);

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
  Wt::WMatrix4x4 proj;
  proj.perspective(45, ((double)width)/height, 1, 40);
  glUniformMatrix4(pMatrixUniform_, proj);
  //  updateWWidget();
}


void G4OpenGLImmediateWtViewer::paintGL() {
// FIXME : for test
  glClearColor(0.1234, 0.1234, 0.1234, 0.1234);

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
  drawScene();

  ComputeView();

  fHasToRepaint = false; // could be set to false by ComputeView

#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateQtViewer::paintGL ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ready %d\n\n\n",fReadyToPaint);
#endif
  fIsRepainting = false;
  // FIXME : for test
  glClearColor(0.1234, 0.1234, 0.1234, 0.1234);
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


void G4OpenGLImmediateWtViewer::SetView () {
  
  if (!fSceneHandler.GetScene()) {
    return;
  }
  // Calculates view representation based on extent of object being
  // viewed and (initial) viewpoint.  (Note: it can change later due
  // to user interaction via visualization system's GUI.)
  
  /*
  // Lighting.
  GLfloat lightPosition [4];
  lightPosition [0] = fVP.GetActualLightpointDirection().x();
  lightPosition [1] = fVP.GetActualLightpointDirection().y();
  lightPosition [2] = fVP.GetActualLightpointDirection().z();
  lightPosition [3] = 0.;
  // Light position is "true" light direction, so must come after gluLookAt.
  GLfloat ambient [] = { 0.2, 0.2, 0.2, 1.};
  GLfloat diffuse [] = { 0.8, 0.8, 0.8, 1.};
//  enable (Wt::WGLWidget::LIGHT0);
  
 G4double ratioX = 1;
  G4double ratioY = 1;
  if (getWinHeight() > getWinWidth()) {
    ratioX = ((G4double)getWinHeight()) / ((G4double)getWinWidth());
  }
  if (getWinWidth() > getWinHeight()) {
    ratioY = ((G4double)getWinWidth()) / ((G4double)getWinHeight());
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
  */
  // FIXME
//  ResizeGLView();
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLWtViewer::SetView() resize viewport to %d %d\n",getWinWidth(),getWinHeight());
#endif
  glViewport(0, 0, getWinWidth(),getWinHeight());
  //SHOULD SetWindowsSizeHint()...
  
/*
 glMatrixMode (GL_PROJECTION); // set up Frustum.
  glLoadIdentity();
  
  const G4Vector3D scaleFactor = fVP.GetScaleFactor();
  glScaled(scaleFactor.x(),scaleFactor.y(),scaleFactor.z());
  if (fVP.GetFieldHalfAngle() == 0.) {
    glOrtho (left, right, bottom, top, pnear, pfar);
  }
  else {
    glFrustum (left, right, bottom, top, pnear, pfar);
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
  gluLookAt (pCamera.x(),  pCamera.y(),  pCamera.z(),       // Viewpoint.
             gltarget.x(), gltarget.y(), gltarget.z(),      // Target point.
             upVector.x(), upVector.y(), upVector.z());     // Up vector.
  
  // Light position is "true" light direction, so must come after gluLookAt.
  glLightfv (GL_LIGHT0, GL_POSITION, lightPosition);
  
  // OpenGL no longer seems to reconstruct clipped edges, so, when the
  // BooleanProcessor is up to it, abandon this and use generic
  // clipping in G4OpenGLSceneHandler::CreateSectionPolyhedron.  Also,
  // force kernel visit on change of clipping plane in
  // G4OpenGLStoredViewer::CompareForKernelVisit.
  //if (fVP.IsSection () ) {  // pair of back to back clip planes.
  if (false) {  // pair of back to back clip planes.
    const G4Plane3D& sp = fVP.GetSectionPlane ();
    double sArray[4];
    sArray[0] = sp.a();
    sArray[1] = sp.b();
    sArray[2] = sp.c();
    sArray[3] = sp.d() + radius * 1.e-05;
    glClipPlane (CLIP_PLANE0, sArray);
    enable (Wt::WGLWidget::CLIP_PLANE0);
    sArray[0] = -sp.a();
    sArray[1] = -sp.b();
    sArray[2] = -sp.c();
    sArray[3] = -sp.d() + radius * 1.e-05;
    glClipPlane (CLIP_PLANE1, sArray);
    enable (Wt::WGLWidget::CLIP_PLANE1);
  } else {
    disable (Wt::WGLWidget::CLIP_PLANE0);
    disable (Wt::WGLWidget::CLIP_PLANE1);
  }

*/
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
/*    glClipPlane (CLIP_PLANE2, a); 
    enable (Wt::WGLWidget::CLIP_PLANE2); */
    if (nPlanes > 1) {
      a[0] = cutaways[1].a();
      a[1] = cutaways[1].b();
      a[2] = cutaways[1].c();
      a[3] = cutaways[1].d();
/*      glClipPlane (CLIP_PLANE3, a);
      enable (Wt::WGLWidget::CLIP_PLANE3); */
    }
    if (nPlanes > 2) {
      a[0] = cutaways[2].a();
      a[1] = cutaways[2].b();
      a[2] = cutaways[2].c();
      a[3] = cutaways[2].d();
/*      glClipPlane (CLIP_PLANE4, a);
      enable (Wt::WGLWidget::CLIP_PLANE4); */
    }
  } else {
/*    disable (Wt::WGLWidget::CLIP_PLANE2);
    disable (Wt::WGLWidget::CLIP_PLANE3);
    disable (Wt::WGLWidget::CLIP_PLANE4); */
  }
  
  // Background.
  background = fVP.GetBackgroundColour ();
  
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

void G4OpenGLImmediateWtViewer::wtDrawArrays(GLenum mode,int first,int nPoints, std::vector<double> a_vertices){
  std::vector<double> data2;
  for (int a=0; a< nPoints*3; a+=3) {
    data2.push_back(a_vertices[a]);
    data2.push_back(a_vertices[a+1]);
    data2.push_back(a_vertices[a+2]);
    data2.push_back(0);
    data2.push_back(0);
    data2.push_back(1);
  }

  //----------------------------
  // Fill VBO buffer
  //----------------------------
  
  // Create a new Buffer
  Buffer objBuffer_2 = glCreateBuffer(); //glGenBuffers(1,&objBuffer_2)

  VBO_Buffer.push_back(objBuffer_2);

  // Bind this buffer
  glBindBuffer(GL_ARRAY_BUFFER, objBuffer_2);

  // Load data into VBO
  glBufferDatafv(GL_ARRAY_BUFFER, data2.begin(), data2.end(), GL_STATIC_DRAW);
  
  
  //----------------------------
  // Draw VBO
  //----------------------------
  glBindBuffer(GL_ARRAY_BUFFER, objBuffer_2);

  // Configure the vertex attributes:
  vertexAttribPointer(vertexPositionAttribute_,
                      3,     // size: Every vertex has an X, Y anc Z component
                      GL_FLOAT, // type: They are floats
                      false, // normalized: Please, do NOT normalize the vertices
                      2*3*4, // stride: The first byte of the next vertex is located this
                      //         amount of bytes further. The format of the VBO is
                      //         vx, vy, vz, nx, ny, nz and every element is a
                      //         Float32, hence 4 bytes large
                      0);    // offset: The byte position of the first vertex in the buffer
  //         is 0.

  glDrawArrays(mode, first, data2.size()/6);
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

void G4OpenGLImmediateWtViewer::drawScene()
{
  // Clear color an depth buffers
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  // Configure the shader: set the uniforms
  // Uniforms are 'configurable constants' in a shader: they are
  // identical for every point that has to be drawn.
  // Set the camera transformation to the value of a client-side JS matrix
  glUniformMatrix4(cMatrixUniform_, jsMatrix_);
  // Often, a model matrix is used to move the model around. We're happy
  // with the location of the model, so we leave it as the unit matrix
  Wt::WMatrix4x4 modelMatrix;
  glUniformMatrix4(mvMatrixUniform_, modelMatrix);
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
  glUniformMatrix4(nMatrixUniform_, (jsMatrix_ * modelMatrix).inverted().transposed());
  
  // Configure the shaders: set the attributes.
  // Attributes are 'variables' within a shader: they vary for every point
  // that has to be drawn. All are stored in one VBO.
  
  // Create a Vertex Buffer Object (VBO) and load all polygon's data
  // (points, normals) into it. In this case we use one VBO that contains
  // all data (6 per point: vx, vy, vz, nx, ny, nz); alternatively you
  // can use multiple VBO's (e.g. one VBO for normals, one for points,
  // one for texture coordinates).
  // Note that if you use indexed buffers, you cannot have indexes
  // larger than 65K, due to the limitations of WebGL.
  
  
  
  /// LG : test
  std::vector<double> vertices;
  vertices.push_back(0);
  vertices.push_back(0);
  vertices.push_back(0);
  vertices.push_back(10.);
  vertices.push_back(0);
  vertices.push_back(0);
  vertices.push_back(0);
  vertices.push_back(0);
  vertices.push_back(0);
  vertices.push_back(0);
  vertices.push_back(10.);
  vertices.push_back(0);
  vertices.push_back(0);
  vertices.push_back(0);
  vertices.push_back(0);
  vertices.push_back(0);
  vertices.push_back(0);
  vertices.push_back(10.);
  wtDrawArrays(GL_LINES,0,6,vertices);
  
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLWtViewer drawScene Call ComputeView\n");
#endif
}


#endif
