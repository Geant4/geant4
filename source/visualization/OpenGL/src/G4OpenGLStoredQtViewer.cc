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
// $Id: G4OpenGLStoredQtViewer.cc 103926 2017-05-03 13:43:27Z gcosmo $
//
//
// Class G4OpenGLStoredQtViewer : a class derived from G4OpenGLQtViewer and
//                                G4OpenGLStoredViewer.

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#include "G4OpenGLStoredQtViewer.hh"

#include "G4OpenGLStoredSceneHandler.hh"
#include "G4ios.hh"
#ifdef G4MULTITHREADED
#include "G4Threading.hh"
#endif

#include <qapplication.h>
#include <qtabwidget.h>

G4OpenGLStoredQtViewer::G4OpenGLStoredQtViewer
(G4OpenGLStoredSceneHandler& sceneHandler,
 const G4String&  name):
  G4VViewer (sceneHandler, sceneHandler.IncrementViewCount (), name),
  G4OpenGLViewer (sceneHandler),
  G4OpenGLQtViewer (sceneHandler),
  G4OpenGLStoredViewer (sceneHandler),             // FIXME : gerer le pb du parent !
  QGLWidget()
{
  fQGLWidgetInitialiseCompleted = false;

    // Indicates that the widget has no background, i.e. when the widget receives paint events, the background is not automatically repainted. Note: Unlike WA_OpaquePaintEvent, newly exposed areas are never filled with the background (e.g., after showing a window for the first time the user can see "through" it until the application processes the paint events). This flag is set or cleared by the widget's author.
  QGLWidget::setAttribute (Qt::WA_NoSystemBackground);

  setFocusPolicy(Qt::StrongFocus); // enable keybord events
  fHasToRepaint = false;
  fPaintEventLock = false;
  fUpdateGLLock = false;

  if (fViewId < 0) return;  // In case error in base class instantiation.
}

G4OpenGLStoredQtViewer::~G4OpenGLStoredQtViewer() {
  makeCurrent();
  // this is connect to the Dialog for deleting it properly
  // when close event.
  //   ((QDialog*)window())->reject();
}

void G4OpenGLStoredQtViewer::Initialise() {
  makeCurrent();

  fQGLWidgetInitialiseCompleted = false;
  CreateMainWindow (this,QString(GetName()));

  glDrawBuffer (GL_BACK);

   // set the good tab active
  if (QGLWidget::parentWidget()) {
    QTabWidget *parentTab = dynamic_cast<QTabWidget*> (QGLWidget::parentWidget()->parent()) ;
    if (parentTab) {
      parentTab->setCurrentIndex(parentTab->count()-1);
    }
  }
  
  fQGLWidgetInitialiseCompleted = true;
}

void G4OpenGLStoredQtViewer::initializeGL () {

  InitializeGLView ();

  if (fSceneHandler.GetScene() == 0) {
    fHasToRepaint =false;
  } else {
    fHasToRepaint =true;
  }
  
   // Set the component visible
  //   setVisible(true) ;

}

G4bool G4OpenGLStoredQtViewer::CompareForKernelVisit(G4ViewParameters& lastVP)
{
  // Identical to G4OpenGLStoredViewer::CompareForKernelVisit except
  // for checking of VisAttributesModifiers, because
  // G4OpenGLStoredQtViewer keeps track of its own touchable
  // modifiers (fTreeItemModels, etc.).
  if (
      (lastVP.GetDrawingStyle ()    != fVP.GetDrawingStyle ())    ||
      (lastVP.IsAuxEdgeVisible ()   != fVP.IsAuxEdgeVisible ())   ||
      (lastVP.IsCulling ()          != fVP.IsCulling ())          ||
      (lastVP.IsCullingInvisible () != fVP.IsCullingInvisible ()) ||
      (lastVP.IsDensityCulling ()   != fVP.IsDensityCulling ())   ||
      (lastVP.IsCullingCovered ()   != fVP.IsCullingCovered ())   ||
      (lastVP.IsSection ()          != fVP.IsSection ())          ||
      // Section (DCUT) implemented locally.  But still need to visit
      // kernel if status changes so that back plane culling can be
      // switched.
      (lastVP.IsCutaway ()          != fVP.IsCutaway ())          ||
      // Cutaways implemented locally.  But still need to visit kernel
      // if status changes so that back plane culling can be switched.
      (lastVP.IsExplode ()          != fVP.IsExplode ())          ||
      (lastVP.GetNoOfSides ()       != fVP.GetNoOfSides ())       ||
      (lastVP.GetDefaultVisAttributes()->GetColour() !=
       fVP.GetDefaultVisAttributes()->GetColour())                ||
      (lastVP.GetDefaultTextVisAttributes()->GetColour() !=
       fVP.GetDefaultTextVisAttributes()->GetColour())            ||
      (lastVP.GetBackgroundColour ()!= fVP.GetBackgroundColour ())||
      (lastVP.IsPicking ()          != fVP.IsPicking ()))
    return true;

  // Don't check VisAttributesModifiers if this comparison has been
  // initiated by a mouse interaction on the scene tree.
  if (fMouseOnSceneTree) {
    // Reset the flag.
    fMouseOnSceneTree = false;
  } else {
    // Not initiated by a mouse so compare for kernel visit.
    if (lastVP.GetVisAttributesModifiers() != fVP.GetVisAttributesModifiers()) {
      return true;
    }
  }

  if (lastVP.IsDensityCulling () &&
      (lastVP.GetVisibleDensity () != fVP.GetVisibleDensity ()))
    return true;

  /**************************************************************
   Section (DCUT) implemented locally.  No need to visit kernel if
   section plane itself changes.
   if (lastVP.IsSection () &&
   (lastVP.GetSectionPlane () != fVP.GetSectionPlane ()))
   return true;
   ***************************************************************/

  /**************************************************************
   Cutaways implemented locally.  No need to visit kernel if cutaway
   planes themselves change.
   if (lastVP.IsCutaway ()) {
   if (lastVP.GetCutawayPlanes ().size () !=
   fVP.GetCutawayPlanes ().size ()) return true;
   for (size_t i = 0; i < lastVP.GetCutawayPlanes().size(); ++i)
   if (lastVP.GetCutawayPlanes()[i] != fVP.GetCutawayPlanes()[i])
   return true;
   }
   ***************************************************************/

  if (lastVP.IsExplode () &&
      (lastVP.GetExplodeFactor () != fVP.GetExplodeFactor ()))
    return true;
  
  return false;
}

G4bool G4OpenGLStoredQtViewer::POSelected(size_t POListIndex)
{
  return isTouchableVisible(POListIndex);
}

G4bool G4OpenGLStoredQtViewer::TOSelected(size_t)
{
  return true;
}

void G4OpenGLStoredQtViewer::DrawView () {
  updateQWidget();
}

void G4OpenGLStoredQtViewer::ComputeView () {

  makeCurrent();
  G4ViewParameters::DrawingStyle dstyle = GetViewParameters().GetDrawingStyle();

  //Make sure current viewer is attached and clean...

  //See if things have changed from last time and remake if necessary...
  // The fNeedKernelVisit flag might have been set by the user in
  // /vis/viewer/rebuild, but if not, make decision and set flag only
  // if necessary...
  if (!fNeedKernelVisit) {
    KernelVisitDecision ();
  }
  G4bool kernelVisitWasNeeded = fNeedKernelVisit; // Keep (ProcessView resets).
  ProcessView ();
   
  if (kernelVisitWasNeeded) {
    displaySceneTreeComponent();
  }

  if(dstyle!=G4ViewParameters::hlr &&
     haloing_enabled) {

    HaloingFirstPass ();
    DrawDisplayLists ();
    glFlush ();

    HaloingSecondPass ();

    DrawDisplayLists ();
    FinishView ();

  } else {
     
    // If kernel visit was needed, drawing and FinishView will already
    // have been done, so...
    if (!kernelVisitWasNeeded) {
      DrawDisplayLists ();
      FinishView ();
    } else {
      // However, union cutaways are implemented in DrawDisplayLists, so make
      // an extra pass...
      if (fVP.IsCutaway() &&
          fVP.GetCutawayMode() == G4ViewParameters::cutawayUnion) {
        ClearView();
        DrawDisplayLists ();
        FinishView ();
      } else { // ADD TO AVOID KernelVisit=1 and nothing to display
        DrawDisplayLists ();
        FinishView ();
      }
    }
  }

  if (isRecording()) {
    savePPMToTemp();
  }

  fHasToRepaint = true;
}


/**
   - Lors du resize de la fenetre, on doit non pas redessiner le detecteur, mais aussi les evenements
*/
void G4OpenGLStoredQtViewer::resizeGL(
                                      int aWidth
                                      ,int aHeight)
{  
  // Set new size, it will be update when next Repaint()->SetView() called
  if ((aWidth > 0) && (aHeight > 0)) {
    ResizeWindow(aWidth,aHeight);
    fHasToRepaint = sizeHasChanged();
  }
}


// We have to get several case :
// - Only activate the windows (mouse click for example) -> Do not redraw
// - resize window -> redraw
// - try to avoid recompute everything if we do not rescale picture (side is the same)
 
void G4OpenGLStoredQtViewer::paintGL()
{
  updateToolbarAndMouseContextMenu();

  if (fPaintEventLock) {
//    return ;
  }
  fPaintEventLock = true;
  if ((getWinWidth() == 0) && (getWinHeight() == 0)) {
    return;
  }

  if (!fQGLWidgetInitialiseCompleted) {
    fPaintEventLock = false;
    return;
  }
  // DO NOT RESIZE IF SIZE HAS NOT CHANGE :
  //    WHEN CLICK ON THE FRAME FOR EXAMPLE
  //    EXECEPT WHEN MOUSE MOVE EVENT
  if ( !fHasToRepaint) {
    // L. Garnier : Trap to get the size with mac OSX 10.6 and Qt 4.6(devel)
    // Tested on Qt4.5 on mac, 4.4 on windows, 4.5 on unbuntu
    int sw = 0;
    int sh = 0;
    if (!isMaximized() && !isFullScreen()) {
      sw = normalGeometry().width();
      sh = normalGeometry().height();
    } else {
      sw = frameGeometry().width();
      sh = frameGeometry().height();
    }
    if ((getWinWidth() == (unsigned int)sw) &&(getWinHeight() == (unsigned int)sh)) {
      return;
    }
  }
  // Ensure that we really draw the BACK buffer
  glDrawBuffer (GL_BACK);

  SetView();
  
  ClearView (); //ok, put the background correct
  ComputeView();

  fHasToRepaint = false;

  fPaintEventLock = false;
}

void G4OpenGLStoredQtViewer::paintEvent(QPaintEvent *) {
  if (! fQGLWidgetInitialiseCompleted) {
    return;
  }
  // Force a repaint next time if the FRAMEBUFFER is not READY
  fHasToRepaint = isFramebufferReady();
  if ( fHasToRepaint) {
    // Will really update the widget by calling CGLFlushDrawable
    // The widget's rendering context will become the current context and initializeGL()
    // will be called if it hasn't already been called.
    // Copies the back buffer of a double-buffered context to the front buffer.
    updateGL();
  }
}

void G4OpenGLStoredQtViewer::mousePressEvent(QMouseEvent *event)
{
  G4MousePressEvent(event);
}

void G4OpenGLStoredQtViewer::keyPressEvent (QKeyEvent * event) 
{
  G4keyPressEvent(event);
}

void G4OpenGLStoredQtViewer::keyReleaseEvent (QKeyEvent * event)
{
  G4keyReleaseEvent(event);
}

void G4OpenGLStoredQtViewer::wheelEvent (QWheelEvent * event)
{
  G4wheelEvent(event);
}

void G4OpenGLStoredQtViewer::showEvent (QShowEvent *)
{
  if (fQGLWidgetInitialiseCompleted) {
    fHasToRepaint = true;
  }
}

/**
 * This function was build in order to make a zoom on double clic event.
 * It was think to build a rubberband on the zoom area, but never work fine
 */
void G4OpenGLStoredQtViewer::mouseDoubleClickEvent(QMouseEvent *)
{
  G4MouseDoubleClickEvent();
}

void G4OpenGLStoredQtViewer::mouseReleaseEvent(QMouseEvent *event)
{
  G4MouseReleaseEvent(event);
}

void G4OpenGLStoredQtViewer::mouseMoveEvent(QMouseEvent *event)
{
  G4MouseMoveEvent(event);
}


void G4OpenGLStoredQtViewer::contextMenuEvent(QContextMenuEvent *e)
{
  G4manageContextMenuEvent(e);
}

void G4OpenGLStoredQtViewer::updateQWidget() {
  if (fUpdateGLLock) {
    return;
  }
  
  if (! isCurrentWidget()){
    return;
  }

  fUpdateGLLock = true;
  fHasToRepaint= true;
  // Will really update the widget by calling CGLFlushDrawable
  // The widget's rendering context will become the current context and initializeGL()
  // will be called if it hasn't already been called.
  // Copies the back buffer of a double-buffered context to the front buffer.
  repaint(); // will read scene tree state
  // updateGL() // From J.Allison picking branch
  updateViewerPropertiesTableWidget();
  updateSceneTreeWidget();
  fUpdateGLLock = false;
}

void G4OpenGLStoredQtViewer::ShowView (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  //  glFlush ();  // Tentativley offered by JA 29/04/16.
  
  // Some X servers fail to draw all trajectories, particularly Mac
  // XQuartz.  Revisit this at a future date.  Meanwhile, issue an
  // extra...
  //  ClearView();  // Necessary?  JA 29/04/16 
  //  DrawView();   // Necessary?  JA 29/04/16
  activateWindow();
  //  glFlush(); // NO NEED and as drawView will already cause a flush
  // that could do a double flush
  
}


void G4OpenGLStoredQtViewer::DisplayTimePOColourModification (
G4Colour& c,
size_t poIndex) {
  c = getColorForPoIndex(poIndex);
}

#endif
