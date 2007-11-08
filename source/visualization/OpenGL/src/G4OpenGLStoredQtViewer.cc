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
// $Id: G4OpenGLStoredQtViewer.cc,v 1.2 2007-11-08 17:00:51 lgarnier Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Class G4OpenGLStoredQtViewer : a class derived from G4OpenGLQtViewer and
//                                G4OpenGLStoredViewer.

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#include "G4OpenGLStoredQtViewer.hh"
#include "G4VisManager.hh"

#include "G4ios.hh"

#include <QtGui/QMouseEvent>
#include <QtGui/QContextMenuEvent>

G4OpenGLStoredQtViewer::G4OpenGLStoredQtViewer
(G4OpenGLStoredSceneHandler& sceneHandler,
 const G4String&  name):
 G4VViewer (sceneHandler, sceneHandler.IncrementViewCount (), name),
 G4OpenGLViewer (sceneHandler),
 G4OpenGLQtViewer (sceneHandler),
 G4OpenGLStoredViewer (sceneHandler),
 QGLWidget()                      // FIXME : gerer le pb du parent !
 {
   nbPaint =0;
   hasToRepaint =false;
  if (fViewId < 0) return;  // In case error in base class instantiation.
}

G4OpenGLStoredQtViewer::~G4OpenGLStoredQtViewer() {
   printf("GLWidget::~GLWidget \n");
   makeCurrent();
   // this is connect to the Dialog for deleting it properly
   // when close event.
   //   ((QDialog*)window())->reject();
   printf("GLWidget::~GLWidget END\n");
}

void G4OpenGLStoredQtViewer::Initialise() {
   printf("GLWidget::Initialise \n");
   readyToPaint = false;
   CreateGLQtContext ();
   printf("G4OpenGLStoredQtViewer::Initialise () 2\n");

  CreateMainWindow (this);
  printf("G4OpenGLStoredQtViewer::Initialise () 3\n");

  CreateFontLists ();  // FIXME Does nothing!
  
  printf("readyToPaint = true \n");
  readyToPaint = true;
  
  // First Draw
  SetView();
  printf("    ClearView\n");
  ClearView (); //ok, put the background correct
  ShowView();
  FinishView();
}

void G4OpenGLStoredQtViewer::initializeGL () {

   InitializeGLView ();

   printf("G4OpenGLStoredQtViewer::InitialiseGL () 1\n");

   // clear the buffers and window.
   ClearView ();
   //   printf("G4OpenGLStoredQtViewer::InitialiseGL () 2\n");
   FinishView ();
   
   glDepthFunc (GL_LEQUAL);
   glDepthMask (GL_TRUE);

   hasToRepaint =true;

   printf("G4OpenGLStoredQtViewer::InitialiseGL  -------------------------------------------------------------------------------------\n");
}


void G4OpenGLStoredQtViewer::DrawView () {

  printf("G4OpenGLStoredQtViewer::DrawView %d %d   VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV\n",WinSize_x, WinSize_y);
  printf("G4OpenGLStoredQtViewer::DrawView Dialog adress : %d\n",GLWindow);
   G4ViewParameters::DrawingStyle style = GetViewParameters().GetDrawingStyle();

   //Make sure current viewer is attached and clean...
   //Qt version needed
   glViewport (0, 0, WinSize_x, WinSize_y);

   //See if things have changed from last time and remake if necessary...
   // The fNeedKernelVisit flag might have been set by the user in
   // /vis/viewer/rebuild, but if not, make decision and set flag only
   // if necessary...
   if (!fNeedKernelVisit) 

   if (!fNeedKernelVisit) KernelVisitDecision ();
   
   G4bool kernelVisitWasNeeded = fNeedKernelVisit; // Keep (ProcessView resets).
   ProcessView ();
   

   if(style!=G4ViewParameters::hlr &&
      haloing_enabled) {
     printf("G4OpenGLStoredQtViewer::DrawView DANS LE IF\n");

     HaloingFirstPass ();
     DrawDisplayLists ();
     glFlush ();

     HaloingSecondPass ();

     DrawDisplayLists ();
     FinishView ();

   } else {
     printf("***************************  CASE 1 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n");
     
     // If kernel visit was needed, drawing and FinishView will already
     // have been done, so...
     if (!kernelVisitWasNeeded) {
       printf("***************************  CASE 2 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n");
       DrawDisplayLists ();
       FinishView ();
     } else {
       printf("***************************  CASE 3 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n");
       // However, union cutaways are implemented in DrawDisplayLists, so make
       // an extra pass...
       if (fVP.IsCutaway() &&
           fVP.GetCutawayMode() == G4ViewParameters::cutawayUnion) {
         ClearView();
         DrawDisplayLists ();
         FinishView ();
         printf("***************************  CASE 4 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n");
       } else { // ADD TO AVOID KernelVisit=1 and nothing to display
         DrawDisplayLists ();
         FinishView ();
       }
     }
   }

   printf("G4OpenGLStoredQtViewer::DrawView %d %d ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n",WinSize_x, WinSize_y);
   hasToRepaint =true;
}


//////////////////////////////////////////////////////////////////////////////
void G4OpenGLStoredQtViewer::FinishView (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  printf("G4OpenGLStoredQtViewer::FinishView VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV\n");

  glFlush ();
  swapBuffers ();
  printf("G4OpenGLStoredQtViewer::FinishView ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n");

}


/**
 - Lors du resize de la fenetre, on doit non pas redessiner le detecteur, mais aussi les evenements
 */
void G4OpenGLStoredQtViewer::resizeGL(
 int aWidth
,int aHeight)
{  
  glViewport(0, 0, aWidth, aHeight);
  glMatrixMode(GL_PROJECTION);
  glMatrixMode(GL_MODELVIEW);

   if (((WinSize_x != (G4int)aWidth)) || (WinSize_y != (G4int) aHeight)) {
     hasToRepaint =true;
   }
   WinSize_x = (G4int) aWidth;
   WinSize_y = (G4int) aHeight;

  printf("G4OpenGLStoredQtViewer::resizeGL ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %d %d=%d %d=%d\n",hasToRepaint,width(),aWidth,height(),aHeight);
}



void G4OpenGLStoredQtViewer::paintGL()
 {
   if (!readyToPaint) {
     printf("G4OpenGLStoredQtViewer::paintGL ============  Not ready %d\n",readyToPaint);
     readyToPaint= true;
     return;
   }
   // DO NOT RESIZE IF SIZE HAS NOT CHANGE :
   //    WHEN CLICK ON THE FRAME FOR EXAMPLE
   //    EXECEPT WHEN MOUSE MOVE EVENT
   if ( !hasToRepaint) {
     if (((WinSize_x == (G4int)width())) &&(WinSize_y == (G4int) height())) {
       printf("G4OpenGLStoredQtViewer::paintGL ============  Dont repaint\n");
       return;
     }
   }
   nbPaint++;
   printf("G4OpenGLStoredQtViewer::paintGL VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV %d ready %d\n",nbPaint,readyToPaint);
   WinSize_x = (G4int) width();
   WinSize_y = (G4int) height();
   
   glViewport (0, 0, width(), height());
   //   glLoadIdentity();
   

   SetView();

//   //  printf("before ClearView\n");
   printf("    ClearView\n");
   
   ClearView (); //ok, put the background correct
   DrawView();

   hasToRepaint =false;

   printf("G4OpenGLStoredQtViewer::paintGL ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %d ready %d\n",nbPaint,readyToPaint);
 }

void G4OpenGLStoredQtViewer::mousePressEvent(QMouseEvent *event)
{
  printf("G4OpenGLStoredQtViewer::mousePressEvent\n");
  G4MousePressEvent(event->pos());
}

void G4OpenGLStoredQtViewer::mouseMoveEvent(QMouseEvent *event)
{
  G4MouseMoveEvent(event->x(),event->y(),event->buttons());
  //  DrawView();
}


void G4OpenGLStoredQtViewer::contextMenuEvent(QContextMenuEvent *e)
{
  manageContextMenuEvent(e);
}

void G4OpenGLStoredQtViewer::updateQWidget() {
  hasToRepaint= true;
  updateGL();
  hasToRepaint= false;
}

#endif
