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
// $Id: G4OpenGLImmediateQtViewer.cc,v 1.6 2008-10-24 13:49:19 lgarnier Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Class G4OpenGLImmediateQtViewer : a class derived from G4OpenGLQtViewer and
//                                G4OpenGLImmediateViewer.

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#include "G4OpenGLImmediateQtViewer.hh"
#include "G4VisManager.hh"

#include "G4ios.hh"

G4OpenGLImmediateQtViewer::G4OpenGLImmediateQtViewer
(G4OpenGLImmediateSceneHandler& sceneHandler,
 const G4String&  name):
 G4VViewer (sceneHandler, sceneHandler.IncrementViewCount (), name),
 G4OpenGLViewer (sceneHandler),
 G4OpenGLQtViewer (sceneHandler),
 G4OpenGLImmediateViewer (sceneHandler),
 QGLWidget()                      // FIXME : gerer le pb du parent !
 {
   nbPaint =0;
  if (fViewId < 0) return;  // In case error in base class instantiation.
}

G4OpenGLImmediateQtViewer::~G4OpenGLImmediateQtViewer() {
     makeCurrent();
}

void G4OpenGLImmediateQtViewer::Initialise() {
#ifdef G4DEBUG
   printf("G4OpenGLImmediateQtViewer::Initialise \n");
#endif
   readyToPaint = false;
   CreateGLQtContext ();

   CreateMainWindow (this,QString(fName));

  CreateFontLists ();  // FIXME Does nothing!
  
  readyToPaint = true;
  
  // First Draw
  SetView();
#ifdef G4DEBUG
  printf("G4OpenGLImmediateQtViewer::Initialise    ClearView\n");
#endif
  ClearView (); //ok, put the background correct
  ShowView();
  FinishView();
}

void G4OpenGLImmediateQtViewer::initializeGL () {

   InitializeGLView ();

#ifdef G4DEBUG
   printf("G4OpenGLImmediateQtViewer::InitialiseGL ()\n");
#endif

   // If a double buffer context has been forced upon us, ignore the
   // back buffer for this OpenGLImmediate view.
   glDrawBuffer (GL_FRONT); // FIXME : Ne marche pas avec cette ligne, mais affiche le run correctement...
   // clear the buffers and window.
   ClearView ();
#ifdef G4DEBUG
   //   printf("G4OpenGLImmediateQtViewer::InitialiseGL () 2\n");
#endif
   FinishView ();
   


   glDepthFunc (GL_LEQUAL);
   glDepthMask (GL_TRUE);

#ifdef G4DEBUG
   printf("G4OpenGLImmediateQtViewer::InitialiseGL END\n");
#endif
}


void G4OpenGLImmediateQtViewer::DrawView () {

#ifdef G4DEBUG
  printf("G4OpenGLImmediateQtViewer::DrawView %d %d   VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV\n",WinSize_x, WinSize_y);
#endif
   // If a double buffer context has been forced upon us, ignore the
   // back buffer for this OpenGLImmediate view.
  glDrawBuffer (GL_FRONT);

   G4ViewParameters::DrawingStyle style = GetViewParameters().GetDrawingStyle();

   //Make sure current viewer is attached and clean...
   //Qt version needed
   //   glViewport (0, 0, WinSize_x, WinSize_y);

   if(style!=G4ViewParameters::hlr &&
      haloing_enabled) {

     HaloingFirstPass ();
     NeedKernelVisit ();
     ProcessView ();
     glFlush ();

     HaloingSecondPass ();

   }

   NeedKernelVisit ();  // Always need to visit G4 kernel.
   ProcessView ();
   FinishView ();
#ifdef G4DEBUG
  printf("G4OpenGLImmediateQtViewer::DrawView %d %d ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n",WinSize_x, WinSize_y);
#endif
  readyToPaint = false;
}


//////////////////////////////////////////////////////////////////////////////
void G4OpenGLImmediateQtViewer::FinishView (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
#ifdef G4DEBUG
  printf("G4OpenGLImmediateQtViewer::FinishView VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV\n");
#endif

   glFlush ();
#ifdef G4DEBUG
  printf("G4OpenGLImmediateQtViewer::FinishView ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n");
#endif

}


/**
 - Lors du resize de la fenetre, on doit non pas redessiner le detecteur, mais aussi les evenements
 */
void G4OpenGLImmediateQtViewer::resizeGL(
 int width
,int height)
{  
  setupViewport(width,height);
#ifdef G4DEBUG
  printf("G4OpenGLImmediateQtViewer::resizeGL ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n");
#endif
}


void G4OpenGLImmediateQtViewer::paintGL()
 {
   if (!readyToPaint) {
     readyToPaint= true;
     return;
   }
   // DO NOT RESIZE IF SIZE HAS NOT CHANGE
   if (((WinSize_x == (G4int)width())) &&(WinSize_y == (G4int) height())) {
     return;
   }
   nbPaint++;
#ifdef G4DEBUG
   printf("\n\nG4OpenGLImmediateQtViewer::paintGL VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV %d ready %d\n",nbPaint,readyToPaint);
#endif
   WinSize_x = (G4int) width();
   WinSize_y = (G4int) height();

   setupViewport(width(),height());

   SetView();
#ifdef G4DEBUG
//   //  printf("before ClearView\n");
#endif
#ifdef G4DEBUG
   printf("    ClearView\n");
#endif
   
   ClearView (); //ok, put the background correct
   DrawView();
   readyToPaint = true; // could be set to false by DrawView


#ifdef G4DEBUG
   printf("G4OpenGLImmediateQtViewer::paintGL ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %d ready %d\n\n\n",nbPaint,readyToPaint);
#endif
 }

void G4OpenGLImmediateQtViewer::updateQWidget() {
  updateGL();
}

#endif
