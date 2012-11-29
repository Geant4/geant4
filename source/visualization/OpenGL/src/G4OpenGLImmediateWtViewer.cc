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
// $Id$
//
//
// Class G4OpenGLImmediateWtViewer : a class derived from G4OpenGLWtViewer and
//                                G4OpenGLImmediateViewer.

#ifdef G4VIS_BUILD_OPENGLWT_DRIVER

#include "G4OpenGLImmediateWtViewer.hh"
#include "G4OpenGLImmediateSceneHandler.hh"

// Qt class
#include <qevent.h>

#include "G4ios.hh"

G4OpenGLImmediateWtViewer::G4OpenGLImmediateWtViewer
(//G4OpenGLImmediateSceneHandler& sceneHandler,
 /*Wt::WContainerWidget *aParent,*/ const G4String&  name)/*:
  G4VViewer (sceneHandler, sceneHandler.IncrementViewCount (), name),
  G4OpenGLViewer (sceneHandler),
  G4OpenGLImmediateViewer (sceneHandler)//,
//  Wt::WPaintedWidget(aParent)
*/
{

#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateWtViewer INIT\n");
#endif

  /**
     fQtViewer = new G4OpenGLImmediateQtViewer(sceneHandler, name);
     
     this->mouseWentDown().connect(this,&G4OpenGLImmediateWtViewer::WtMousePressEvent);
     this->keyPressed().connect(this,&G4OpenGLImmediateWtViewer::WtKeyPressEvent);
     this->mouseWheel().connect(this,&G4OpenGLImmediateWtViewer::WtWheelEvent);
     this->doubleClicked().connect(this,&G4OpenGLImmediateWtViewer::WtMouseDoubleClickEvent);
     this->mouseMoved().connect(this,&G4OpenGLImmediateWtViewer::WtMouseMoveEvent);
  */
  //  this->mouseWentUp().connect(this,&G4OpenGLImmediateWtViewer::WtMouseReleaseEvent);
}

G4OpenGLImmediateWtViewer::~G4OpenGLImmediateWtViewer() {
}



void  G4OpenGLImmediateWtViewer::DrawView() {
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateWtViewer DrawView\n");
#endif
  fQtViewer->updateQWidget();
  // FIXME et printEPS
}


void G4OpenGLImmediateWtViewer::WtMousePressEvent(Wt::WMouseEvent event)
{
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateWtViewer mousePress\n");
#endif
  // boutons et position
  fQtViewer->G4MousePressEvent(ConvertWtMouseEventToQt(event));
}

void G4OpenGLImmediateWtViewer::WtKeyPressEvent (Wt::WKeyEvent event) 
{
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateWtViewer keyPressEvent\n");
#endif
  fQtViewer->G4keyPressEvent(ConvertWtKeyEventToQt(event));
}

void G4OpenGLImmediateWtViewer::WtWheelEvent (Wt::WMouseEvent event) 
{
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateWtViewer wheelEvent\n");
#endif
  fQtViewer->G4wheelEvent(ConvertWtWheelEventToQt(event));
}

/**
   void G4OpenGLImmediateWtViewer::WtShowEvent (QShowEvent *) 
   {
   fHasToRepaint = true;
   }
*/

/**
 * This function was build in order to make a zoom on double clic event.
 * It was think to build a rubberband on the zoom area, but never work fine
 */
void G4OpenGLImmediateWtViewer::WtMouseDoubleClickEvent(Wt::WMouseEvent )
{
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateWtViewer mouseDoubleClickEvent\n");
#endif
  fQtViewer->G4MouseDoubleClickEvent();
}

/**
   void G4OpenGLImmediateWtViewer::WtMouseReleaseEvent(Wt::WMouseEvent )
   {
   G4MouseReleaseEvent();
   }
*/

void G4OpenGLImmediateWtViewer::WtMouseMoveEvent(Wt::WMouseEvent event)
{
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateWtViewer mouseMoveEvent\n");
#endif
  fQtViewer->G4MouseMoveEvent(ConvertWtMouseEventToQt(event));
}


/**
   void G4OpenGLImmediateWtViewer::contextMenuEvent(QContextMenuEvent *e)
   {
   G4manageContextMenuEvent(e);
   }
*/


void G4OpenGLImmediateWtViewer::paintEvent(Wt::WPaintDevice * painter) {
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateWtViewer paintEvent\n");
#endif
    fQtViewer->updateQWidget();
}





void G4OpenGLImmediateWtViewer::ShowView (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateWtViewer ShowView\n");
#endif
  fQtViewer->activateWindow();
}


QMouseEvent *  G4OpenGLImmediateWtViewer::ConvertWtMouseEventToQt(Wt::WMouseEvent event) {
  Qt::MouseButton bt;
  if (event.button() == (Wt::WMouseEvent::NoButton)) {
    bt = Qt::NoButton;
  } else if (event.button() == (Wt::WMouseEvent::LeftButton)) {
    bt = Qt::LeftButton;
  } else if (event.button() == (Wt::WMouseEvent::RightButton)) {
    bt = Qt::RightButton;
  } else if (event.button() == (Wt::WMouseEvent::MiddleButton)) {
    bt = Qt::MidButton;
  }

  Qt::KeyboardModifiers km;
  if (event.modifiers() == Wt::NoModifier) {
    km = Qt::NoModifier;
  } else if (event.modifiers() == Wt::ControlModifier) {
    km = Qt::ControlModifier;
  } else if (event.modifiers() == Wt::ShiftModifier) {
    km = Qt::ShiftModifier;
  } else if (event.modifiers() == Wt::AltModifier) {
    km = Qt::AltModifier;
  } else if (event.modifiers() == Wt::MetaModifier) {
    km = Qt::MetaModifier;
  }
  QPoint pt =  QPoint(event.widget().x,event.widget().y);
  QMouseEvent *e = new QMouseEvent(QEvent::None, pt, bt, bt, km);
  return e;
}


QWheelEvent *  G4OpenGLImmediateWtViewer::ConvertWtWheelEventToQt(Wt::WMouseEvent event) {
  Qt::MouseButton bt;
  if (event.button() == (Wt::WMouseEvent::NoButton)) {
    bt = Qt::NoButton;
  } else if (event.button() == (Wt::WMouseEvent::LeftButton)) {
    bt = Qt::LeftButton;
  } else if (event.button() == (Wt::WMouseEvent::RightButton)) {
    bt = Qt::RightButton;
  } else if (event.button() == (Wt::WMouseEvent::MiddleButton)) {
    bt = Qt::MidButton;
  }

  Qt::KeyboardModifiers km;
  if (event.modifiers() == Wt::NoModifier) {
    km = Qt::NoModifier;
  } else if (event.modifiers() == Wt::ControlModifier) {
    km = Qt::ControlModifier;
  } else if (event.modifiers() == Wt::ShiftModifier) {
    km = Qt::ShiftModifier;
  } else if (event.modifiers() == Wt::AltModifier) {
    km = Qt::AltModifier;
  } else if (event.modifiers() == Wt::MetaModifier) {
    km = Qt::MetaModifier;
  }
  QPoint pt =  QPoint(event.widget().x,event.widget().y);
  QWheelEvent *e = new QWheelEvent( pt, event.wheelDelta (), bt, km);
  return e;
}


QKeyEvent *  G4OpenGLImmediateWtViewer::ConvertWtKeyEventToQt(Wt::WKeyEvent event) {

  Qt::KeyboardModifiers km;
  if (event.modifiers() == Wt::NoModifier) {
    km = Qt::NoModifier;
  } else if (event.modifiers() == Wt::ControlModifier) {
    km = Qt::ControlModifier;
  } else if (event.modifiers() == Wt::ShiftModifier) {
    km = Qt::ShiftModifier;
  } else if (event.modifiers() == Wt::AltModifier) {
    km = Qt::AltModifier;
  } else if (event.modifiers() == Wt::MetaModifier) {
    km = Qt::MetaModifier;
  }
  QKeyEvent *e = new QKeyEvent( QEvent::None, event.charCode(), km );
  return e;
}
#endif
