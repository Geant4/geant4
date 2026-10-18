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
//
//
// Alex Vergara Gil  17th October 2026
// Ported to Vulkan from the Geant4 OpenGL visualization driver.
// Class G4VulkanImmediateQtViewer : a class derived from G4VulkanQtViewer and
//                                G4VulkanImmediateViewer.

#include "G4VulkanImmediateQtViewer.hh"

#include "G4VulkanImmediateSceneHandler.hh"
#include "G4Threading.hh"
#include "G4ios.hh"

#include <qapplication.h>
#include <qtabwidget.h>
#if 0x060000 <= QT_VERSION
#  include "G4Qt.hh"
#endif

G4VulkanImmediateQtViewer::G4VulkanImmediateQtViewer(G4VulkanImmediateSceneHandler& sceneHandler,
                                                     const G4String& name)
  : G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name),
    G4VulkanViewer(sceneHandler),
    G4VulkanQtViewer(sceneHandler),
    G4VulkanImmediateViewer(sceneHandler)
{
#if QT_VERSION < 0x060000
  fQVKWidgetInitialiseCompleted = false;

  setFocusPolicy(Qt::StrongFocus);  // enable keybord events
  fHasToRepaint = false;
  fPaintEventLock = false;
  fUpdateVKLock = false;

  if (fViewId < 0) return;  // In case error in base class instantiation.
#else
  setFocusPolicy(Qt::StrongFocus);  // enable keybord events
#endif
}

G4VulkanImmediateQtViewer::~G4VulkanImmediateQtViewer() {}

void G4VulkanImmediateQtViewer::Initialise()
{
#if QT_VERSION < 0x060000

  fQVKWidgetInitialiseCompleted = false;
  CreateMainWindow(this, QString(GetName()));

  makeCurrent();
  vkDrawBuffer(VK_BACK);

  // set the good tab active
  if (G4QVulkanWidgetType::parentWidget())
  {
    auto* parentTab = dynamic_cast<QTabWidget*>(G4QVulkanWidgetType::parentWidget()->parent());
    if (parentTab)
    {
      parentTab->setCurrentIndex(parentTab->count() - 1);
    }
  }

  fQVKWidgetInitialiseCompleted = true;
#else
  CreateMainWindow(this, QString(GetName()));
  // Set jpg as default export format for Qt viewer
  setExportImageFormat("jpg");
#endif
}

#if QT_VERSION < 0x060000
void G4VulkanImmediateQtViewer::initializeVK()
{
  InitializeVKView();

  // If a double buffer context has been forced upon us, ignore the
  // back buffer for this VulkanImmediate view.
  //  vkDrawBuffer (VK_FRONT); // FIXME : Ne marche pas avec cette ligne, mais affiche le run
  //  correctement...

  if (fSceneHandler.GetScene() == 0)
  {
    fHasToRepaint = false;
  }
  else
  {
    fHasToRepaint = true;
  }

  // Set the component visible

  // and update it immediatly before wait for SessionStart() (batch mode)
  //  QCoreApplication::sendPostedEvents () ;

  // Set jpg as default export format for Qt viewer
  setExportImageFormat("jpg");
}
#endif

void G4VulkanImmediateQtViewer::DrawView()
{
#if QT_VERSION < 0x060000
#else
  if (IsGettingPickInfos())
  {
    paintVK();
    return;
  }
#endif
  if (G4Threading::IsMasterThread())
  {
    updateQWidget();
  }
}

void G4VulkanImmediateQtViewer::ComputeView()
{
#if QT_VERSION < 0x060000
  makeCurrent();
  // If a double buffer context has been forced upon us, ignore the
  // back buffer for this VulkanImmediate view.
  //  vkDrawBuffer (VK_FRONT);
#endif

  G4ViewParameters::DrawingStyle dstyle = GetViewParameters().GetDrawingStyle();

  if (dstyle != G4ViewParameters::hlr && haloing_enabled)
  {
    HaloingFirstPass();
    NeedKernelVisit();
    ProcessView();
    FinishView();
    HaloingSecondPass();
  }

  NeedKernelVisit();  // Always need to visit G4 kernel.
  ProcessView();

  if (isRecording())
  {
    savePPMToTemp();
  }

#if QT_VERSION < 0x060000
  fHasToRepaint = true;
#endif
}

/**
   - Lors du resize de la fenetre, on doit non pas redessiner le detecteur, mais aussi les
   evenements
*/
void G4VulkanImmediateQtViewer::resizeVK(int aWidth, int aHeight)
{
  if ((aWidth > 0) && (aHeight > 0))
  {
#if QT_VERSION < 0x060000
    ResizeWindow(aWidth, aHeight);
    fHasToRepaint = sizeHasChanged();
#else
    ResizeWindow(devicePixelRatio() * aWidth, devicePixelRatio() * aHeight);
#endif
  }
}

void G4VulkanImmediateQtViewer::paintVK()
{
#if QT_VERSION < 0x060000
  updateToolbarAndMouseContextMenu();
#else
  // G.Barrand: don't do any change in the GUI here, just "paint" this widget!
#endif

#if QT_VERSION < 0x060000
  if (fPaintEventLock)
  {
    //    return ;
  }
  if (!fQVKWidgetInitialiseCompleted)
  {
    fPaintEventLock = false;
    return;
  }
  if ((getWinWidth() == 0) && (getWinHeight() == 0))
  {
    return;
  }

  // DO NOT RESIZE IF SIZE HAS NOT CHANGE
  if (!fHasToRepaint)
  {
    // L. Garnier : Trap to get the size with mac OSX 10.6 and Qt 4.6(devel)
    // Tested on Qt4.5 on mac, 4.4 on windows, 4.5 on unbuntu
    int sw = 0;
    int sh = 0;
    if (!isMaximized() && !isFullScreen())
    {
      sw = normalGeometry().width();
      sh = normalGeometry().height();
    }
    else
    {
      sw = frameGeometry().width();
      sh = frameGeometry().height();
    }
    if ((getWinWidth() == (unsigned int)sw) && (getWinHeight() == (unsigned int)sh))
    {
      return;
    }
    else if ((sw == 0) && (sh == 0))
    {  // NOT A TOP LEVEL WIDGET
      if (((getWinWidth() == (unsigned int)width())) && (getWinHeight() == (unsigned int)height()))
      {
        return;
      }
    }
  }
#else
  if ((getWinWidth() == 0) && (getWinHeight() == 0)) return;  // G.Barrand: needed?
#endif

#if QT_VERSION < 0x060000
#else
  InitializeVKView();
  vkDrawBuffer(VK_BACK);
#endif

  SetView();

  ClearView();  // ok, put the background correct
  ComputeView();

#if QT_VERSION < 0x060000
  fHasToRepaint = false;  // could be set to false by ComputeView

  fPaintEventLock = false;
#endif
}

void G4VulkanImmediateQtViewer::mousePressEvent(QMouseEvent* event)
{
  G4MousePressEvent(event);
}

void G4VulkanImmediateQtViewer::keyPressEvent(QKeyEvent* event)
{
  G4keyPressEvent(event);
}

void G4VulkanImmediateQtViewer::keyReleaseEvent(QKeyEvent* event)
{
  G4keyReleaseEvent(event);
}

void G4VulkanImmediateQtViewer::wheelEvent(QWheelEvent* event)
{
  G4wheelEvent(event);
}

#if QT_VERSION < 0x060000
void G4VulkanImmediateQtViewer::showEvent(QShowEvent*)
{
  if (fQVKWidgetInitialiseCompleted)
  {
    fHasToRepaint = true;
  }
}
#endif

/**
 * This function was build in order to make a zoom on double clic event.
 * It was think to build a rubberband on the zoom area, but never work fine
 */
void G4VulkanImmediateQtViewer::mouseDoubleClickEvent(QMouseEvent*)
{
  G4MouseDoubleClickEvent();
}

void G4VulkanImmediateQtViewer::mouseReleaseEvent(QMouseEvent* event)
{
  G4MouseReleaseEvent(event);
}

void G4VulkanImmediateQtViewer::mouseMoveEvent(QMouseEvent* event)
{
  G4MouseMoveEvent(event);
}

void G4VulkanImmediateQtViewer::contextMenuEvent(QContextMenuEvent* e)
{
  G4manageContextMenuEvent(e);
}

#if QT_VERSION < 0x060000
void G4VulkanImmediateQtViewer::paintEvent(QPaintEvent*)
{
  if (!fQVKWidgetInitialiseCompleted)
  {
    return;
  }
  // Force a repaint next time if the FRAMEBUFFER is not READY
  fHasToRepaint = isFramebufferReady();
  if (fHasToRepaint)
  {
#  if (QT_VERSION < QT_VERSION_CHECK(6, 0, 0))
    updateVK();
#  else
    // Not sure this is correct....
    paintVK();
#  endif
  }
}
#endif

void G4VulkanImmediateQtViewer::updateQWidget()
{
#if QT_VERSION < 0x060000
  if (fUpdateVKLock)
  {
    return;
  }

  if (!isCurrentWidget())
  {
    return;
  }

  fUpdateVKLock = true;
  fHasToRepaint = true;
  repaint();
  updateViewerPropertiesTableWidget();
  updateSceneTreeWidget();
  fUpdateVKLock = false;
#else
  // if (!isCurrentWidget()) return; //G.Barrand: Qt must know if it has to activate paintVK() if
  // the widget is not visible. G.Barrand: don't do any change in the GUI here, just ask to "paint"
  // this widget!
  update();
#endif
}

void G4VulkanImmediateQtViewer::ShowView()
{
#if QT_VERSION < 0x060000
  fHasToRepaint = true;
  activateWindow();
#else
  activateWindow();
  ((QApplication*)G4Qt::getInstance())->processEvents();
#endif
}
