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
// Class G4VulkanStoredQtViewer : a class derived from G4VulkanQtViewer and
//                                G4VulkanStoredViewer.

#include "G4VulkanStoredQtViewer.hh"

#include "G4VulkanStoredSceneHandler.hh"
#include "G4Threading.hh"
#include "G4UIQt.hh"
#include "G4ios.hh"
#if 0x060000 <= QT_VERSION
#  include "G4Qt.hh"
#endif

#include <qapplication.h>
#include <qtabwidget.h>

G4VulkanStoredQtViewer::G4VulkanStoredQtViewer(G4VulkanStoredSceneHandler& sceneHandler,
                                               const G4String& name)
  : G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name),
    G4VulkanViewer(sceneHandler),
    G4VulkanQtViewer(sceneHandler),
    G4VulkanStoredViewer(sceneHandler),  // FIXME : gerer le pb du parent !
    G4QVulkanWidgetType()
{
  if (fViewId < 0) return;  // In case error in base class instantiation.

#if QT_VERSION < 0x060000
  fQVKWidgetInitialiseCompleted = false;

  // Indicates that the widget has no background, i.e. when the widget receives paint events, the
  // background is not automatically repainted. Note: Unlike WA_OpaquePaintEvent, newly exposed
  // areas are never filled with the background (e.g., after showing a window for the first time the
  // user can see "through" it until the application processes the paint events). This flag is set
  // or cleared by the widget's author.
  G4QVulkanWidgetType::setAttribute(Qt::WA_NoSystemBackground);

  setFocusPolicy(Qt::StrongFocus);  // enable keybord events
  fHasToRepaint = false;
  fPaintEventLock = false;
  fUpdateVKLock = false;
#else
  // Indicates that the widget has no background, i.e. when the widget receives paint events, the
  // background is not automatically repainted. Note: Unlike WA_OpaquePaintEvent, newly exposed
  // areas are never filled with the background (e.g., after showing a window for the first time the
  // user can see "through" it until the application processes the paint events). This flag is set
  // or cleared by the widget's author.
  G4QVulkanWidgetType::setAttribute(Qt::WA_NoSystemBackground);

  setFocusPolicy(Qt::StrongFocus);  // enable keybord events
#endif
}

G4VulkanStoredQtViewer::~G4VulkanStoredQtViewer() {}

void G4VulkanStoredQtViewer::Initialise()
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
void G4VulkanStoredQtViewer::initializeVK()
{
  InitializeVKView();

  if (fSceneHandler.GetScene() == 0)
  {
    fHasToRepaint = false;
  }
  else
  {
    fHasToRepaint = true;
  }

  // Set the component visible
  //   setVisible(true) ;

  // Set jpg as default export format for Qt viewer
  setExportImageFormat("jpg");
}
#endif

G4bool G4VulkanStoredQtViewer::CompareForKernelVisit(G4ViewParameters& lastVP)
{
  // Identical to G4VulkanStoredViewer::CompareForKernelVisit except
  // for checking of VisAttributesModifiers, because
  // G4VulkanStoredQtViewer keeps track of its own touchable
  // modifiers (fTreeItemModels, etc.).
  if ((lastVP.GetDrawingStyle() != fVP.GetDrawingStyle())
      || (lastVP.GetNumberOfCloudPoints() != fVP.GetNumberOfCloudPoints())
      || (lastVP.IsAuxEdgeVisible() != fVP.IsAuxEdgeVisible())
      || (lastVP.IsCulling() != fVP.IsCulling())
      || (lastVP.IsCullingInvisible() != fVP.IsCullingInvisible())
      || (lastVP.IsDensityCulling() != fVP.IsDensityCulling())
      || (lastVP.IsCullingCovered() != fVP.IsCullingCovered())
      || (lastVP.GetCBDAlgorithmNumber() != fVP.GetCBDAlgorithmNumber()) ||
      // Note: Section and Cutaway can reveal back-facing faces. If
      // backface culling is implemented, the image can look strange because
      // the back-facing faces are not there. For the moment, we have disabled
      // (commented out) backface culling (it seems not to affect performance -
      // in fact, performance seems to improve), so there is no problem.
      (lastVP.IsSection() != fVP.IsSection()) ||
      // Section (DCUT) is NOT implemented locally so we need to visit the kernel.
      // (lastVP.IsCutaway ()          != fVP.IsCutaway ())          ||
      // Cutaways are implemented locally so we do not need to visit the kernel.
      (lastVP.IsExplode() != fVP.IsExplode()) || (lastVP.GetNoOfSides() != fVP.GetNoOfSides())
      || (lastVP.GetGlobalMarkerScale() != fVP.GetGlobalMarkerScale())
      || (lastVP.GetGlobalLineWidthScale() != fVP.GetGlobalLineWidthScale())
      || (lastVP.IsMarkerNotHidden() != fVP.IsMarkerNotHidden())
      || (lastVP.GetDefaultVisAttributes()->GetColour()
          != fVP.GetDefaultVisAttributes()->GetColour())
      || (lastVP.GetDefaultTextVisAttributes()->GetColour()
          != fVP.GetDefaultTextVisAttributes()->GetColour())
      || (lastVP.GetBackgroundColour() != fVP.GetBackgroundColour())
      || (lastVP.IsPicking() != fVP.IsPicking())
      || (lastVP.IsSpecialMeshRendering() != fVP.IsSpecialMeshRendering())
      || (lastVP.GetSpecialMeshRenderingOption() != fVP.GetSpecialMeshRenderingOption())
      || (lastVP.GetTransparencyByDepth() != fVP.GetTransparencyByDepth())
      || (lastVP.IsDotsSmooth() != fVP.IsDotsSmooth())
      || (lastVP.GetDotsSize() != fVP.GetDotsSize()))
    return true;

  // Don't check VisAttributesModifiers if this comparison has been
  // initiated by a mouse interaction on the scene tree.
  if (fMouseOnSceneTree)
  {
    // Reset the flag.
    fMouseOnSceneTree = false;
  }
  else
  {
    // Not initiated by a mouse so compare for kernel visit.
    if (lastVP.GetVisAttributesModifiers() != fVP.GetVisAttributesModifiers())
    {
      return true;
    }
  }

  if (lastVP.IsDensityCulling() && (lastVP.GetVisibleDensity() != fVP.GetVisibleDensity()))
    return true;

  //  /**************************************************************
  //   If section (DCUT) is implemented locally, comment this out.
  if (lastVP.IsSection() && (lastVP.GetSectionPlane() != fVP.GetSectionPlane())) return true;
  //   ***************************************************************/

  /**************************************************************
   If cutaways are implemented locally, comment this out.
   if (lastVP.IsCutaway ()) {
   if (vp.GetCutawayMode() != fVP.GetCutawayMode()) return true;
   if (lastVP.GetCutawayPlanes ().size () !=
   fVP.GetCutawayPlanes ().size ()) return true;
   for (size_t i = 0; i < lastVP.GetCutawayPlanes().size(); ++i)
   if (lastVP.GetCutawayPlanes()[i] != fVP.GetCutawayPlanes()[i])
   return true;
   }
   ***************************************************************/

  if (lastVP.GetCBDAlgorithmNumber() > 0)
  {
    if (lastVP.GetCBDParameters().size() != fVP.GetCBDParameters().size())
      return true;
    else if (lastVP.GetCBDParameters() != fVP.GetCBDParameters())
      return true;
  }

  if (lastVP.IsExplode() && (lastVP.GetExplodeFactor() != fVP.GetExplodeFactor())) return true;

  if (lastVP.IsSpecialMeshRendering()
      && (lastVP.GetSpecialMeshVolumes() != fVP.GetSpecialMeshVolumes()))
    return true;

  if (lastVP.GetTransparencyByDepth() > 0.
      && lastVP.GetTransparencyByDepthOption() != fVP.GetTransparencyByDepthOption())
    return true;

  return false;
}

G4bool G4VulkanStoredQtViewer::POSelected(size_t POListIndex)
{
  return isTouchableVisible((int)POListIndex);
}

G4bool G4VulkanStoredQtViewer::TOSelected(size_t)
{
  return true;
}

void G4VulkanStoredQtViewer::DrawView()
{
#if QT_VERSION < 0x060000
#else
  if (IsGettingPickInfos())
  {
    paintVK();
    return;
  }
#endif
#if (QT_VERSION < 0x060000) || !defined(G4MULTITHREADED)
  updateQWidget();
#else
  if (G4Threading::IsMasterThread())
  {
    updateQWidget();
  }
  else
  {
    update();  // G.Barrand: updateQWidget() induces a crash on run beamOn.
  }
#endif
}

void G4VulkanStoredQtViewer::ComputeView()
{
#if QT_VERSION < 0x060000
  makeCurrent();
#endif

  G4ViewParameters::DrawingStyle dstyle = GetViewParameters().GetDrawingStyle();

  // Make sure current viewer is attached and clean...

  // See if things have changed from last time and remake if necessary...
  //  The fNeedKernelVisit flag might have been set by the user in
  //  /vis/viewer/rebuild, but if not, make decision and set flag only
  //  if necessary...
  if (!fNeedKernelVisit)
  {
    KernelVisitDecision();
  }
  fLastVP = fVP;
  G4bool kernelVisitWasNeeded = fNeedKernelVisit;  // Keep (ProcessView resets).
  ProcessView();

  if (kernelVisitWasNeeded)
  {
    displaySceneTreeComponent();
  }

  if (dstyle != G4ViewParameters::hlr && haloing_enabled)
  {
    HaloingFirstPass();
    DrawDisplayLists();
    vkFlush();

    HaloingSecondPass();

    DrawDisplayLists();
    FinishView();
  }
  else
  {
    // If kernel visit was needed, drawing and FinishView will already
    // have been done, so...
    if (!kernelVisitWasNeeded)
    {
      DrawDisplayLists();
      FinishView();
    }
    else
    {
      // However, union cutaways are implemented in DrawDisplayLists, so make
      // an extra pass...
      if (fVP.IsCutaway() && fVP.GetCutawayMode() == G4ViewParameters::cutawayUnion)
      {
        ClearView();
        DrawDisplayLists();
        FinishView();
      }
      else
      {  // ADD TO AVOID KernelVisit=1 and nothing to display
        DrawDisplayLists();
        FinishView();
      }
    }
  }

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
void G4VulkanStoredQtViewer::resizeVK(int aWidth, int aHeight)
{
  // Set new size, it will be update when next Repaint()->SetView() called
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

// We have to get several case :
// - Only activate the windows (mouse click for example) -> Do not redraw
// - resize window -> redraw
// - try to avoid recompute everything if we do not rescale picture (side is the same)

void G4VulkanStoredQtViewer::paintVK()
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
  fPaintEventLock = true;
  if ((getWinWidth() == 0) && (getWinHeight() == 0))
  {
    return;
  }

  if (!fQVKWidgetInitialiseCompleted)
  {
    fPaintEventLock = false;
    return;
  }

  // DO NOT RESIZE IF SIZE HAS NOT CHANGE :
  //    WHEN CLICK ON THE FRAME FOR EXAMPLE
  //    EXECEPT WHEN MOUSE MOVE EVENT
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
  }
#else
  if ((getWinWidth() == 0) && (getWinHeight() == 0)) return;  // G.Barrand: needed?
#endif

#if QT_VERSION < 0x060000
#else
  InitializeVKView();
#endif

  // Ensure that we really draw the BACK buffer
  vkDrawBuffer(VK_BACK);

  SetView();

  ClearView();  // ok, put the background correct
  ComputeView();

#if QT_VERSION < 0x060000
  fHasToRepaint = false;

  fPaintEventLock = false;
#endif
}

#if QT_VERSION < 0x060000
void G4VulkanStoredQtViewer::paintEvent(QPaintEvent*)
{
  if (!fQVKWidgetInitialiseCompleted)
  {
    return;
  }
  // Force a repaint next time if the FRAMEBUFFER is not READY
  fHasToRepaint = isFramebufferReady();
  if (fHasToRepaint)
  {
    // Will really update the widget by calling CGLFlushDrawable
    // The widget's rendering context will become the current context and initializeVK()
    // will be called if it hasn't already been called.
    // Copies the back buffer of a double-buffered context to the front buffer.
#  if (QT_VERSION < QT_VERSION_CHECK(6, 0, 0))
    updateVK();
#  else
    // Not sure this is correct....
    paintVK();
#  endif
  }
}
#endif

void G4VulkanStoredQtViewer::mousePressEvent(QMouseEvent* event)
{
  G4MousePressEvent(event);
}

void G4VulkanStoredQtViewer::keyPressEvent(QKeyEvent* event)
{
  G4keyPressEvent(event);
}

void G4VulkanStoredQtViewer::keyReleaseEvent(QKeyEvent* event)
{
  G4keyReleaseEvent(event);
}

void G4VulkanStoredQtViewer::wheelEvent(QWheelEvent* event)
{
  G4wheelEvent(event);
}

#if QT_VERSION < 0x060000
void G4VulkanStoredQtViewer::showEvent(QShowEvent*)
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
void G4VulkanStoredQtViewer::mouseDoubleClickEvent(QMouseEvent*)
{
  G4MouseDoubleClickEvent();
}

void G4VulkanStoredQtViewer::mouseReleaseEvent(QMouseEvent* event)
{
  G4MouseReleaseEvent(event);
}

void G4VulkanStoredQtViewer::mouseMoveEvent(QMouseEvent* event)
{
  G4MouseMoveEvent(event);
}

void G4VulkanStoredQtViewer::contextMenuEvent(QContextMenuEvent* e)
{
  G4manageContextMenuEvent(e);
}

void G4VulkanStoredQtViewer::updateQWidget()
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
  // Will really update the widget by calling CGLFlushDrawable
  // The widget's rendering context will become the current context and initializeVK()
  // will be called if it hasn't already been called.
  // Copies the back buffer of a double-buffered context to the front buffer.
  repaint();  // will read scene tree state
  // updateVK() // From J.Allison picking branch
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

void G4VulkanStoredQtViewer::ShowView()
{
  activateWindow();
#if 0x060000 <= QT_VERSION
  ((QApplication*)G4Qt::getInstance())->processEvents();
#endif
}

void G4VulkanStoredQtViewer::DisplayTimePOColourModification(G4Colour& c, size_t poIndex)
{
  c = getColorForPoIndex((int)poIndex);
}
