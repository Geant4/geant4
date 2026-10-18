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
// Class G4VulkanImmediateQtViewer : a class derived from
//   G4VulkanQtViewer and G4VulkanImmediateViewer.

#ifndef G4VULKANIMMEDIATEQTVIEWER_HH
#define G4VULKANIMMEDIATEQTVIEWER_HH

#include "G4VulkanImmediateViewer.hh"
#include "G4VulkanQtViewer.hh"
#include "globals.hh"

class G4VulkanImmediateSceneHandler;

class G4VulkanImmediateQtViewer : public G4VulkanQtViewer,
                                  public G4VulkanImmediateViewer,
                                  public G4QVulkanWidgetType
{
  public:

    G4VulkanImmediateQtViewer(G4VulkanImmediateSceneHandler& scene, const G4String& name = "");
    ~G4VulkanImmediateQtViewer();
    void Initialise();
#if QT_VERSION < 0x060000
    void initializeVK();
#endif
    void DrawView();
    void resizeVK(int width, int height);
    void paintVK();
    void updateQWidget();
    void ShowView();
#if QT_VERSION < 0x060000
#else
    // G.Barrand: macOS: to avoid a crash at startup at first gl call.
    virtual void ClearView()
    {
      if (!G4QVulkanWidgetType::isValid()) return;
      G4VulkanViewer::ClearView();
    }
    virtual void SetView()
    {
      if (!G4QVulkanWidgetType::isValid()) return;
      G4VulkanViewer::SetView();
    }
#endif

  protected:

#if QT_VERSION < 0x060000
    void showEvent(QShowEvent* event);
#endif
    void wheelEvent(QWheelEvent* event);
    void mousePressEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
    void mouseDoubleClickEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    void contextMenuEvent(QContextMenuEvent* e);
    void keyPressEvent(QKeyEvent* event);
    void keyReleaseEvent(QKeyEvent* event);
#if QT_VERSION < 0x060000
    void paintEvent(QPaintEvent* event);
#endif

  private:

    void ComputeView();
};

#endif
