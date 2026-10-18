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
// Class G4VulkanStoredWin32Viewer : a class derived from G4VulkanWin32Viewer and
//                             G4VulkanStoredViewer.

#include "G4VulkanStoredWin32Viewer.hh"

#include "G4VulkanStoredSceneHandler.hh"
#include "G4ios.hh"

G4VulkanStoredWin32Viewer::G4VulkanStoredWin32Viewer(G4VulkanStoredSceneHandler& sceneHandler,
                                                     const G4String& name)
  : G4VulkanViewer(sceneHandler),
    G4VulkanWin32Viewer(sceneHandler),
    G4VulkanStoredViewer(sceneHandler),
    G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name)
{
  if (fViewId < 0) return;  // In case error in base class instantiation.
}

void G4VulkanStoredWin32Viewer::Initialise()
{
  // Check that G4VulkanWin32Viewer got a double buffered colour visual

  CreateVKWin32Context();
  CreateMainWindow();
  CreateFontLists();

  // clear the buffers and window.
  ClearView();
  FinishView();

  vkDepthFunc(VK_LEQUAL);
  vkDepthMask(VK_TRUE);
}

void G4VulkanStoredWin32Viewer::DrawView()
{
  vkViewport(0, 0, getWinWidth(), getWinHeight());

  G4ViewParameters::DrawingStyle style = GetViewParameters().GetDrawingStyle();

  // See if things have changed from last time and remake if necessary...
  //  The fNeedKernelVisit flag might have been set by the user in
  //  /vis/viewer/rebuild, but if not, make decision and set flag only
  //  if necessary...
  if (!fNeedKernelVisit) KernelVisitDecision();
  fLastVP = fVP;
  G4bool kernelVisitWasNeeded = fNeedKernelVisit;  // Keep (ProcessView resets).
  ProcessView();

  if (style != G4ViewParameters::hlr && haloing_enabled)
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
}

//////////////////////////////////////////////////////////////////////////////
void G4VulkanStoredWin32Viewer::FinishView()
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if (!fHDC) return;

  vkFlush();
  VKint renderMode;
  vkGetIntegerv(VK_RENDER_MODE, &renderMode);
  if (renderMode == VK_RENDER) PresentVulkanFramebuffer();

  // Empty the Windows message queue :
  MSG event;
  while (::PeekMessage(&event, NULL, 0, 0, PM_REMOVE))
  {
    ::TranslateMessage(&event);
    ::DispatchMessage(&event);
  }
}
