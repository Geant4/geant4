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
// Ported to Vulkan from the OpenGL work by Andrew Walkden (7th February 1997).
// Class G4VulkanStoredXViewer : a class derived from G4VulkanXViewer and
//                             G4VulkanStoredViewer.

#include "G4VulkanStoredXViewer.hh"

#include "G4VulkanStoredSceneHandler.hh"
#include "G4Threading.hh"
#include "G4ios.hh"

G4VulkanStoredXViewer::G4VulkanStoredXViewer(G4VulkanStoredSceneHandler& sceneHandler,
                                             const G4String& name)
  : G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name),
    G4VulkanViewer(sceneHandler),
    G4VulkanXViewer(sceneHandler),
    G4VulkanStoredViewer(sceneHandler)
{
  if (fViewId < 0) return;  // In case error in base class instantiation.

  if (!vi_stored)
  {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4VulkanStoredXViewer::G4VulkanStoredXViewer -"
              " G4VulkanXViewer couldn't get a visual."
           << G4endl;
    return;
  }
}

G4VulkanStoredXViewer::~G4VulkanStoredXViewer() {}

void G4VulkanStoredXViewer::Initialise()
{
#ifdef G4DEBUG_VIS_VULKAN
  printf("G4VulkanStoredXViewer::Initialise\n");
#endif
  CreateVKXContext(vi_stored);
  CreateMainWindow();
  CreateFontLists();

  InitializeVKView();

  vkDrawBuffer(VK_BACK);
}

void G4VulkanStoredXViewer::DrawView()
{
#ifdef G4DEBUG_VIS_VULKAN
  printf("G4VulkanStoredXViewer::DrawView\n");
#endif

  G4ViewParameters::DrawingStyle style = GetViewParameters().GetDrawingStyle();

  // See if things have changed from last time and remake if necessary...
  // The fNeedKernelVisit flag might have been set by the user in
  // /vis/viewer/rebuild, but if not, make decision and set flag only
  // if necessary...
  if (!fNeedKernelVisit) KernelVisitDecision();
  fLastVP = fVP;
  G4bool kernelVisitWasNeeded = fNeedKernelVisit;  // Keep (ProcessView resets).
  ProcessView();

  if (style != G4ViewParameters::hlr && haloing_enabled)
  {
    HaloingFirstPass();
    DrawDisplayLists();
#ifdef G4DEBUG_VIS_VULKAN
    printf("G4VulkanStoredXViewer::DrawView flush \n");
#endif
    vkFlush();

    HaloingSecondPass();

    DrawDisplayLists();
  }
  else
  {
    if (!kernelVisitWasNeeded)
    {
#ifdef G4DEBUG_VIS_VULKAN
      printf("G4VulkanStoredXViewer::DrawView NO need kernel visit\n");
#endif
      DrawDisplayLists();
    }
    else
    {
#ifdef G4DEBUG_VIS_VULKAN
      printf("G4VulkanStoredXViewer::DrawView NEED kernel visit\n");
#endif
      // However, union cutaways are implemented in DrawDisplayLists, so make
      // an extra pass...
      if (fVP.IsCutaway() && fVP.GetCutawayMode() == G4ViewParameters::cutawayUnion)
      {
        ClearView();
        DrawDisplayLists();
      }
      else
      {  // ADD TO AVOID KernelVisit=1 and nothing to display
        DrawDisplayLists();
      }
    }
  }

  FinishView();
}

void G4VulkanStoredXViewer::FinishView()
{
#ifdef G4DEBUG_VIS_VULKAN
  printf("G4VulkanStoredXViewer::FinishView\n");
#endif

  vkFlush();  // FIXME

  VKint renderMode;
  vkGetIntegerv(VK_RENDER_MODE, &renderMode);
  if (renderMode == VK_RENDER) PresentVulkanFramebuffer();
}
