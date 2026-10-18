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
// Ported to Vulkan from the OpenGL work by Andrew Walkden (10th February 1997).
// Class G4VulkanImmediateXmViewer : a class derived from G4VulkanXmViewer
//                                     and G4VulkanImmediateViewer.

#include "G4VulkanImmediateXmViewer.hh"

#include "G4VulkanImmediateSceneHandler.hh"
#include "G4Threading.hh"
#include "G4ios.hh"

G4VulkanImmediateXmViewer::G4VulkanImmediateXmViewer(G4VulkanImmediateSceneHandler& sceneHandler,
                                                     const G4String& name)
  : G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name),
    G4VulkanViewer(sceneHandler),
    G4VulkanXmViewer(sceneHandler),
    G4VulkanImmediateViewer(sceneHandler)
{
  if (fViewId < 0) return;  // In case error in base class instantiation.

  // ensure a suitable window was found
  if (!vi_immediate)
  {
    G4cerr << "G4VulkanImmediateXmViewer::G4VulkanImmediateXmViewer -"
              " G4VulkanXmViewer couldn't get a visual."
           << G4endl;
    fViewId = -1;  // This flags an error.
    return;
  }
}

G4VulkanImmediateXmViewer::~G4VulkanImmediateXmViewer() {}

void G4VulkanImmediateXmViewer::Initialise()
{
  CreateVKXContext(vi_immediate);
  CreateMainWindow();
  CreateFontLists();

  InitializeVKView();

  // If a double buffer context has been forced upon us, ignore the
  // back buffer for this VulkanImmediate view.
  vkDrawBuffer(VK_FRONT);
}

void G4VulkanImmediateXmViewer::DrawView()
{
  G4ViewParameters::DrawingStyle style = GetViewParameters().GetDrawingStyle();

#ifdef G4DEBUG_VIS_VULKAN
  printf("G4VulkanImmediateXmViewer::DrawView : \n");
#endif

  if (style != G4ViewParameters::hlr && haloing_enabled)
  {
    HaloingFirstPass();
    NeedKernelVisit();
#ifdef G4DEBUG_VIS_VULKAN
    printf("G4VulkanImmediateXmViewer::DrawView : change param\n");
#endif
    ProcessView();
    vkFlush();

    HaloingSecondPass();
  }

#ifdef G4DEBUG_VIS_VULKAN
  printf("G4VulkanImmediateXmViewer::DrawView : need Kernel/Process/Finish\n");
#endif
  NeedKernelVisit();  // Always need to visit G4 kernel.
  ProcessView();
  FinishView();
}

void G4VulkanImmediateXmViewer::FinishView()
{
  vkFlush();
  PresentVulkanFramebuffer();
}
