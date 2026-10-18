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
// Class G4VulkanImmediateXViewer : a class derived from G4VulkanXViewer and
//                                G4VulkanImmediateViewer.

#include "G4VulkanImmediateXViewer.hh"

#include "G4VulkanImmediateSceneHandler.hh"
#include "G4Threading.hh"
#include "G4ios.hh"

G4VulkanImmediateXViewer::G4VulkanImmediateXViewer(G4VulkanImmediateSceneHandler& sceneHandler,
                                                   const G4String& name)
  : G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name),
    G4VulkanViewer(sceneHandler),
    G4VulkanXViewer(sceneHandler),
    G4VulkanImmediateViewer(sceneHandler)
{
  if (fViewId < 0) return;  // In case error in base class instantiation.

  // ensure a suitable window was found
  if (!vi_immediate)
  {
    G4cerr << "G4VulkanImmediateXViewer::G4VulkanImmediateXViewer -"
              " G4VulkanXViewer couldn't get a visual."
           << G4endl;
    fViewId = -1;  // This flags an error.
    return;
  }
}

G4VulkanImmediateXViewer::~G4VulkanImmediateXViewer() {}

void G4VulkanImmediateXViewer::Initialise()
{
  CreateVKXContext(vi_immediate);
  CreateMainWindow();
  CreateFontLists();

  InitializeVKView();

  // If a double buffer context has been forced upon us, ignore the
  // back buffer for this VulkanImmediate view.
  vkDrawBuffer(VK_FRONT);

  vkDepthFunc(VK_LEQUAL);
  vkDepthMask(VK_TRUE);
}

void G4VulkanImmediateXViewer::DrawView()
{
  G4ViewParameters::DrawingStyle style = GetViewParameters().GetDrawingStyle();

  if (style != G4ViewParameters::hlr && haloing_enabled)
  {
    HaloingFirstPass();
    NeedKernelVisit();
    ProcessView();
    vkFlush();

    HaloingSecondPass();
  }

  NeedKernelVisit();  // Always need to visit G4 kernel.
  ProcessView();
  FinishView();
}

void G4VulkanImmediateXViewer::FinishView()
{
  vkFlush();
  PresentVulkanFramebuffer();
}
