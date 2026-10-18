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
// Vulkan graphics system factory.

#include "G4VulkanStoredX.hh"

#include "G4VulkanStoredSceneHandler.hh"
#include "G4VulkanStoredXViewer.hh"
#include "G4VulkanViewer.hh"
#include "G4VulkanViewerMessenger.hh"
#include "G4VSceneHandler.hh"
#include "G4VisFeaturesOfVulkan.hh"

G4VulkanStoredX::G4VulkanStoredX()
  : G4VGraphicsSystem("VulkanStoredX", "VKSX", G4VisFeaturesOfVulkanSX(),
                      G4VGraphicsSystem::threeD)
{
  G4VulkanViewerMessenger::GetInstance();
}

G4VulkanStoredX::~G4VulkanStoredX() {}

G4VSceneHandler* G4VulkanStoredX::CreateSceneHandler(const G4String& name)
{
  G4VSceneHandler* pScene = new G4VulkanStoredSceneHandler(*this, name);
  return pScene;
}

G4VViewer* G4VulkanStoredX::CreateViewer(G4VSceneHandler& scene, const G4String& name)
{
  G4VViewer* pView = new G4VulkanStoredXViewer((G4VulkanStoredSceneHandler&)scene, name);
  if (pView)
  {
    if (pView->GetViewId() < 0)
    {
      delete pView;
      pView = 0;
      G4cerr << "G4VulkanStoredX::CreateViewer: error flagged by"
                " negative view id in G4VulkanStoredXViewer creation."
                "\n Destroying view and returning null pointer."
             << G4endl;
    }
  }
  else
  {
    G4cerr << "G4VulkanStoredX::CreateViewer: null pointer on"
              " new G4VulkanStoredXViewer."
           << G4endl;
  }
  return pView;
}
