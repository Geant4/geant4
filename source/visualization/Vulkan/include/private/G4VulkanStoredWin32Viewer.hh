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
// Class G4VulkanStoredWin32Viewer : a class derived from
//   G4VulkanWin32Viewer and G4VulkanStoredViewer.

#ifndef G4VULKANSTOREDWIN32VIEWER_HH
#define G4VULKANSTOREDWIN32VIEWER_HH

#include "G4VulkanStoredViewer.hh"
#include "G4VulkanWin32Viewer.hh"
#include "G4VViewer.hh"

class G4VulkanStoredSceneHandler;

class G4VulkanStoredWin32Viewer : public G4VulkanWin32Viewer, public G4VulkanStoredViewer
{
  public:

    G4VulkanStoredWin32Viewer(G4VulkanStoredSceneHandler& scene, const G4String& name = "");
    void Initialise();
    void DrawView();
    void FinishView();
};

#endif
