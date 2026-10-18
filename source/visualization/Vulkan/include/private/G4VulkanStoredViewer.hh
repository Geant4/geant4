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
// Class G4VulkanStoredViewer : Encapsulates the `storedness' of
//                              an Vulkan viewer, for inheritance by
//                              derived (X, Xm...) classes.

#ifndef G4VULKANSTOREDVIEWER_HH
#define G4VULKANSTOREDVIEWER_HH

#include "G4VulkanViewer.hh"

class G4VulkanStoredSceneHandler;
class G4Colour;
class G4Text;
class G4Circle;

class G4VulkanStoredViewer : virtual public G4VulkanViewer
{
  public:

    G4VulkanStoredViewer(G4VulkanStoredSceneHandler& scene);
    virtual ~G4VulkanStoredViewer();

  protected:

    void KernelVisitDecision();
    virtual G4bool CompareForKernelVisit(G4ViewParameters&);

    void DrawDisplayLists();

    virtual void DisplayTimePOColourModification(G4Colour&, size_t /*currentPOListIndex*/) {}

    void AddPrimitiveForASingleFrame(const G4Text& text);
    void AddPrimitiveForASingleFrame(const G4Circle& circle);

    G4VulkanStoredSceneHandler& fG4VulkanStoredSceneHandler;

    G4ViewParameters fLastVP;  // Memory for making kernel visit decisions.

    // Two virtual functions to return sub-class selection.
    virtual G4bool POSelected(size_t) { return true; }
    virtual G4bool TOSelected(size_t) { return true; }

    G4bool fDepthTestEnable;
    G4Colour fOldDisplayListColor;
};

#endif
