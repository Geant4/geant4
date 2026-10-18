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
// Push button class. Inherits from G4VulkanXmVWidgetComponent

#ifndef G4VULKANXMPUSHBUTTON_HH
#define G4VULKANXMPUSHBUTTON_HH

#include "G4VulkanXmVWidgetComponent.hh"

class G4VulkanXmPushButton : public G4VulkanXmVWidgetComponent
{
  public:

    G4VulkanXmPushButton(const char* = NULL,
                         XtCallbackRec* = NULL);  // constructor
    virtual ~G4VulkanXmPushButton();  // destructor

    void SetName(const char*);
    const char* GetName();

    void AddYourselfTo(G4VulkanXmVWidgetContainer*);

    Widget* GetPointerToParent();
    Widget* GetPointerToWidget();

  private:

    G4VulkanXmPushButton(const G4VulkanXmPushButton&);
    G4VulkanXmPushButton& operator=(const G4VulkanXmPushButton&);
    const char* name;
    XtCallbackRec* callback;
    Widget button;
    Widget* parent;
};

#endif
