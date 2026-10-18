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
// Text field class. Inherits from G4VulkanXmVWidgetComponent

#ifndef G4VULKANXMTEXTFIELD_HH
#define G4VULKANXMTEXTFIELD_HH

#include "G4VulkanXmVWidgetComponent.hh"
// #include "G4VulkanXmConvenienceRoutines.hh"

class G4VulkanXmTextField : public G4VulkanXmVWidgetComponent
{
  public:

    G4VulkanXmTextField(const char*, G4double*);  // constructor
    G4VulkanXmTextField(const char*, const char*);  // constructor
    virtual ~G4VulkanXmTextField();  // destructor

    void SetName(const char*);
    const char* GetName();

    void SetValue(G4double);
    void SetValue(const char*);
    const char* GetValue();

    void AddYourselfTo(G4VulkanXmVWidgetContainer*);

    Widget* GetPointerToParent();
    Widget* GetPointerToWidget();

  private:

    G4VulkanXmTextField(const G4VulkanXmTextField&);
    G4VulkanXmTextField& operator=(const G4VulkanXmTextField&);
    const char* name;
    void* value;
    G4bool text;
    char* initial;
    Widget text_label;
    Widget text_field;
    Widget* parent;
};

#endif
