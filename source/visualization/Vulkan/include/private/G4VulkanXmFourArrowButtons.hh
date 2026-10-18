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
// Four arrow buttons class. Inherits from G4VulkanXmVWidgetComponent

#ifndef G4VULKANXMFOURARROWBUTTONS_HH
#define G4VULKANXMFOURARROWBUTTONS_HH

#include "G4VulkanXmVWidgetComponent.hh"

class G4VulkanXmFourArrowButtons : public G4VulkanXmVWidgetComponent
{
  public:

    G4VulkanXmFourArrowButtons(XtCallbackRec** = NULL);  // array of 4 callbacks
                                                         // constructor
    virtual ~G4VulkanXmFourArrowButtons();  // destructor

    void SetName(char*);

    char* GetName();

    void AddYourselfTo(G4VulkanXmVWidgetContainer*);

    Widget* GetPointerToParent();
    Widget* GetPointerToWidget();

  private:

    XtCallbackRec** callback;
    Widget arrow_form;
    Widget arrow;
    Widget* parent;

  private:

    G4VulkanXmFourArrowButtons(const G4VulkanXmFourArrowButtons&);
    G4VulkanXmFourArrowButtons& operator=(const G4VulkanXmFourArrowButtons&);
};

#endif
