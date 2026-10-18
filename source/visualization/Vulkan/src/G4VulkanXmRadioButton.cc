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
// Alex Vergara Gil  17th October 2026
// Ported to Vulkan from the Geant4 OpenGL visualization driver.
// Radio button class. Inherits from G4VulkanXmVWidgetComponent

#include "G4VulkanXmRadioButton.hh"

#include "G4VulkanXmVWidgetComponent.hh"
#include "G4VulkanXmVWidgetContainer.hh"
#include "globals.hh"

#include <X11/Intrinsic.h>
#include <Xm/ToggleB.h>

G4VulkanXmRadioButton::G4VulkanXmRadioButton(const char* n, XtCallbackRec* c, G4bool d, G4int num)
  : button(0), parent(0)
{
  name = n;
  callback = c;
  default_button = d;
  number = num;
}

G4VulkanXmRadioButton::~G4VulkanXmRadioButton() {}

void G4VulkanXmRadioButton::SetName(const char* n)
{
  name = n;
  XmString button_string = XmStringCreateLocalized((char*)name);
  XtVaSetValues(button, XmNlabelString, button_string, NULL);
  XmStringFree(button_string);
}

const char* G4VulkanXmRadioButton::GetName()
{
  return name;
}

void G4VulkanXmRadioButton::AddYourselfTo(G4VulkanXmVWidgetContainer* container)
{
  pView = container->GetView();
  ProcesspView();
  parent = container->GetPointerToWidget();
  XmString button_string = XmStringCreateLocalized((char*)name);
  button = XtVaCreateManagedWidget(name, xmToggleButtonWidgetClass, *parent,

                                   XmNlabelString, button_string, XmNset, default_button,
                                   XmNuserData, number,

                                   XtNvisual, visual, XtNdepth, depth, XtNcolormap, cmap,
                                   XtNborderColor, borcol, XtNbackground, bgnd,

                                   NULL);

  XtAddCallbacks(button, XmNarmCallback, callback);

  XmStringFree(button_string);
}

Widget* G4VulkanXmRadioButton::GetPointerToParent()
{
  return parent;
}

Widget* G4VulkanXmRadioButton::GetPointerToWidget()
{
  return &button;
}
