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
// Framed box container class

#include "G4VulkanXmFramedBox.hh"

#include "G4VulkanXmBox.hh"
#include "G4VulkanXmVWidgetComponent.hh"
#include "G4VulkanXmVWidgetShell.hh"

#include <Xm/Frame.h>
#include <Xm/RowColumn.h>

G4VulkanXmFramedBox::G4VulkanXmFramedBox(const char* n, G4bool r) : G4VulkanXmBox(n, r)
{
  frame = NULL;
}

G4VulkanXmFramedBox::~G4VulkanXmFramedBox() {}

void G4VulkanXmFramedBox::AddChild(G4VulkanXmVWidgetComponent* component)
{
  component->AddYourselfTo(this);
  Cardinal num_children;
  XtVaGetValues(box_row_col, XmNnumChildren, &num_children, NULL);
  //  G4cout << name << " now parents " << num_children << " children." << G4endl;
}

void G4VulkanXmFramedBox::AddYourselfTo(G4VulkanXmVWidgetShell* window)
{
  pView = window->GetView();
  ProcesspView();
  char framename[50];
  strcpy(framename, name);
  strcat(framename, "_frame");

  parent = window->GetPointerToWidget();
  frame = XtVaCreateManagedWidget(framename, xmFrameWidgetClass, *parent,

                                  XtNvisual, visual, XtNdepth, depth, XtNcolormap, cmap,
                                  XtNborderColor, borcol, XtNbackground, bgnd,

                                  NULL);

  box_row_col =
    XtVaCreateManagedWidget(name, xmRowColumnWidgetClass, frame,

                            XmNadjustMargin, True, XmNisHomogeneous, False, XmNlabelString,
                            (XmString)name, XmNradioAlwaysOne, radio, XmNradioBehavior, radio,

                            XtNvisual, visual, XtNdepth, depth, XtNcolormap, cmap, XtNborderColor,
                            borcol, XtNbackground, bgnd,

                            NULL);
}
