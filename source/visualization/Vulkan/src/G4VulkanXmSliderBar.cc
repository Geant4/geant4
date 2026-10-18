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
// Slider bar class. Inherits from G4VulkanXmVWidgetComponent

#include "G4VulkanXmSliderBar.hh"

#include "G4VulkanXmVWidgetComponent.hh"
#include "G4VulkanXmVWidgetContainer.hh"
#include "globals.hh"

#include <X11/Intrinsic.h>
#include <Xm/Scale.h>

G4VulkanXmSliderBar::G4VulkanXmSliderBar(const char* n, XtCallbackRec* c, G4bool sh, short dp,
                                         G4double v, G4double max, G4double min, unsigned char o,
                                         unsigned char d)
  : sliderbar(0), parent(0)
{
  name = n;
  callback = c;
  show = sh;
  decimal_places = dp;
  initial_value = int(v * std::pow(10.0, (G4double)dp));
  max_value = int(max * std::pow(10.0, (G4double)dp));
  min_value = int(min * std::pow(10.0, (G4double)dp));
  orientation = o;
  direction = d;
}

G4VulkanXmSliderBar::~G4VulkanXmSliderBar() {}

const char* G4VulkanXmSliderBar::GetName()
{
  return name;
}

G4bool G4VulkanXmSliderBar::GetShow()
{
  return show;
}

short G4VulkanXmSliderBar::GetDecimalPlaces()
{
  return decimal_places;
}

G4double G4VulkanXmSliderBar::GetInitialValue()
{
  return (G4double)initial_value / std::pow(10.0, (G4double)GetDecimalPlaces());
}

G4double G4VulkanXmSliderBar::GetMaxValue()
{
  return (G4double)max_value / std::pow(10.0, (G4double)GetDecimalPlaces());
}

G4double G4VulkanXmSliderBar::GetMinValue()
{
  return (G4double)min_value / std::pow(10.0, (G4double)GetDecimalPlaces());
}

unsigned char G4VulkanXmSliderBar::GetOrientation()
{
  return orientation;
}

unsigned char G4VulkanXmSliderBar::GetDirection()
{
  return direction;
}

void G4VulkanXmSliderBar::SetName(const char* n)
{
  name = n;
  XmString sliderbar_string = XmStringCreateLocalized((char*)name);
  XtVaSetValues(sliderbar, XmNlabelString, sliderbar_string, NULL);
  XmStringFree(sliderbar_string);
}

void G4VulkanXmSliderBar::SetShow(G4bool sh)
{
  show = sh;
  XtVaSetValues(sliderbar, XmNshowValue, show, NULL);
}

void G4VulkanXmSliderBar::SetDecimalPlaces(short dp)
{
  decimal_places = dp;
  XtVaSetValues(sliderbar, XmNdecimalPoints, decimal_places, NULL);
}

void G4VulkanXmSliderBar::SetInitialValue(G4double v)
{
  initial_value = int(v * std::pow(10.0, (G4double)GetDecimalPlaces()));
  XtVaSetValues(sliderbar, XmNvalue, initial_value, NULL);
}

void G4VulkanXmSliderBar::SetMaxValue(G4double v)
{
  max_value = int(v * std::pow(10.0, (G4double)GetDecimalPlaces()));
  XtVaSetValues(sliderbar, XmNmaximum, max_value, NULL);
}

void G4VulkanXmSliderBar::SetMinValue(G4double v)
{
  min_value = int(v * std::pow(10.0, (G4double)GetDecimalPlaces()));
  XtVaSetValues(sliderbar, XmNminimum, min_value, NULL);
}

void G4VulkanXmSliderBar::SetOrientation(unsigned char o)
{
  orientation = o;
  XtVaSetValues(sliderbar, XmNorientation, orientation, NULL);
}

void G4VulkanXmSliderBar::SetDirection(unsigned char d)
{
  direction = d;
  XtVaSetValues(sliderbar, XmNprocessingDirection, direction, NULL);
}

void G4VulkanXmSliderBar::AddYourselfTo(G4VulkanXmVWidgetContainer* container)
{
  pView = container->GetView();
  ProcesspView();

  parent = container->GetPointerToWidget();
  XmString name_string = XmStringCreateLocalized((char*)name);
  sliderbar = XtVaCreateManagedWidget(name, xmScaleWidgetClass, *parent,

                                      XmNtitleString, name_string, XmNmaximum, max_value,
                                      XmNminimum, min_value, XmNvalue, initial_value, XmNshowValue,
                                      show, XmNdecimalPoints, decimal_places, XmNorientation,
                                      orientation, XmNprocessingDirection, direction,

                                      XtNvisual, visual, XtNdepth, depth, XtNcolormap, cmap,
                                      XtNborderColor, borcol, XtNbackground, bgnd,

                                      NULL);

  XtAddCallbacks(sliderbar, XmNvalueChangedCallback, callback);

  XtAddCallbacks(sliderbar, XmNdragCallback, callback);
  XmStringFree(name_string);
}

Widget* G4VulkanXmSliderBar::GetPointerToParent()
{
  return parent;
}

Widget* G4VulkanXmSliderBar::GetPointerToWidget()
{
  return &sliderbar;
}
