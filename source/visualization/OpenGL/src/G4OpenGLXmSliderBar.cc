// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmSliderBar.cc,v 1.1 1999-01-07 16:15:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//Slider bar class. Inherits from G4OpenGLXmVWidgetComponent

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmVWidgetComponent.hh"
#include "G4OpenGLXmVWidgetContainer.hh"
#include "G4OpenGLXmSliderBar.hh"
#include <X11/Intrinsic.h>
#include <Xm/Scale.h>
#include "globals.hh"

G4OpenGLXmSliderBar::G4OpenGLXmSliderBar (char* n,
					  XtCallbackRec* c,
					  G4bool s,
					  short dp,
					  G4double v,
					  G4double max,
					  G4double min,
					  unsigned char o,
					  unsigned char d) 
{
  name = n;
  callback = c;
  show = s;
  decimal_places = dp;
  initial_value = int(v * pow(10.0, (G4double)dp));
  max_value = int(max * pow(10.0, (G4double)dp));
  min_value = int(min * pow(10.0, (G4double)dp));
  orientation = o;
  direction = d;
}

G4OpenGLXmSliderBar::~G4OpenGLXmSliderBar ()
{}

char* G4OpenGLXmSliderBar::GetName () 
{
  return name;
}

G4bool G4OpenGLXmSliderBar::GetShow () 
{
  return show;
}

short G4OpenGLXmSliderBar::GetDecimalPlaces () 
{
  return decimal_places;
}

G4double G4OpenGLXmSliderBar::GetInitialValue () 
{
  return (G4double)initial_value / pow(10.0, (G4double)GetDecimalPlaces());
}

G4double G4OpenGLXmSliderBar::GetMaxValue () 
{
  return (G4double)max_value / pow(10.0, (G4double)GetDecimalPlaces());
}

G4double G4OpenGLXmSliderBar::GetMinValue () 
{
  return (G4double)min_value / pow(10.0, (G4double)GetDecimalPlaces());
}

unsigned char G4OpenGLXmSliderBar::GetOrientation () 
{
  return orientation;
}

unsigned char G4OpenGLXmSliderBar::GetDirection () 
{
  return direction;
}

void G4OpenGLXmSliderBar::SetName (char* n) 
{
  name = n;
  XmString sliderbar_string = XmStringCreateLocalized (name);
  XtVaSetValues (sliderbar,
		 XmNlabelString, sliderbar_string,
		 NULL);
 XmStringFree (sliderbar_string);
}

void G4OpenGLXmSliderBar::SetShow (G4bool s) 
{
  show = s;
  XtVaSetValues (sliderbar,
		 XmNshowValue, show,
		 NULL);
  
}

void G4OpenGLXmSliderBar::SetDecimalPlaces (short dp) 
{
  decimal_places = dp;
  XtVaSetValues (sliderbar,
		 XmNdecimalPoints, decimal_places,
		 NULL);
  
}

void G4OpenGLXmSliderBar::SetInitialValue (G4double v) 
{
  initial_value = int(v * pow(10.0, (G4double)GetDecimalPlaces()));
  XtVaSetValues (sliderbar,
		 XmNvalue, initial_value,
		 NULL);
  
}

void G4OpenGLXmSliderBar::SetMaxValue (G4double v) 
{
  max_value = int(v * pow(10.0, (G4double)GetDecimalPlaces()));
  XtVaSetValues (sliderbar,
		 XmNmaximum, max_value,
		 NULL);
  
}

void G4OpenGLXmSliderBar::SetMinValue (G4double v) 
{
  min_value = int(v * pow(10.0, (G4double)GetDecimalPlaces()));
  XtVaSetValues (sliderbar,
		 XmNminimum, min_value,
		 NULL);
  
}

void G4OpenGLXmSliderBar::SetOrientation (unsigned char o) 
{
  orientation = o;
  XtVaSetValues (sliderbar,
		 XmNorientation, orientation,
		 NULL);
  
}

void G4OpenGLXmSliderBar::SetDirection (unsigned char d) 
{
  direction = d;
  XtVaSetValues (sliderbar,
		 XmNprocessingDirection, direction,
		 NULL);
  
}

void G4OpenGLXmSliderBar::AddYourselfTo (G4OpenGLXmVWidgetContainer* container)
{

  pView = container->GetView ();
  ProcesspView ();

  parent = container->GetPointerToWidget ();
  XmString name_string = XmStringCreateLocalized (name);
  sliderbar = XtVaCreateManagedWidget (name,
				       xmScaleWidgetClass,
				       *parent,
				       
				       XmNtitleString, name_string,
				       XmNmaximum, max_value,
				       XmNminimum, min_value,
				       XmNvalue, initial_value,
				       XmNshowValue, show,
				       XmNdecimalPoints, decimal_places,
				       XmNorientation, orientation,
				       XmNprocessingDirection, direction,
  
				       XtNvisual, visual,
				       XtNdepth, depth,
				       XtNcolormap, cmap,
				       XtNborderColor, borcol,
				       XtNbackground, bgnd,

				       NULL);
				       
  XtAddCallbacks (sliderbar,
		  XmNvalueChangedCallback,
		  callback);

  XtAddCallbacks (sliderbar,
		  XmNdragCallback,
		  callback);
  XmStringFree (name_string);
}

Widget* G4OpenGLXmSliderBar::GetPointerToParent ()
{
  return parent;
}

Widget* G4OpenGLXmSliderBar::GetPointerToWidget () 
{
  return &sliderbar;
}

#endif
