// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmSliderBar.hh,v 1.2 1999-01-09 16:23:01 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//Slider bar class. Inherits from G4OpenGLXmVWidgetComponent

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMSLIDERBAR_HH
#define G4OPENGLXMSLIDERBAR_HH

#include "G4OpenGLXmVWidgetComponent.hh"

class G4OpenGLXmSliderBar : public G4OpenGLXmVWidgetComponent
{

public:
  G4OpenGLXmSliderBar (char* = NULL,           // name of slider bar
		       XtCallbackRec* = NULL,  // callbacks for slider bar
		       G4bool = False,         // show current value if True
		       short = 0,              // decimal places for show value
		       G4double = 0.,          // initial value
		       G4double = 0.,          // max value
		       G4double = 0.,          // min value
		       unsigned char = XmHORIZONTAL,
		       unsigned char = XmMAX_ON_RIGHT); 
                                               //constructor
  ~G4OpenGLXmSliderBar ();                     //destructor

  void SetName (char*);
  void SetShow (G4bool);
  void SetDecimalPlaces (short);
  void SetInitialValue (G4double);
  void SetMaxValue (G4double);
  void SetMinValue (G4double);
  void SetOrientation (unsigned char);
  void SetDirection (unsigned char);

 
  char* GetName ();
  G4bool GetShow ();
  short GetDecimalPlaces ();
  G4double GetInitialValue ();
  G4double GetMaxValue ();
  G4double GetMinValue ();
  unsigned char GetOrientation ();
  unsigned char GetDirection ();

  void AddYourselfTo (G4OpenGLXmVWidgetContainer*);

  Widget* GetPointerToParent ();
  Widget* GetPointerToWidget ();

private:
  char* name;
  XtCallbackRec* callback;
  Widget sliderbar;
  Widget* parent;
  G4bool show;
  short decimal_places;
  G4int initial_value;
  G4int max_value;
  G4int min_value;
  unsigned char orientation;
  unsigned char direction;
};

#endif

#endif
