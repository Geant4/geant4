// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmRadioButton.hh,v 1.2 1999-01-09 16:22:59 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//Radio button class. Inherits from G4OpenGLXmVWidgetComponent

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMRADIOBUTTON_HH
#define G4OPENGLXMRADIOBUTTON_HH

#include "G4OpenGLXmVWidgetComponent.hh"

class G4OpenGLXmRadioButton : public G4OpenGLXmVWidgetComponent
{

public:
  G4OpenGLXmRadioButton (char*,
			 XtCallbackRec*,
			 G4bool,
			 G4int);                    //constructor
  ~G4OpenGLXmRadioButton ();                        //destructor

  void SetName (char*);
  char* GetName ();

  void AddYourselfTo (G4OpenGLXmVWidgetContainer*);

  Widget* GetPointerToParent ();
  Widget* GetPointerToWidget ();

private:
  char* name;
  XtCallbackRec* callback;
  Widget button;
  Widget* parent;
  G4bool default_button;
  G4int number;
};

#endif

#endif
