// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmPushButton.hh,v 1.4 2001-02-03 18:39:11 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//Push button class. Inherits from G4OpenGLXmVWidgetComponent

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMPUSHBUTTON_HH
#define G4OPENGLXMPUSHBUTTON_HH

#include "G4OpenGLXmVWidgetComponent.hh"

class G4OpenGLXmPushButton : public G4OpenGLXmVWidgetComponent
{

public:
  G4OpenGLXmPushButton (char* = NULL,
			XtCallbackRec* = NULL); //constructor
  virtual ~G4OpenGLXmPushButton ();             //destructor

  void SetName (char*);
  char* GetName ();

  void AddYourselfTo (G4OpenGLXmVWidgetContainer*);

  Widget* GetPointerToParent ();
  Widget* GetPointerToWidget ();

private:
  G4OpenGLXmPushButton (const G4OpenGLXmPushButton&);
  G4OpenGLXmPushButton& operator = (const G4OpenGLXmPushButton&);
  char* name;
  XtCallbackRec* callback;
  Widget button;
  Widget* parent;
};

#endif

#endif
