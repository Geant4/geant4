// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmFourArrowButtons.hh,v 1.2 1999-01-09 16:22:56 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//Four arrow buttons class. Inherits from G4OpenGLXmVWidgetComponent

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMFOURARROWBUTTONS_HH
#define G4OPENGLXMFOURARROWBUTTONS_HH

#include "G4OpenGLXmVWidgetComponent.hh"

class G4OpenGLXmFourArrowButtons : public G4OpenGLXmVWidgetComponent
{

public:
  G4OpenGLXmFourArrowButtons (XtCallbackRec** = NULL); // array of 4 callbacks
                                                       //constructor
  ~G4OpenGLXmFourArrowButtons ();                      //destructor

  void SetName (char*);
 
  char* GetName ();

  void AddYourselfTo (G4OpenGLXmVWidgetContainer*);

  Widget* GetPointerToParent ();
  Widget* GetPointerToWidget ();

private:
  XtCallbackRec** callback;
  Widget arrow_form;
  Widget arrow;
  Widget* parent;
};

#endif

#endif
