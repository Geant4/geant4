// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmVWidgetComponent.hh,v 1.3 1999-12-15 14:54:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//Base class for all Motif component widgets

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMVWIDGETCOMPONENT_HH
#define G4OPENGLXMVWIDGETCOMPONENT_HH

#include "G4OpenGLXmVWidgetObject.hh"

class G4OpenGLXmVWidgetContainer;

class G4OpenGLXmVWidgetComponent : public G4OpenGLXmVWidgetObject
{

public:
  G4OpenGLXmVWidgetComponent();   //constructor
  ~G4OpenGLXmVWidgetComponent();  //destructor

  virtual void AddYourselfTo (G4OpenGLXmVWidgetContainer*) = 0;

  virtual Widget* GetPointerToParent () = 0;
  virtual Widget* GetPointerToWidget () = 0;

private:
  
};

#endif

#endif
