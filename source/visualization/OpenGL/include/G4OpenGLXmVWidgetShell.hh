// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmVWidgetShell.hh,v 1.4 2001-02-03 18:39:25 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//Base class for all Motif window widgets (shells)

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMVWIDGETSHELL_HH
#define G4OPENGLXMVWIDGETSHELL_HH

#include "G4OpenGLXmVWidgetObject.hh"

class G4OpenGLXmVWidgetContainer;

class G4OpenGLXmVWidgetShell : public G4OpenGLXmVWidgetObject
{

public:
  G4OpenGLXmVWidgetShell();   //constructor
  virtual ~G4OpenGLXmVWidgetShell();  //destructor

  virtual Widget* GetPointerToWidget() = 0;
  virtual void AddChild (G4OpenGLXmVWidgetContainer*) = 0;
  virtual void Realize () = 0;

private:
};

#endif

#endif
