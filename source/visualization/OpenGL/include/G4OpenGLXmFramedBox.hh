// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmFramedBox.hh,v 1.2 1999-01-09 16:22:56 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//Framed box container class
//Inherits from G4OpenGLXmBox

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMFRAMEDBOX_HH
#define G4OPENGLXMFRAMEDBOX_HH

#include "G4OpenGLXmVWidgetContainer.hh"
#include "G4OpenGLXmBox.hh"
#include "globals.hh"
#include <Xm/Frame.h>
#include <Xm/RowColumn.h>

class G4OpenGLXmVWidgetComponent;
class G4OpenGLXmVWidgetShell;

class G4OpenGLXmFramedBox : public G4OpenGLXmBox
{

public:
  G4OpenGLXmFramedBox (char* = NULL, 
		       G4bool = False);   //constructor
  ~G4OpenGLXmFramedBox ();               //destructor

  void AddChild (G4OpenGLXmVWidgetComponent*);
  void AddYourselfTo (G4OpenGLXmVWidgetShell*);

private:
  Widget frame;
};

#endif

#endif
