// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmFramedBox.hh,v 1.4 2001-02-03 18:39:10 johna Exp $
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
  virtual ~G4OpenGLXmFramedBox ();               //destructor

  void AddChild (G4OpenGLXmVWidgetComponent*);
  void AddYourselfTo (G4OpenGLXmVWidgetShell*);

private:
  G4OpenGLXmFramedBox (const G4OpenGLXmFramedBox&);
  G4OpenGLXmFramedBox& operator = (const G4OpenGLXmFramedBox&);
  Widget frame;
};

#endif

#endif
