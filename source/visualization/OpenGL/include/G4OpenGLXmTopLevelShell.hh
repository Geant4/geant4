// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmTopLevelShell.hh,v 1.1 1999-01-07 16:14:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//Top level shell class

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMTOPLEVELSHELL_HH
#define G4OPENGLXMTOPLEVELSHELL_HH

#include "G4OpenGLXmVWidgetShell.hh"

class G4OpenGLXmVWidgetContainer;

class G4OpenGLXmTopLevelShell : public G4OpenGLXmVWidgetShell
{

public:
  G4OpenGLXmTopLevelShell(G4OpenGLXmView*, char*);   //constructor
  ~G4OpenGLXmTopLevelShell();  //destructor

  void AddChild (G4OpenGLXmVWidgetContainer*);
  void Realize ();

  Widget* GetPointerToWidget ();
  char* GetName ();

private:
  char* name;
  Widget toplevel;
  Widget top_box;
  Widget frame;
};

#endif

#endif
