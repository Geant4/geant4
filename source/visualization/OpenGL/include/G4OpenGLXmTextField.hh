// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmTextField.hh,v 1.5 2001-03-07 14:56:17 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//Text field class. Inherits from G4OpenGLXmVWidgetComponent

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMTEXTFIELD_HH
#define G4OPENGLXMTEXTFIELD_HH

#include "G4OpenGLXmVWidgetComponent.hh"
//#include "G4OpenGLXmConvenienceRoutines.hh"

class G4OpenGLXmTextField : public G4OpenGLXmVWidgetComponent
{

public:
  G4OpenGLXmTextField (const char*,G4double*);   //constructor
  G4OpenGLXmTextField (const char*,const char*); //constructor
  virtual ~G4OpenGLXmTextField ();               //destructor

  void SetName (const char*);
  const char* GetName ();

  void SetValue (G4double);
  void SetValue (const char*);
  const char* GetValue ();

  void AddYourselfTo (G4OpenGLXmVWidgetContainer*);

  Widget* GetPointerToParent ();
  Widget* GetPointerToWidget ();

private:
  G4OpenGLXmTextField (const G4OpenGLXmTextField&);
  G4OpenGLXmTextField& operator = (const G4OpenGLXmTextField&);
  const char* name;
  void* value;
  G4bool text;
  char* initial;
  Widget text_label;
  Widget text_field;
  Widget* parent;
};

#endif

#endif
