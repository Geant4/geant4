//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4OpenGLXmVWidgetShell.hh,v 1.5 2001-07-11 10:08:52 gunter Exp $
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
