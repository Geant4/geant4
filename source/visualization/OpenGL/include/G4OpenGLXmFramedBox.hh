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
// $Id: G4OpenGLXmFramedBox.hh,v 1.6 2001-07-11 10:08:51 gunter Exp $
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
  G4OpenGLXmFramedBox (const char* = NULL, 
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
