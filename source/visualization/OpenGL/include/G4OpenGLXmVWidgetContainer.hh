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
// $Id: G4OpenGLXmVWidgetContainer.hh,v 1.5 2001-07-11 10:08:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//Base class for all Motif container widgets

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMVWIDGETCONTAINER_HH
#define G4OPENGLXMVWIDGETCONTAINER_HH

#include "G4OpenGLXmVWidgetObject.hh"

class G4OpenGLXmVWidgetShell;
class G4OpenGLXmVWidgetComponent;

class G4OpenGLXmVWidgetContainer : public G4OpenGLXmVWidgetObject
{

public:
  G4OpenGLXmVWidgetContainer();           //constructor
  virtual ~G4OpenGLXmVWidgetContainer();  //destructor

  virtual void AddChild (G4OpenGLXmVWidgetComponent*) = 0;
  virtual void AddYourselfTo (G4OpenGLXmVWidgetShell*) = 0;

  virtual Widget* GetPointerToParent () = 0;
  virtual Widget* GetPointerToWidget () = 0;
  
private:

};

#endif

#endif
