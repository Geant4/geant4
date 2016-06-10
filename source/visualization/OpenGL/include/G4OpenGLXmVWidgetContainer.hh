//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4OpenGLXmVWidgetContainer.hh 66373 2012-12-18 09:41:34Z gcosmo $
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
