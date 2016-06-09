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
// $Id: G4OpenGLXmFramedBox.hh,v 1.7 2006/06/29 21:18:28 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
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
