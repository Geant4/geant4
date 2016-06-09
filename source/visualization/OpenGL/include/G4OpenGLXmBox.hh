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
// $Id: G4OpenGLXmBox.hh,v 1.7 2006/06/29 21:18:24 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// 
//Box container class

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMBOX_HH
#define G4OPENGLXMBOX_HH

#include "G4OpenGLXmVWidgetContainer.hh"
#include "globals.hh"
#include <Xm/Frame.h>
#include <Xm/RowColumn.h>

class G4OpenGLXmVWidgetComponent;
class G4OpenGLXmVWidgetShell;

class G4OpenGLXmBox : public G4OpenGLXmVWidgetContainer
{

public:
  G4OpenGLXmBox (const char* = NULL,
		 G4bool = False);   //constructor
  virtual ~G4OpenGLXmBox ();  //destructor

  void AddChild (G4OpenGLXmVWidgetComponent*);
  void AddYourselfTo (G4OpenGLXmVWidgetShell*);

  Widget* GetPointerToParent ();
  Widget* GetPointerToWidget ();
  
  const char* GetName ();
  void SetName (const char*);

protected:
  const char* name;
  Widget* parent;
  Widget box_row_col;
  G4bool radio;

private:
  G4OpenGLXmBox (const G4OpenGLXmBox&);
  G4OpenGLXmBox& operator = (const G4OpenGLXmBox&);
};

#endif

#endif
