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
//
// 
//Radio button class. Inherits from G4OpenGLXmVWidgetComponent

#if defined (G4VIS_BUILD_OPENGLXM_DRIVER) || defined (G4VIS_USE_OPENGLXM)

#ifndef G4OPENGLXMRADIOBUTTON_HH
#define G4OPENGLXMRADIOBUTTON_HH

#include "G4OpenGLXmVWidgetComponent.hh"

class G4OpenGLXmRadioButton : public G4OpenGLXmVWidgetComponent
{

public:
  G4OpenGLXmRadioButton (const char*,
			 XtCallbackRec*,
			 G4bool,
			 G4int);                    //constructor
  virtual ~G4OpenGLXmRadioButton ();                //destructor

  void SetName (const char*);
  const char* GetName ();

  void AddYourselfTo (G4OpenGLXmVWidgetContainer*);

  Widget* GetPointerToParent ();
  Widget* GetPointerToWidget ();

private:
  G4OpenGLXmRadioButton (const G4OpenGLXmRadioButton&);
  G4OpenGLXmRadioButton& operator = (const G4OpenGLXmRadioButton&);
  const char* name;
  XtCallbackRec* callback;
  Widget button;
  Widget* parent;
  G4bool default_button;
  G4int number;
};

#endif

#endif
