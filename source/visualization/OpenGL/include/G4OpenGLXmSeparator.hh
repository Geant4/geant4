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
//Separator class. Inherits from G4OpenGLXmVWidgetComponent

#if defined (G4VIS_BUILD_OPENGLXM_DRIVER) || defined (G4VIS_USE_OPENGLXM)

#ifndef G4OPENGLXMSEPARATOR_HH
#define G4OPENGLXMSEPARATOR_HH

#include "G4OpenGLXmVWidgetComponent.hh"

class G4OpenGLXmSeparator : public G4OpenGLXmVWidgetComponent
{

public:
  G4OpenGLXmSeparator (unsigned char = XmSINGLE_LINE);           //constructor
  virtual ~G4OpenGLXmSeparator ();                               //destructor

  void AddYourselfTo (G4OpenGLXmVWidgetContainer*);

  Widget* GetPointerToParent ();
  Widget* GetPointerToWidget ();

private:
  G4OpenGLXmSeparator (const G4OpenGLXmSeparator&);
  G4OpenGLXmSeparator& operator = (const G4OpenGLXmSeparator&);
  unsigned char line_type;
  Widget line;
  Widget* parent;
};

#endif

#endif
