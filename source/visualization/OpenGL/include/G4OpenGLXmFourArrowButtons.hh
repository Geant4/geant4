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
// $Id: G4OpenGLXmFourArrowButtons.hh,v 1.5 2001-07-11 10:08:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//Four arrow buttons class. Inherits from G4OpenGLXmVWidgetComponent

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMFOURARROWBUTTONS_HH
#define G4OPENGLXMFOURARROWBUTTONS_HH

#include "G4OpenGLXmVWidgetComponent.hh"

class G4OpenGLXmFourArrowButtons : public G4OpenGLXmVWidgetComponent
{

public:
  G4OpenGLXmFourArrowButtons (XtCallbackRec** = NULL); // array of 4 callbacks
                                                       //constructor
  virtual ~G4OpenGLXmFourArrowButtons ();              //destructor

  void SetName (char*);
 
  char* GetName ();

  void AddYourselfTo (G4OpenGLXmVWidgetContainer*);

  Widget* GetPointerToParent ();
  Widget* GetPointerToWidget ();

private:
  XtCallbackRec** callback;
  Widget arrow_form;
  Widget arrow;
  Widget* parent;

private:
  G4OpenGLXmFourArrowButtons (const G4OpenGLXmFourArrowButtons&);
  G4OpenGLXmFourArrowButtons& operator = (const G4OpenGLXmFourArrowButtons&);
};

#endif

#endif
