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
// $Id: G4OpenGLXmSliderBar.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// 
//Slider bar class. Inherits from G4OpenGLXmVWidgetComponent

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMSLIDERBAR_HH
#define G4OPENGLXMSLIDERBAR_HH

#include "G4OpenGLXmVWidgetComponent.hh"

class G4OpenGLXmSliderBar : public G4OpenGLXmVWidgetComponent
{

public:
  G4OpenGLXmSliderBar (const char* = NULL,     // name of slider bar
		       XtCallbackRec* = NULL,  // callbacks for slider bar
		       G4bool = False,         // show current value if True
		       short = 0,              // decimal places for show value
		       G4double = 0.,          // initial value
		       G4double = 0.,          // max value
		       G4double = 0.,          // min value
		       unsigned char = XmHORIZONTAL,
		       unsigned char = XmMAX_ON_RIGHT); 
                                               //constructor
  virtual ~G4OpenGLXmSliderBar ();             //destructor

  void SetName (const char*);
  void SetShow (G4bool);
  void SetDecimalPlaces (short);
  void SetInitialValue (G4double);
  void SetMaxValue (G4double);
  void SetMinValue (G4double);
  void SetOrientation (unsigned char);
  void SetDirection (unsigned char);

 
  const char* GetName ();
  G4bool GetShow ();
  short GetDecimalPlaces ();
  G4double GetInitialValue ();
  G4double GetMaxValue ();
  G4double GetMinValue ();
  unsigned char GetOrientation ();
  unsigned char GetDirection ();

  void AddYourselfTo (G4OpenGLXmVWidgetContainer*);

  Widget* GetPointerToParent ();
  Widget* GetPointerToWidget ();

private:
  G4OpenGLXmSliderBar (const G4OpenGLXmSliderBar&);
  G4OpenGLXmSliderBar& operator = (const G4OpenGLXmSliderBar&);
  const char* name;
  XtCallbackRec* callback;
  Widget sliderbar;
  Widget* parent;
  G4bool show;
  short decimal_places;
  G4int initial_value;
  G4int max_value;
  G4int min_value;
  unsigned char orientation;
  unsigned char direction;
};

#endif

#endif
