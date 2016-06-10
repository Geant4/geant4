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
// $Id: G4OpenGLXmFourArrowButtons.cc 68043 2013-03-13 14:27:49Z gcosmo $
//
//Four arrow buttons class. Inherits from G4OpenGLXmVWidgetComponent

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmVWidgetComponent.hh"
#include "G4OpenGLXmVWidgetContainer.hh"
#include "G4OpenGLXmFourArrowButtons.hh"
#include <X11/Intrinsic.h>
#include <Xm/Form.h>
#include <Xm/ArrowBG.h>

#include "globals.hh"

G4OpenGLXmFourArrowButtons::G4OpenGLXmFourArrowButtons (XtCallbackRec** c)
: arrow_form(0)
, arrow(0)
, parent(0)
{
  callback = c;
}

G4OpenGLXmFourArrowButtons::~G4OpenGLXmFourArrowButtons ()
{}

void G4OpenGLXmFourArrowButtons::AddYourselfTo (G4OpenGLXmVWidgetContainer* container)
{

  pView = container->GetView ();
  ProcesspView ();

  parent = container->GetPointerToWidget ();

  arrow_form = XtVaCreateManagedWidget 
    ("arrow_form",
     xmFormWidgetClass,
     *parent,
     XmNfractionBase, 3,

     XtNvisual, visual, 
     XtNdepth, depth, 
     XtNcolormap, cmap, 
     XtNborderColor, borcol,
     XtNbackground, bgnd,

     NULL);


///////////////`up' arrow///////////////  
  arrow = XtVaCreateManagedWidget 
    ("up_arrow",
     xmArrowButtonGadgetClass, 
     arrow_form,
     
     XmNtopAttachment, XmATTACH_POSITION,
     XmNtopPosition, 0,
     
     XmNbottomAttachment, XmATTACH_POSITION,
     XmNbottomPosition, 1,
     
     XmNleftAttachment, XmATTACH_POSITION,
     XmNleftPosition, 1,

     XmNrightAttachment, XmATTACH_POSITION,
     XmNrightPosition, 2,

     XmNarrowDirection, XmARROW_UP,
     XmNuserData, True,
     NULL);
  
  XtAddCallbacks (arrow, 
		  XmNactivateCallback, 
		  callback[0]);
  
  XtAddCallbacks (arrow, 
		  XmNarmCallback, 
		  callback[0]);
  
  XtAddCallbacks (arrow, 
		  XmNdisarmCallback, 
		  callback[0]);
  
///////////////`down' arrow///////////////  
  arrow = XtVaCreateManagedWidget 
    ("down_arrow",
     xmArrowButtonGadgetClass, 
     arrow_form,
     
     XmNtopAttachment, XmATTACH_POSITION,
     XmNtopPosition, 2,

     XmNbottomAttachment, XmATTACH_POSITION,
     XmNbottomPosition, 3,

     XmNleftAttachment, XmATTACH_POSITION,
     XmNleftPosition, 1,

     XmNrightAttachment, XmATTACH_POSITION,
     XmNrightPosition, 2,

     XmNarrowDirection, XmARROW_DOWN,
     XmNuserData, False,
     NULL);
  
  XtAddCallbacks (arrow, 
		  XmNactivateCallback, 
		  callback[1]);
  
  XtAddCallbacks (arrow, 
		  XmNarmCallback, 
		  callback[1]);
  
  XtAddCallbacks (arrow, 
		  XmNdisarmCallback, 
		  callback[1]);
  
///////////////`left' arrow///////////////  
  arrow = XtVaCreateManagedWidget 
    ("left_arrow",
     xmArrowButtonGadgetClass, 
     arrow_form,
     
     XmNtopAttachment, XmATTACH_POSITION,
     XmNtopPosition, 1,
     
     XmNbottomAttachment, XmATTACH_POSITION,
     XmNbottomPosition, 2,
     
     XmNleftAttachment, XmATTACH_POSITION,
     XmNleftPosition, 0,
     
     XmNrightAttachment, XmATTACH_POSITION,
     XmNrightPosition, 1,
     
     XmNarrowDirection, XmARROW_LEFT,
     XmNuserData, False,
     NULL);
  
  XtAddCallbacks (arrow, 
		  XmNactivateCallback, 
		  callback[2]);
  
  XtAddCallbacks (arrow, 
		  XmNarmCallback, 
		  callback[2]);
  
  XtAddCallbacks (arrow, 
		  XmNdisarmCallback, 
		  callback[2]);
  
///////////////`right' arrow///////////////  
  arrow = XtVaCreateManagedWidget 
    ("right_arrow",
     xmArrowButtonGadgetClass, 
     arrow_form,
     
     XmNtopAttachment, XmATTACH_POSITION,
     XmNtopPosition, 1,
     
     XmNbottomAttachment, XmATTACH_POSITION,
     XmNbottomPosition, 2,

     XmNleftAttachment, XmATTACH_POSITION,
     XmNleftPosition, 2,

     XmNrightAttachment, XmATTACH_POSITION,
     XmNrightPosition, 3,

     XmNarrowDirection, XmARROW_RIGHT,
     XmNuserData, True,
     NULL);
  
  XtAddCallbacks (arrow, 
		  XmNactivateCallback, 
		  callback[3]);
  
  XtAddCallbacks (arrow, 
		  XmNarmCallback, 
		  callback[3]);
  
  XtAddCallbacks (arrow, 
		  XmNdisarmCallback, 
		  callback[3]);

}

Widget* G4OpenGLXmFourArrowButtons::GetPointerToParent ()
{
  return parent;
}

Widget* G4OpenGLXmFourArrowButtons::GetPointerToWidget () 
{
  return &arrow_form;
}

#endif
