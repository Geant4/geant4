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
// Andrew Walkden  16th April 1997
// G4OpenGLXmConvenienceRoutines : 
//                       Collection of routines to facilitate
//                       the addition of simple push button boxes,
//                       and slider bars to the control panel.

#include "G4OpenGLXmViewer.hh"

#include <Xm/Form.h>
#include <Xm/ToggleB.h>
#include <Xm/ArrowBG.h>
#include <Xm/RowColumn.h>
#include <Xm/TextF.h>
#include <Xm/Separator.h>
#include <Xm/Scale.h>

#include <sstream>

const G4String G4OpenGLXmViewer::e_str = "";

void G4OpenGLXmViewer::Add_four_arrow_buttons (G4OpenGLXmViewer* pView,
					     XtCallbackRec** arrow_callbacks,
					     Widget* parent_widget) {
  
  Widget arrow_form = XtVaCreateWidget 
    ("arrow_form",
     xmFormWidgetClass,
     *parent_widget,
     XmNfractionBase, 3,
     XtNvisual, pView->vi->visual, 
     XtNdepth, pView->vi->depth, 
     XtNcolormap, pView->cmap, 
     XtNborderColor, pView->borcol,
     XtNbackground, pView->bgnd,
     NULL);
  
  Widget arrow = XtVaCreateManagedWidget 
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
     NULL);
  
  XtVaSetValues (arrow,
		 XmNuserData, True,
		 NULL);
  
  XtAddCallbacks (arrow, 
		 XmNactivateCallback, 
		 arrow_callbacks[0]);
  
  XtAddCallbacks (arrow, 
		 XmNarmCallback, 
		 arrow_callbacks[0]);
  
  XtAddCallbacks (arrow, 
		 XmNdisarmCallback, 
		 arrow_callbacks[0]);
  
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
     NULL);
  
  XtVaSetValues (arrow,
		 XmNuserData, False,
		 NULL);
  
  XtAddCallbacks (arrow, 
		 XmNactivateCallback, 
		 arrow_callbacks[1]);
  
  XtAddCallbacks (arrow, 
		 XmNarmCallback,
		 arrow_callbacks[1]);
  
  XtAddCallbacks (arrow, 
		 XmNdisarmCallback, 
		 arrow_callbacks[1]);
  
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
     NULL);
  
  XtVaSetValues (arrow,
		 XmNuserData, False,
		 NULL);

  XtAddCallbacks (arrow, 
		 XmNactivateCallback, 
		 arrow_callbacks[2]);

  XtAddCallbacks (arrow, 
		 XmNarmCallback, 
		 arrow_callbacks[2]);
  
  XtAddCallbacks (arrow, 
		 XmNdisarmCallback, 
		 arrow_callbacks[2]);
      
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
     NULL);
  
  XtVaSetValues (arrow,
		 XmNuserData, True,
		 NULL);
  
  XtAddCallbacks (arrow, 
		 XmNactivateCallback, 
		 arrow_callbacks[3]);
  
  XtAddCallbacks (arrow, 
		 XmNarmCallback, 
		 arrow_callbacks[3]);

  XtAddCallbacks (arrow, 
		 XmNdisarmCallback, 
		 arrow_callbacks[3]);
  
  XtManageChild (arrow_form);
  
}

void G4OpenGLXmViewer::Add_radio_box (char* label_string,
				    Widget* parent_widget,
				    XtCallbackRec* radio_box_callback,
				    G4int num_buttons,
				    G4int default_button,
				    char* radio_box_name,
				    char** button_names,
				    G4OpenGLXmViewer* pView)
{
  XmString button_str = XmStringCreateLocalized((char*) e_str.c_str());
  // initialise to something to avoid pedantic warning.
  Arg** args;
  args = new Arg* [num_buttons];
  Widget button;

  G4int i;
  for (i = 0; i < num_buttons; i++) {

    args[i] = new Arg[7];
    button_str = XmStringCreateLocalized (button_names[i]);

    XtSetArg (args[i][0], XtNvisual, pView->vi->visual);
    XtSetArg (args[i][1], XtNdepth, pView->vi->depth);
    XtSetArg (args[i][2], XtNcolormap, pView->cmap);
    XtSetArg (args[i][3], XtNborderColor, pView->borcol);
    XtSetArg (args[i][4], XtNbackground, pView->bgnd);
    XtSetArg (args[i][5], XmNlabelString, button_str); 

    if (i == default_button) {
      XtSetArg (args[i][6], XmNset, True);
    } else {
      XtSetArg (args[i][6], XmNset, False);
    }
  }
  
  Widget radio_box = XtVaCreateWidget (radio_box_name,
				       xmRowColumnWidgetClass,
				       *parent_widget,
				       XmNisHomogeneous, False,
				       XmNradioBehavior, True,
				       XmNradioAlwaysOne, True,
				       XmNuserData, pView,
				       XtNvisual, pView->vi->visual,
				       XtNdepth, pView->vi->depth,
				       XtNcolormap, pView->cmap,
				       XtNborderColor, pView->borcol,
				       XtNbackground, pView->bgnd,
				       NULL);

  XmString lab = XmStringCreateLocalized (label_string);

  // Unused!
  //Widget label = XtVaCreateManagedWidget ("radio_label",
  //				  xmLabelWidgetClass,
  //				  radio_box,
  //				  XmNalignment, XmALIGNMENT_CENTER,
  //				  XmNlabelString, lab,
  //				  XtNvisual, pView->vi->visual,
  //				  XtNdepth, pView->vi->depth,
  //				  XtNcolormap, pView->cmap,
  //				  XtNborderColor, pView->borcol,
  //				  XtNbackground, pView->bgnd,
  //				  NULL);

  XmStringFree (lab);

  for (i = 0; i < num_buttons; i++) {
    button = XtCreateManagedWidget (button_names[i],
				    xmToggleButtonWidgetClass,
				    radio_box,
				    args[i],
				    7);
    XtVaSetValues (button,
		   XmNuserData, i,
		   NULL);
    
    XtAddCallbacks (button,
		   XmNarmCallback,
		   radio_box_callback);
  }

  XtManageChild (radio_box);

  XmStringFree (button_str);
  
  for (i = 0; i < num_buttons; i++) {

    delete[] args[i];

  }

  delete[] args;
}  

void G4OpenGLXmViewer::Add_set_field (char* w_name, 
				    char* w_text,
				    Widget* row_col_box,
				    Widget* wid,
				    G4double* val,
				    G4OpenGLXmViewer* pView)
{

  char local_w_text[50];
  strcpy (local_w_text, w_text);

  char label_name[50];
  strcpy (label_name, w_name);
  strcat (label_name, "_label");
  
  char text_field_name[50];
  strcpy (text_field_name, w_name);
  strcat (text_field_name, "_text_field");
  
  XmString local_text = XmStringCreateLocalized (local_w_text);

  // Unused!
  //  Widget label = XtVaCreateManagedWidget (label_name, 
  //				  xmLabelWidgetClass,
  //				  *row_col_box,
  //				  XmNlabelString, local_text,
  //				  XtNvisual, pView->vi->visual, 
  //				  XtNdepth, pView->vi->depth, 
  //				  XtNcolormap, pView->cmap, 
  //				  XtNborderColor, pView->borcol,
  //				  XtNbackground, pView->bgnd,
  //				  NULL);

  XmStringFree (local_text);

  char initial[50];
  snprintf (initial, sizeof initial, "%6.2f", *val);
  
  *wid = XtVaCreateManagedWidget (text_field_name,
				  xmTextFieldWidgetClass,
				  *row_col_box,
				  XmNvalue, (String)initial,
				  XtNvisual, pView->vi->visual, 
				  XtNdepth, pView->vi->depth, 
				  XtNcolormap, pView->cmap, 
				  XtNborderColor, pView->borcol,
				  XtNbackground, pView->bgnd,
				  NULL);

  XtAddCallback (*wid, 
		 XmNvalueChangedCallback,
		 get_double_value_callback,
		 val);

  /* Not actually used - comment out to prevent compiler warnings.
     Instead, just in case it matters, just invoke
     XtVaCreateManagedWidget (JA)
  Widget sep = XtVaCreateManagedWidget ("sep",
					xmSeparatorWidgetClass,
					*row_col_box,
					XmNorientation, XmHORIZONTAL,
					XtNvisual, pView->vi->visual, 
					XtNdepth, pView->vi->depth, 
					XtNcolormap, pView->cmap, 
					XtNborderColor, pView->borcol,
					XtNbackground, pView->bgnd,
					NULL);
  sep = XtVaCreateManagedWidget ("sep",
				 xmSeparatorWidgetClass,
				 *row_col_box,
				 XmNseparatorType, XmNO_LINE,
				 XmNmargin, 1,
				 XmNorientation, XmHORIZONTAL,
				 XtNvisual, pView->vi->visual, 
				 XtNdepth, pView->vi->depth, 
				 XtNcolormap, pView->cmap, 
				 XtNborderColor, pView->borcol,
				 XtNbackground, pView->bgnd,
				 NULL);
  */
  XtVaCreateManagedWidget ("sep",
				 xmSeparatorWidgetClass,
				 *row_col_box,
				 XmNseparatorType, XmNO_LINE,
				 XmNmargin, 1,
				 XmNorientation, XmHORIZONTAL,
				 XtNvisual, pView->vi->visual, 
				 XtNdepth, pView->vi->depth, 
				 XtNcolormap, pView->cmap, 
				 XtNborderColor, pView->borcol,
				 XtNbackground, pView->bgnd,
				 NULL);
}

void G4OpenGLXmViewer::Add_slider_box (char* label_string,
				     G4int num_sliders,
				     char** slider_names,
				     G4OpenGLXmViewer* pView,
				     G4double* min_array,
				     G4double* max_array,
				     G4double* value_array,
				     G4bool* show,
				     short* decimals,
				     unsigned char* orientation,
				     unsigned char* direction,
				     XtCallbackRec** slider_box_callbacks,
				     Widget* parent_widget)
{
  XmString slider_name_str = XmStringCreateLocalized((char*) e_str.c_str());
  // initialise to something to avoid pedantic warning.
  Arg** slider_args;
  slider_args = new Arg*[num_sliders];
  Widget slider;
  G4int j = 0;

  G4int i;
  for (i = 0; i < num_sliders; i++) {
    j = 0; 
    slider_args[i] = new Arg[13];
    slider_name_str = XmStringCreateLtoR (slider_names[i], 
					  XmFONTLIST_DEFAULT_TAG);
    
    XtSetArg (slider_args[i][j], 
	      XtNvisual, pView->vi->visual); j++;
    XtSetArg (slider_args[i][j], 
	      XtNdepth, pView->vi->depth); j++;
    XtSetArg (slider_args[i][j], 
	      XtNcolormap, pView->cmap); j++;
    XtSetArg (slider_args[i][j], 
	      XtNborderColor, pView->borcol); j++;
    XtSetArg (slider_args[i][j], 
	      XtNbackground, pView->bgnd); j++;
    
    XtSetArg (slider_args[i][j], 
	      XmNtitleString, slider_name_str);  j++;
    
    XtSetArg (slider_args[i][j], 
	      XmNmaximum, G4int(max_array[i] * std::pow(10.0, (G4double)decimals[i]))); j++;
    XtSetArg (slider_args[i][j], 
	      XmNminimum, G4int(min_array[i] * std::pow(10.0, (G4double)decimals[i]))); j++;
    XtSetArg (slider_args[i][j], 
	      XmNvalue, G4int(value_array[i] * std::pow(10.0, (G4double)decimals[i]))); j++;
    XtSetArg (slider_args[i][j], 
	      XmNshowValue, show[i]); j++;
    XtSetArg (slider_args[i][j], 
	      XmNdecimalPoints, decimals[i]); j++;
    
    XtSetArg (slider_args[i][j], 
	      XmNorientation, orientation[i]);  j++;
    XtSetArg (slider_args[i][j], 
	      XmNprocessingDirection, direction[i]); j++;

  }

  Widget slider_box = XtVaCreateWidget ("slider_box",
					xmRowColumnWidgetClass,
					*parent_widget,
					XmNisHomogeneous, False,
					XtNvisual, pView->vi->visual,
					XtNdepth, pView->vi->depth,
					XtNcolormap, pView->cmap,
					XtNborderColor, pView->borcol,
					XtNbackground, pView->bgnd,
					NULL);

  XmString lab = XmStringCreateLocalized (label_string);

  // Unused!
  //Widget label = XtVaCreateManagedWidget ("slider_label",
  //				  xmLabelWidgetClass,
  //				  slider_box,
  //				  XmNlabelString, lab,
  //				  XmNalignment, XmALIGNMENT_CENTER,
  //				  XtNvisual, pView->vi->visual,
  //				  XtNdepth, pView->vi->depth,
  //				  XtNcolormap, pView->cmap,
  //				  XtNborderColor, pView->borcol,
  //				  XtNbackground, pView->bgnd,
  //				  NULL);

  XmStringFree (lab);
  
  for (i = 0; i < num_sliders; i++) {
    
    slider = XtCreateManagedWidget (slider_names[i],
				    xmScaleWidgetClass,
				    slider_box,
				    slider_args[i],
				    j);

    XtAddCallbacks (slider,
		   XmNvalueChangedCallback,
		   slider_box_callbacks[i]);
    
    XtAddCallbacks (slider,
		   XmNdragCallback,
		   slider_box_callbacks[i]);
    
  }

  XtManageChild (slider_box);
  XmStringFree (slider_name_str);

  for (i = 0; i < num_sliders; i++) {

    delete[] slider_args[i];

  }

  delete[] slider_args;

}

void G4OpenGLXmViewer::get_double_value_callback (Widget w, 
						XtPointer clientData, 
						XtPointer) 
{
  G4double* val = (G4double*) clientData;
  String string;

  XtVaGetValues (w,
		 XmNvalue, &string,
		 NULL);

//  sscanf (string, "%lg", val);
  std::istringstream iss(string);
  iss >> *val;
}

void G4OpenGLXmViewer::get_text_callback (Widget w, 
					XtPointer clientData, 
					XtPointer) 
{
  char* txt = (char*)clientData;
  String string;

  XtVaGetValues (w,
		 XmNvalue, &string,
		 NULL);

  strcpy(txt, string);
}

G4bool G4OpenGLXmViewer::get_boolean_userData (Widget w)
{
  XtPointer userData;
  XtVaGetValues (w,XmNuserData,&userData,NULL);
  return (G4bool)(((unsigned long)userData)&0xffff);
}

G4int G4OpenGLXmViewer::get_int_userData (Widget w)
{
  XtPointer userData;
  XtVaGetValues (w,XmNuserData,&userData,NULL);
  return (G4int)(unsigned long)userData;
}
