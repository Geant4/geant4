// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmViewer.hh,v 1.2 1999-11-11 15:38:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// G4OpenGLXmViewer : Class derived from G4OpenGLXViewer, to provide
//                    (Motif) widget OpenGL functionality for GEANT4.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMVIEWER_HH
#define G4OPENGLXMVIEWER_HH

#include "G4OpenGLXViewer.hh"
#include "G4OpenGLSceneHandler.hh"
#include "globals.hh"
#include "g4rw/tvordvec.h"

#include <Xm/Xm.h>
#include <Xm/MainW.h>
#include <Xm/Label.h>
#include <Xm/Form.h>
#include <Xm/RowColumn.h>
#include <Xm/PushB.h>
#include <Xm/ToggleB.h>
#include <Xm/ArrowBG.h>
#include <Xm/CascadeB.h>
#include <Xm/DrawingA.h>
#include <Xm/Frame.h>
#include <X11/keysym.h>
#include <Xm/XmStrDefs.h>
#include <X11/Xmu/StdCmap.h>
#include <Xm/TextF.h>
#include <Xm/ToggleB.h>
#include <X11/Shell.h>
#include <Xm/Separator.h>
#include <Xm/Scale.h>

class G4OpenGLXmTopLevelShell;
class G4OpenGLXmRadioButton;
class G4OpenGLXmPushButton;
class G4OpenGLXmSliderBar;
class G4OpenGLXmBox;
class G4OpenGLXmTextField;
class G4OpenGLXmFramedBox;
class G4OpenGLXmFourArrowButtons;
class G4OpenGLXmSeparator;

class G4OpenGLXmViewer: public G4OpenGLXViewer {
  
public:
  G4OpenGLXmViewer (G4OpenGLSceneHandler& scene);
  ~G4OpenGLXmViewer ();
  
protected:
  virtual void ShowView ();
  void GetXmConnection ();
  virtual void CreateMainWindow ();

  XtAppContext                      app;
  XtWorkProcId                      workId;
  Widget                            toplevel, 
                                    shell,
                                    main_win, 
                                    menubar,
                                    style_cascade,
                                    actions_cascade,
                                    misc_cascade,
                                    spec_cascade,
                                    rep_style_pullright, 
                                    drawing_style_pullright, 
                                    background_color_pullright, 
                                    transparency_pullright, 
                                    antialias_pullright, 
                                    haloing_pullright, 
                                    frame, 
                                    glxarea;
  
  XmString                          style_str,
                                    actions_str,
                                    misc_str,
                                    spec_str,
                                    rep_str,
                                    draw_str,
                                    polyhedron_str,
                                    nurbs_str,
                                    wireframe_str,
                                    hlr_str,
                                    hsr_str,
                                    hlhsr_str,
                                    set_str,
				    rot_str,
				    pan_str,
				    exit_str,
				    quit_str,
                                    print_str,
				    white_str,
				    black_str,
				    anti_str,
				    trans_str,
				    halo_str,
				    bgnd_str,
				    off_str,
				    on_str;

  G4double                          zoom_high,
                                    zoom_low,
                                    pan_low,
                                    pan_high,
                                    dolly_low,
                                    dolly_high,
                                    fov,
                                    rot_sens_limit,
                                    pan_sens_limit,
                                    rot_sens,
                                    pan_sens,
                                    wob_sens;

  Pixel                             bgnd, 
                                    borcol;

  G4bool                            pan_right,
                                    rotate_right,
                                    pan_up,
                                    rotate_up;

  XtIntervalId                      rotation_timer,
                                    pan_timer,
                                    wobble_timer;

  G4Vector3D                        original_vp;

  G4int                             frameNo;

  G4OpenGLXmTopLevelShell*          fprotation_top;
  G4OpenGLXmBox*                    fprotation_button_box;
  G4OpenGLXmRadioButton*            fprotation_button1;
  G4OpenGLXmRadioButton*            fprotation_button2;
  G4OpenGLXmBox*                    fprotation_slider_box;
  G4OpenGLXmSliderBar*              fprotation_slider;
  G4OpenGLXmBox*                    fprotation_arrow_box;
  G4OpenGLXmFourArrowButtons*       fprotation_arrow;

  G4OpenGLXmTopLevelShell*          fppanning_top; 
  G4OpenGLXmFramedBox*              fppanning_box;
  G4OpenGLXmFourArrowButtons*       fppanning_arrows;
  G4OpenGLXmSliderBar*              fppanning_slider;
  G4OpenGLXmFramedBox*              fpzoom_box;
  G4OpenGLXmSliderBar*              fpzoom_slider;
  G4OpenGLXmFramedBox*              fpdolly_box;
  G4OpenGLXmSliderBar*              fpdolly_slider;

  G4OpenGLXmTopLevelShell*          fpsetting_top;
  G4OpenGLXmFramedBox*              fpsetting_box;
  G4OpenGLXmTextField*              fppan_set;
  G4OpenGLXmTextField*              fprot_set;
  G4OpenGLXmTextField*              fpzoom_upper;
  G4OpenGLXmTextField*              fpzoom_lower;
  G4OpenGLXmTextField*              fpdolly_upper;
  G4OpenGLXmTextField*              fpdolly_lower;
  G4OpenGLXmPushButton*             fpok_button;

  G4OpenGLXmTopLevelShell*          fpmiscellany_top;
  G4OpenGLXmFramedBox*              fpwobble_box;
  G4OpenGLXmPushButton*             fpwobble_button;
  G4OpenGLXmSliderBar*              fpwobble_slider;
  G4OpenGLXmFramedBox*              fpreset_box;
  G4OpenGLXmPushButton*             fpreset_button;
  G4OpenGLXmFramedBox*              fpproj_style_box;
  G4OpenGLXmRadioButton*            fporthogonal_button;
  G4OpenGLXmRadioButton*            fpperspective_button;
  G4OpenGLXmTextField*              fpfov_text;

  G4OpenGLXmTopLevelShell*          fpprint_top;
  G4OpenGLXmFramedBox*              fpprint_box;
  G4OpenGLXmFramedBox*              fpprint_col_box;
  G4OpenGLXmFramedBox*              fpprint_style_box;
  G4OpenGLXmTextField*              fpprint_text;
  G4OpenGLXmPushButton*             fpprint_button;
  G4OpenGLXmSeparator*              fpprint_line;
  G4OpenGLXmRadioButton*            fpprint_col_radio1;
  G4OpenGLXmRadioButton*            fpprint_col_radio2;
  G4OpenGLXmRadioButton*            fpprint_style_radio1;
  G4OpenGLXmRadioButton*            fpprint_style_radio2;

public:

  static void expose_callback (Widget w, 
			       XtPointer clientData, 
			       XtPointer callData);
  
  static void resize_callback (Widget w, 
			       XtPointer clientData, 
			       XtPointer callData);
  
  static void actions_callback (Widget w, 
				XtPointer clientData, 
				XtPointer callData);
  
  static void misc_callback (Widget w, 
			     XtPointer clientData, 
			     XtPointer callData);
  
  static void Add_set_field (char* widget, 
			     char* widget_text, 
			     Widget* row_col_box, 
			     Widget* wid, 
			     G4double* val,
			     G4OpenGLXmViewer* pView);
  
  static void zoom_callback (Widget w, 
			     XtPointer clientData, 
			     XtPointer callData);
  
  static void dolly_callback (Widget w, 
			      XtPointer clientData, 
			      XtPointer callData);
  
  static void pan_left_right_callback (Widget w, 
				       XtPointer clientData, 
				       XtPointer callData);
  
  static void left_right_pan_callback (XtPointer clientData, 
				       XtIntervalId* timer_id); 
  
  static void theta_rotation_callback (Widget w, 
				       XtPointer clientData, 
				       XtPointer callData);
  
  static void phi_rotation_callback (Widget w, 
				     XtPointer clientData, 
				     XtPointer callData);
  
  static void pan_up_down_callback (Widget w, 
				    XtPointer clientData, 
				    XtPointer callData);
  
  static void up_down_pan_callback (XtPointer clientData, 
				    XtIntervalId* timer_id); 
  
  static void drawing_style_callback (Widget w, 
				      XtPointer clientData, 
				      XtPointer callData);
  
  static void rep_style_callback (Widget w, 
				  XtPointer clientData, 
				  XtPointer callData);
  
  static void background_color_callback (Widget w, 
					 XtPointer clientData, 
					 XtPointer callData);
  
  static void set_rot_subject_callback (Widget w, 
					XtPointer clientData, 
					XtPointer callData);
  
  static void set_rot_sens_callback (Widget w, 
				     XtPointer clientData, 
				     XtPointer callData);
  
  static void set_pan_sens_callback (Widget w, 
				     XtPointer clientData, 
				     XtPointer callData);
  
  static void set_wob_sens_callback (Widget w, 
				     XtPointer clientData, 
				     XtPointer callData);
  
  static void projection_callback (Widget w, 
				   XtPointer clientData, 
				   XtPointer callData);
  
  static void wobble_callback (Widget w, 
			       XtPointer clientData, 
			       XtPointer callData);
  
  static void reset_callback (Widget w, 
			      XtPointer clientData, 
			      XtPointer callData);
  
  static void update_panels_callback (Widget w, 
				      XtPointer clientData, 
				      XtPointer callData);
  
  static void wobble_timer_callback (XtPointer clientData, 
				     XtIntervalId* timerid);
  
  static void Add_radio_box (char* label_string,
			     Widget* parent_frame_widget,
			     XtCallbackRec* radio_box_cb,
			     G4int num_buttons,
			     G4int default_button,
			     char* radio_box_name,
			     char** button_names,
			     G4OpenGLXmViewer* pView);
  
  static void Add_four_arrow_buttons (G4OpenGLXmViewer* pView,
				      XtCallbackRec** arrow_callbacks,
				      Widget* parent_widget);
  
  static void Add_slider_box (char* label_string,
			      G4int num_sliders,
			      char** slider_name,
			      G4OpenGLXmViewer* pView,
			      G4double* min_array,
			      G4double* max_array,
			      G4double* value_array,
			      G4bool* show,
			      short* decimals,
			      unsigned char* orientation,
			      unsigned char* direction,
			      XtCallbackRec** slider_box_cb,
			      Widget* parent_frame_widget);
  
  static void rotate_in_theta (XtPointer clientData,
			       XtIntervalId* timer_id);
  
  static void rotate_in_phi (XtPointer clientData,
			     XtIntervalId* timer_id);
  
  static void get_double_value_callback (Widget w, 
					 XtPointer clientData, 
					 XtPointer callData); 
  
  static void get_text_callback (Widget w, 
				 XtPointer clientData, 
				 XtPointer callData); 
  
  static void transparency_callback (Widget w, 
				     XtPointer clientData, 
				     XtPointer callData); 
  
  static void antialias_callback (Widget w, 
				  XtPointer clientData, 
				  XtPointer callData); 

  static void haloing_callback (Widget w, 
				XtPointer clientData, 
				XtPointer callData); 

  static void set_print_colour_callback (Widget w, 
					 XtPointer clientData, 
					 XtPointer callData); 

  static void set_print_style_callback (Widget w, 
					XtPointer clientData, 
					XtPointer callData); 

  static void print_callback (Widget w, 
			      XtPointer clientData, 
			      XtPointer callData); 

  static G4bool get_boolean_userData (Widget w);

  static G4int  get_int_userData (Widget w);

friend class G4OpenGLXmVWidgetObject;
  
};

#endif

#endif
  
