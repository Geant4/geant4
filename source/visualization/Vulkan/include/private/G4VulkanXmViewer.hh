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
// Alex Vergara Gil  17th October 2026
// Ported to Vulkan from the OpenGL work by Andrew Walkden (10th February 1997).
// G4VulkanXmViewer : Class derived from G4VulkanXViewer, to provide
//                    (Motif) widget Vulkan functionality for GEANT4.

#ifndef G4VULKANXMVIEWER_HH
#define G4VULKANXMVIEWER_HH

#include "G4VulkanXViewer.hh"
#include "globals.hh"

#include <Xm/Xm.h>

class G4VulkanXmTopLevelShell;
class G4VulkanXmRadioButton;
class G4VulkanXmPushButton;
class G4VulkanXmSliderBar;
class G4VulkanXmBox;
class G4VulkanXmTextField;
class G4VulkanXmFramedBox;
class G4VulkanXmFourArrowButtons;
class G4VulkanXmSeparator;

class G4VulkanXmViewer : public G4VulkanXViewer
{
  public:

    G4VulkanXmViewer(G4VulkanSceneHandler& scene);
    virtual ~G4VulkanXmViewer();

  protected:

    virtual void ShowView();
    void ResetView();
    void GetXmConnection();
    virtual void CreateMainWindow();

    XtAppContext app;
    XtWorkProcId workId;
    Widget toplevel, shell, main_win, menubar, style_cascade, actions_cascade, misc_cascade,
      spec_cascade, drawing_style_pullright, background_color_pullright, transparency_pullright,
      antialias_pullright, haloing_pullright, aux_edge_pullright, frame, vkxarea;

    XmString style_str, actions_str, misc_str, spec_str, draw_str, polyhedron_str, wireframe_str,
      hlr_str, hsr_str, hlhsr_str, set_str, rot_str, pan_str, exit_str, quit_str, print_str,
      white_str, black_str, anti_str, trans_str, halo_str, aux_edge_str, bgnd_str, off_str, on_str;

    G4double zoom_high, zoom_low, pan_low, pan_high, dolly_low, dolly_high, fov, rot_sens_limit,
      pan_sens_limit, wob_high, wob_low, wob_sens;

    Pixel bgnd, borcol;

    G4bool pan_right, rotate_right, pan_up, rotate_up;

    XtIntervalId rotation_timer, pan_timer, wobble_timer;

    G4Vector3D original_vp;

    G4int frameNo;
    static const G4String e_str;
    G4String menu_str[37] = {"Style",
                             "style",
                             "Actions",
                             "actions",
                             "Miscellany",
                             "miscellany",
                             "Special",
                             "special",
                             "menubar",
                             "Drawing",
                             "Background color",
                             "Wireframe",
                             "Hidden line removal",
                             "Hidden surface removal",
                             "Hidden line and surface removal",
                             "drawing_style",
                             "White",
                             "Black",
                             "background_color",
                             "Rotation control panel",
                             "Panning control panel",
                             "Set control panel limits",
                             "Miscellany control panel",
                             "Exit to G4Vis>",
                             "Create .eps file",
                             "Transparency",
                             "transparency",
                             "Antialiasing",
                             "antialias",
                             "Haloing",
                             "haloing",
                             "Auxiliary edges",
                             "aux_edge",
                             "Off",
                             "On",
                             "frame",
                             "vkxarea"};

    G4VulkanXmTopLevelShell* fprotation_top;
    G4VulkanXmBox* fprotation_button_box;
    G4VulkanXmRadioButton* fprotation_button1;
    G4VulkanXmRadioButton* fprotation_button2;
    G4VulkanXmBox* fprotation_slider_box;
    G4VulkanXmSliderBar* fprotation_slider;
    G4VulkanXmBox* fprotation_arrow_box;
    G4VulkanXmFourArrowButtons* fprotation_arrow;

    G4VulkanXmTopLevelShell* fppanning_top;
    G4VulkanXmFramedBox* fppanning_box;
    G4VulkanXmFourArrowButtons* fppanning_arrows;
    G4VulkanXmSliderBar* fppanning_slider;
    G4VulkanXmFramedBox* fpzoom_box;
    G4VulkanXmSliderBar* fpzoom_slider;
    G4VulkanXmFramedBox* fpdolly_box;
    G4VulkanXmSliderBar* fpdolly_slider;

    G4VulkanXmTopLevelShell* fpsetting_top;
    G4VulkanXmFramedBox* fpsetting_box;
    G4VulkanXmTextField* fppan_set;
    G4VulkanXmTextField* fprot_set;
    G4VulkanXmTextField* fpzoom_upper;
    G4VulkanXmTextField* fpzoom_lower;
    G4VulkanXmTextField* fpdolly_upper;
    G4VulkanXmTextField* fpdolly_lower;
    G4VulkanXmPushButton* fpok_button;

    G4VulkanXmTopLevelShell* fpmiscellany_top;
    G4VulkanXmFramedBox* fpwobble_box;
    G4VulkanXmPushButton* fpwobble_button;
    G4VulkanXmSliderBar* fpwobble_slider;
    G4VulkanXmFramedBox* fpreset_box;
    G4VulkanXmPushButton* fpreset_button;
    G4VulkanXmFramedBox* fpproj_style_box;
    G4VulkanXmRadioButton* fporthogonal_button;
    G4VulkanXmRadioButton* fpperspective_button;
    G4VulkanXmTextField* fpfov_text;

    G4VulkanXmTopLevelShell* fpprint_top;
    G4VulkanXmFramedBox* fpprint_box;
    G4VulkanXmFramedBox* fpprint_col_box;
    G4VulkanXmFramedBox* fpprint_style_box;
    G4VulkanXmTextField* fpprint_text;
    G4VulkanXmPushButton* fpprint_button;
    G4VulkanXmSeparator* fpprint_line;
    G4VulkanXmRadioButton* fpprint_col_radio1;
    G4VulkanXmRadioButton* fpprint_col_radio2;
    G4VulkanXmRadioButton* fpprint_style_radio1;
    G4VulkanXmRadioButton* fpprint_style_radio2;

  public:

    static void expose_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void resize_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void actions_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void misc_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void Add_set_field(char* widget, char* widget_text, Widget* row_col_box, Widget* wid,
                              G4double* val, G4VulkanXmViewer* pView);

    static void zoom_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void dolly_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void pan_left_right_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void left_right_pan_callback(XtPointer clientData, XtIntervalId* timer_id);

    static void theta_rotation_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void phi_rotation_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void pan_up_down_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void up_down_pan_callback(XtPointer clientData, XtIntervalId* timer_id);

    static void drawing_style_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void background_color_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void set_rot_subject_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void set_rot_sens_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void set_pan_sens_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void set_wob_sens_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void projection_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void wobble_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void reset_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void update_panels_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void wobble_timer_callback(XtPointer clientData, XtIntervalId* timerid);

    static void Add_radio_box(char* label_string, Widget* parent_frame_widget,
                              XtCallbackRec* radio_box_cb, G4int num_buttons, G4int default_button,
                              char* radio_box_name, char** button_names, G4VulkanXmViewer* pView);

    static void Add_four_arrow_buttons(G4VulkanXmViewer* pView, XtCallbackRec** arrow_callbacks,
                                       Widget* parent_widget);

    static void Add_slider_box(char* label_string, G4int num_sliders, char** slider_name,
                               G4VulkanXmViewer* pView, G4double* min_array, G4double* max_array,
                               G4double* value_array, G4bool* show, short* decimals,
                               unsigned char* orientation, unsigned char* direction,
                               XtCallbackRec** slider_box_cb, Widget* parent_frame_widget);

    static void rotate_in_theta(XtPointer clientData, XtIntervalId* timer_id);

    static void rotate_in_phi(XtPointer clientData, XtIntervalId* timer_id);

    static void get_double_value_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void get_text_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void transparency_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void antialias_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void haloing_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void aux_edge_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void set_print_colour_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void set_print_style_callback(Widget w, XtPointer clientData, XtPointer callData);

    static void print_callback(Widget w, XtPointer clientData, XtPointer callData);

    static G4bool get_boolean_userData(Widget w);

    static G4int get_int_userData(Widget w);

    friend class G4VulkanXmVWidgetObject;
    friend class G4VulkanXmViewerMessenger;

  private:

    G4VulkanXmViewer(const G4VulkanXmViewer&);
    G4VulkanXmViewer& operator=(const G4VulkanXmViewer&);
    void UpdateControlPanel();
    // Update the content of the control panel
};

#endif
