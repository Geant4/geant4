#ifndef mkplmakeplots_h
#define mkplmakeplots_h

#include "TApplication.h"
#include "TGFrame.h"
#include "TGLayout.h"
#include "TGMenu.h"
#include "TGFileDialog.h" 
#include "TGTab.h"
#include "TRootEmbeddedCanvas.h"
#include "TCanvas.h"
#include "TGTextView.h"
#include "TGButton.h"

#include "mkplfilemanager.h"
#include "mkplexpdatamanager.h"
#include "mkplsimmanager.h"
#include "mkplcomparisonselector.h"

using namespace std;


class mkplmakeplots : public TGMainFrame
{
private:
  TGMenuBar              * mkpl_menu_bar;
  TGLayoutHints          * mkpl_menu_bar_layout;
  // File Menu
  TGPopupMenu            * mkpl_file_menu;
  TGLayoutHints          * mkpl_file_menu_layout;

  // The Tab
  TGTab                  * mkpl_tab;
  TGLayoutHints          * mkpl_tab_layout;

  // Tab1: Options
  TGHorizontalFrame      * mkpl_options;
  TGLayoutHints          * mkpl_options_layout;

  TGVerticalFrame        * mkpl_info_group;
  TGVerticalFrame        * mkpl_actions_group;
  TGLayoutHints          * mkpl_info_group_layout;
  TGLayoutHints          * mkpl_actions_group_layout;
  TGHorizontalFrame      * mkpl_info_exp_group;
  TGHorizontalFrame      * mkpl_info_sim_group;
  TGHorizontalFrame      * mkpl_info_comp_group;
  TGLayoutHints          * mkpl_info_exp_group_layout;
  TGLayoutHints          * mkpl_info_sim_group_layout;
  TGLayoutHints          * mkpl_info_comp_group_layout;
  TGLabel                * mkpl_info_exp_label;
  TGLabel                * mkpl_info_sim_label;
  TGLabel                * mkpl_info_comp_label;
  TGLayoutHints          * mkpl_info_exp_label_layout;
  TGLayoutHints          * mkpl_info_sim_label_layout;
  TGLayoutHints          * mkpl_info_comp_label_layout;
  TGLayoutHints          * mkpl_info_exp_text_layout;
  TGLayoutHints          * mkpl_info_sim_text_layout;
  TGLayoutHints          * mkpl_info_comp_text_layout;
  TGTextView             * mkpl_info_exp_text;
  TGTextView             * mkpl_info_sim_text;
  TGTextView             * mkpl_info_comp_text;
  TGCheckButton          * mkpl_comparisons_button;
  TGCheckButton          * mkpl_tests_button;
  TGCheckButton          * mkpl_components_button;
  TGTextButton           * mkpl_fill_button;

  // Tab2: Graphics
  TRootEmbeddedCanvas    * mkpl_graphics;
  TGLayoutHints          * mkpl_graphics_layout;

  // File Manager
  mkplfilemanager        * mkpl_file_manager;
  mkplexpdatamanager     * mkpl_exp_data_manager;
  mkplsimmanager         * mkpl_sim_data_manager;

  // Comparison selector
  mkplcomparisonselector * mkpl_comparison_selector;

public:
  mkplmakeplots(const TGWindow * p, UInt_t w, UInt_t h);
  virtual ~mkplmakeplots();

  inline virtual void CloseWindow();
  virtual Bool_t ProcessMessage(Long_t msg, Long_t param1, Long_t param2);

private:
  inline void OpenExpFile(const TString& fname);
  inline void CloseExpFile();
  inline void OpenSimFile(const TString& fname);
  inline void CloseSimFile();
  void MenuCloseExpFile();
  void MenuCloseSimFile();
  void MenuOpenExpFile();
  void MenuOpenSimFile();
  void MenuSaveHistograms();
  void MenuSaveEPS();

  void ExpFileShowInfo();
  void SimFileShowInfo();
  void FillHistograms();

  void PlotComparison(const Long_t p1, const Long_t p2);
  void PlotTest(const Long_t p1);

enum CommandIdentifiers
  {
    M_FILE_EXP_OPEN,
    M_FILE_EXP_CLOSE,
    M_FILE_SIM_OPEN,
    M_FILE_SIM_CLOSE,
    M_FILE_SAVE_HISTOGRAMS,
    M_FILE_SAVE_EPS,
    M_FILE_EXIT
  };

  static const char * filetypes[];
};

#include "mkplmakeplots.icc"

#endif
