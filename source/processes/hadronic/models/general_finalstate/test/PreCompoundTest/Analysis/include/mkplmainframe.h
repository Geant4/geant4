#ifndef mkplmainframe_h
#define mkplmainframe_h

#include "TApplication.h"
#include "TGFrame.h"
#include "TGLayout.h"
#include "TGMenu.h"
#include "TGFileDialog.h" 
#include "TGTab.h"
#include "TRootEmbeddedCanvas.h"
#include "TCanvas.h"
#include "TGTextView.h"

#include "mkplfilemanager.h"
#include "mkplexpdatamanager.h"
#include "mkplenergyselector.h"
#include "mkplangleselector.h"
#include "mkplrangeselector.h"

using namespace std;


class mkplmainframe : public TGMainFrame
{
private:
  TGMenuBar           * mkpl_menu_bar;
  TGLayoutHints       * mkpl_menu_bar_layout;
  // File Menu
  TGPopupMenu         * mkpl_file_menu;
  TGLayoutHints       * mkpl_file_menu_layout;

  // The Tab
  TGTab               * mkpl_tab;
  TGLayoutHints       * mkpl_tab_layout;

  // Tab1: Summary 
  TGLayoutHints       * mkpl_summary_layout;
  TGTextView          * mkpl_summary_text_viewer;

  // Tab2: Details
  TGLayoutHints       * mkpl_details_layout;
  TGTextView          * mkpl_details_text_viewer;

  // Tab3: Graphics
  TGLayoutHints       * mkpl_graphics_layout;
  TRootEmbeddedCanvas * mkpl_graphics_canvas;


  // the File Manager 
  mkplfilemanager     * mkpl_file_manager;
  mkplexpdatamanager  * mkpl_exp_data_manager;
  // The energies selector
  mkplenergyselector  * mkpl_energy_selector;
  // The angles selector
  mkplangleselector   * mkpl_angle_selector;
  // The angles selector
  mkplrangeselector   * mkpl_range_selector;

public:
  mkplmainframe(const TGWindow *p, UInt_t w, UInt_t h);
  virtual ~mkplmainframe();

  inline virtual void CloseWindow();
  virtual Bool_t ProcessMessage(Long_t msg, Long_t param1, Long_t);

private:
  inline void OpenExpFile(const TString& fname);
  inline void CloseExpFile();
  void MenuCloseFile();
  void MenuOpenFile();

  void ShowDEdata(const int n);
  void ShowDAdata(const int n);
  void ShowDDangles(const int n);
  void ShowDDdata(const int energy, const int angle);
  void ShowDDdataALL(const int energy);
  void ShowDDAranges(const int n);
  void ShowDDAdata(const int energy, const int range);

enum CommandIdentifiers 
  {
    M_FILE_OPEN,
    M_FILE_CLOSE,
    M_FILE_EXIT
  };

  static const char * filetypes[];
};

#include "mkplmainframe.icc"

#endif

