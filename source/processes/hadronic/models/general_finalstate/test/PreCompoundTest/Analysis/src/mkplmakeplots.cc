#include "mkplmakeplots.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TGLabel.h"
#include "TGMsgBox.h"
#include "TStyle.h"
#include "TFrame.h"

#include "mkplmessages.h"
#include "mkplcomparisonhistograms.h"
#include "mkplcomparisonagreement.h"

const char * mkplmakeplots::filetypes[] = 
  {
    "Root files", "*.root",
    "All files", "*",
    0, 0
  };

mkplmakeplots::mkplmakeplots(const TGWindow * p, UInt_t w, UInt_t h)
  : TGMainFrame(p,w,h), mkpl_exp_data_manager(0), mkpl_sim_data_manager(0),
    mkpl_comparison_selector(0)
{
  mkpl_file_manager = new mkplfilemanager();

  // File Menu
  mkpl_file_menu_layout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
  mkpl_file_menu = new TGPopupMenu(fClient->GetRoot());
  mkpl_file_menu->AddEntry("&Open Experimental File", M_FILE_EXP_OPEN);
  mkpl_file_menu->AddEntry("&Close Experimental File", M_FILE_EXP_CLOSE);
  mkpl_file_menu->AddSeparator();
  mkpl_file_menu->AddEntry("&Open Simulation File", M_FILE_SIM_OPEN);
  mkpl_file_menu->AddEntry("&Close Simulation File", M_FILE_SIM_CLOSE);
  mkpl_file_menu->AddSeparator();
  mkpl_file_menu->AddEntry("&Save Histograms", M_FILE_SAVE_HISTOGRAMS);
  mkpl_file_menu->AddEntry("&Save EPS", M_FILE_SAVE_EPS);
  mkpl_file_menu->AddSeparator();
  mkpl_file_menu->AddEntry("E&xit", M_FILE_EXIT);
  mkpl_file_menu->Associate(this);

  mkpl_file_menu->DisableEntry(M_FILE_SAVE_HISTOGRAMS);
  mkpl_file_menu->DisableEntry(M_FILE_SAVE_EPS);

  // Menu Bar
  mkpl_menu_bar_layout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 1, 1);
  mkpl_menu_bar = new TGMenuBar(this, 1, 1, kHorizontalFrame);
  mkpl_menu_bar->AddPopup("&File", mkpl_file_menu, mkpl_file_menu_layout);

  
  this->AddFrame(mkpl_menu_bar,mkpl_menu_bar_layout);
  this->MapSubwindows();

  // The TAB
  mkpl_tab = new TGTab(this, 300, 300);
  mkpl_tab_layout = new TGLayoutHints(kLHintsBottom | kLHintsExpandX | kLHintsExpandY, 2, 2, 5, 1);

  // Tab: Options
  //=============
  TGCompositeFrame * compFrame = mkpl_tab->AddTab("Options");
  mkpl_options_layout = new TGLayoutHints(kLHintsTop | kLHintsLeft |
					  kLHintsExpandX | kLHintsExpandY,
					  5, 5, 5, 5);
  mkpl_options = new TGHorizontalFrame(compFrame,100,100);
  compFrame->AddFrame(mkpl_options,mkpl_options_layout);

  // info
  //-----
  mkpl_info_group = new TGVerticalFrame(mkpl_options, 100, 100);
  mkpl_info_group_layout = new TGLayoutHints(kLHintsTop | kLHintsLeft | 
					     kLHintsExpandX | kLHintsExpandY, 
					     5, 5, 5, 5);

  mkpl_info_exp_group = new TGHorizontalFrame(mkpl_info_group, 100,100);
  mkpl_info_exp_group_layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX |
						 kLHintsExpandX | kLHintsExpandY,
						 5, 5, 5, 5);

  mkpl_info_exp_label = new TGLabel(mkpl_info_exp_group, new TGHotString("Experimental: "));
  mkpl_info_exp_label_layout = new TGLayoutHints(kLHintsTop | kLHintsLeft |
						 kLHintsExpandY,
						 5, 5, 5, 5);
  mkpl_info_exp_label->SetTextJustify(kTextRight);
  mkpl_info_exp_label->Resize(mkpl_info_exp_label->GetDefaultSize());
  mkpl_info_exp_text = new TGTextView(mkpl_info_exp_group, 250, 150);
  mkpl_info_exp_text_layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX |
						kLHintsExpandX | kLHintsExpandY,
						5, 5, 5, 5);

  mkpl_info_sim_group = new TGHorizontalFrame(mkpl_info_group, 100,100);
  mkpl_info_sim_group_layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX |
						 kLHintsExpandX | kLHintsExpandY,
						 5, 5, 5, 5);

  mkpl_info_sim_label = new TGLabel(mkpl_info_sim_group, new TGHotString("  Simulation: "));
  mkpl_info_sim_label_layout = new TGLayoutHints(kLHintsTop | kLHintsLeft |
						 kLHintsExpandY,
						 5, 5, 5, 5);
  mkpl_info_sim_label->SetTextJustify(kTextRight);
  mkpl_info_sim_label->Resize(mkpl_info_exp_label->GetDefaultSize());
  mkpl_info_sim_text = new TGTextView(mkpl_info_sim_group, 250, 150);
  mkpl_info_sim_text_layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX |
						kLHintsExpandX | kLHintsExpandY,
						5, 5, 5, 5);
  

  mkpl_info_comp_group = new TGHorizontalFrame(mkpl_info_group, 100,100);
  mkpl_info_comp_group_layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX |
						  kLHintsExpandX | kLHintsExpandY,
						  5, 5, 5, 5);


  mkpl_info_comp_label = new TGLabel(mkpl_info_comp_group, new TGHotString(" Comparison: "));
  mkpl_info_comp_label_layout = new TGLayoutHints(kLHintsTop | kLHintsLeft |
						  kLHintsExpandY,
						  5, 5, 5, 5);
  mkpl_info_comp_label->SetTextJustify(kTextRight);
  mkpl_info_comp_label->Resize(mkpl_info_exp_label->GetDefaultSize());
  mkpl_info_comp_text = new TGTextView(mkpl_info_comp_group, 250,150);
  mkpl_info_comp_text_layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX |
						 kLHintsExpandX | kLHintsExpandY,
						 5, 5, 5, 5);
  

  mkpl_info_exp_group->AddFrame(mkpl_info_exp_label,mkpl_info_exp_label_layout);
  mkpl_info_exp_group->AddFrame(mkpl_info_exp_text,mkpl_info_exp_text_layout);

  mkpl_info_sim_group->AddFrame(mkpl_info_sim_label,mkpl_info_sim_label_layout);
  mkpl_info_sim_group->AddFrame(mkpl_info_sim_text,mkpl_info_sim_text_layout);

  mkpl_info_comp_group->AddFrame(mkpl_info_comp_label,mkpl_info_comp_label_layout);
  mkpl_info_comp_group->AddFrame(mkpl_info_comp_text,mkpl_info_comp_text_layout);


  mkpl_info_group->AddFrame(mkpl_info_exp_group,mkpl_info_exp_group_layout);
  mkpl_info_group->AddFrame(mkpl_info_sim_group,mkpl_info_sim_group_layout);
  mkpl_info_group->AddFrame(mkpl_info_comp_group,mkpl_info_comp_group_layout);

  mkpl_options->AddFrame(mkpl_info_group, mkpl_info_group_layout);


  mkpl_info_group->Resize(mkpl_info_group->GetDefaultSize());

  // actions
  //--------
  mkpl_actions_group = new TGVerticalFrame(mkpl_options, 10,10);
  mkpl_actions_group_layout = new TGLayoutHints(kLHintsTop | kLHintsRight | 
						kLHintsExpandY, 
						5, 5, 5, 5);
  mkpl_options->AddFrame(mkpl_actions_group, mkpl_actions_group_layout);


  mkpl_comparisons_button = new TGCheckButton(mkpl_actions_group, "Comparisons",-1);
  mkpl_tests_button = new TGCheckButton(mkpl_actions_group, "Tests",-1);
  mkpl_components_button = new TGCheckButton(mkpl_actions_group, "Components",-1); 
  mkpl_fill_button = new TGTextButton(mkpl_actions_group, "Fill Histograms", 81);
  mkpl_fill_button->Resize(mkpl_fill_button->GetDefaultWidth()+20, 
			   mkpl_fill_button->GetDefaultHeight());
  mkpl_fill_button->Associate(this);

  mkpl_actions_group->AddFrame(mkpl_comparisons_button);
  mkpl_actions_group->AddFrame(mkpl_tests_button);
  mkpl_actions_group->AddFrame(mkpl_components_button);
  mkpl_actions_group->AddFrame(mkpl_fill_button);
  
  mkpl_actions_group->Resize(mkpl_actions_group->GetDefaultSize());


  // Tab: Graphics
  //==================
  compFrame = mkpl_tab->AddTab("Graphics");
  mkpl_graphics_layout = new TGLayoutHints(kLHintsTop | kLHintsLeft |
					   kLHintsExpandX | kLHintsExpandY,
					   5, 5, 5, 5);
  mkpl_graphics = new TRootEmbeddedCanvas("GraphicsCanvas", compFrame, 200, 200);
  compFrame->AddFrame(mkpl_graphics, mkpl_graphics_layout);
  gROOT->SetStyle("Plain");
  gROOT->GetStyle("Plain")->SetOptStat(Int_t(0000000));
  gROOT->GetStyle("Plain")->SetTitleXSize(0.03);
  gROOT->GetStyle("Plain")->SetTitleYSize(0.025);
  gROOT->GetStyle("Plain")->SetTitleXOffset(1.5);
  gROOT->GetStyle("Plain")->SetTitleYOffset(1.5);
  gROOT->GetStyle("Plain")->SetLabelSize(0.03,"X");
  gROOT->GetStyle("Plain")->SetLabelSize(0.03,"Y");
  //--------------

  this->AddFrame(mkpl_tab, mkpl_tab_layout);
  this->MapSubwindows();
  
  this->SetWindowName("GMAKEPLOTS");
  TGDimension dim(GetSize());
  this->Resize(GetDefaultSize());
  this->Resize(dim);
  this->MapWindow();
}

mkplmakeplots::~mkplmakeplots()
{  
  delete mkpl_info_group_layout;
  delete mkpl_actions_group_layout;
  delete mkpl_info_group;
  delete mkpl_actions_group;
  delete mkpl_info_sim_text;
  delete mkpl_info_exp_text;
  delete mkpl_options;
  delete mkpl_options_layout;
  delete mkpl_comparisons_button;
  delete mkpl_tests_button;
  delete mkpl_components_button;
  delete mkpl_fill_button;
  delete mkpl_info_exp_group;
  delete mkpl_info_sim_group;
  delete mkpl_info_comp_group;
  delete mkpl_info_exp_group_layout;
  delete mkpl_info_sim_group_layout;
  delete mkpl_info_comp_group_layout;
  delete mkpl_info_exp_label;
  delete mkpl_info_sim_label;
  delete mkpl_info_comp_label;
  delete mkpl_info_exp_label_layout;
  delete mkpl_info_sim_label_layout;
  delete mkpl_info_comp_label_layout;
  delete mkpl_info_exp_text_layout;
  delete mkpl_info_sim_text_layout;
  delete mkpl_info_comp_text_layout;

  delete mkpl_file_manager;
  delete mkpl_graphics;
  delete mkpl_graphics_layout;
  delete mkpl_tab;
  delete mkpl_tab_layout;
  delete mkpl_file_menu;
  delete mkpl_file_menu_layout;
  delete mkpl_menu_bar;
  delete mkpl_menu_bar_layout;
}


Bool_t mkplmakeplots::ProcessMessage(Long_t msg, Long_t param1, Long_t param2)
{
  switch (GET_MSG(msg))
    {
    case kC_COMMAND:
      switch (GET_SUBMSG(msg))
	{
	case kCM_MENUSELECT:
	  break;
	case kCM_MENU:
	  switch (param1)
	    {
	    case M_FILE_EXP_OPEN:
	      this->MenuOpenExpFile();
	      break;
	    case M_FILE_EXP_CLOSE:
	      this->MenuCloseExpFile();
	      break;
	    case M_FILE_SIM_OPEN:
	      this->MenuOpenSimFile();
	      break;
	    case M_FILE_SIM_CLOSE:
	      this->MenuCloseSimFile();
	      break;
	    case M_FILE_SAVE_HISTOGRAMS:
	      this->MenuSaveHistograms();
	      break;
	    case M_FILE_SAVE_EPS:
	      this->MenuSaveEPS();
	      break;
	    case M_FILE_EXIT:
	      CloseWindow();
	      break;
	    default:
	      break;
	    }
	  break;
	case kCM_BUTTON:
	  switch (param1)
	    {
	    case 81:
	      {
		if (mkpl_comparison_selector) delete mkpl_comparison_selector;
		this->FillHistograms();
		mkpl_comparison_selector = new mkplcomparisonselector(fClient->GetRoot(), 
								      this, 400, 200, 
								      mkpl_sim_data_manager);
	      }
	      break;
	    default:
	      break;
	    }
	  break;
	default:
	  break;
	}
      break;
    case COMPARISON_SELECTOR_MSG:
      switch (GET_SUBMSG(msg))
	{
	case COMPARISON_SELECTOR_DELETED:
	  mkpl_comparison_selector = 0;
	  break;
	case COMPARISON_SELECTOR_TEST:
	  this->PlotTest(param1);
	  break;
	case COMPARISON_SELECTOR_COMP:
	  this->PlotComparison(param1,param2);
	default:
	  break;
	}
      break;
    default:
      break;
    }
  return kTRUE;
}


void mkplmakeplots::MenuOpenExpFile()
{
  if (mkpl_file_menu->IsEntryEnabled(M_FILE_SAVE_HISTOGRAMS)) 
    mkpl_file_menu->DisableEntry(M_FILE_SAVE_HISTOGRAMS);
  if (mkpl_file_menu->IsEntryEnabled(M_FILE_SAVE_EPS)) 
    mkpl_file_menu->DisableEntry(M_FILE_SAVE_EPS);
  static TString expdir(".");
  TGFileInfo file_info;
  file_info.fFileTypes = filetypes;
  file_info.fIniDir = StrDup(expdir);
  new TGFileDialog(fClient->GetRoot(), this, kFDOpen, &file_info);
  expdir = file_info.fIniDir;
  if (file_info.fFilename)
    {
      this->OpenExpFile(file_info.fFilename);
      this->ExpFileShowInfo();
    }
  return;
}

void mkplmakeplots::MenuOpenSimFile()
{
  if (mkpl_file_menu->IsEntryEnabled(M_FILE_SAVE_HISTOGRAMS)) 
    mkpl_file_menu->DisableEntry(M_FILE_SAVE_HISTOGRAMS);
  if (mkpl_file_menu->IsEntryEnabled(M_FILE_SAVE_EPS)) 
    mkpl_file_menu->DisableEntry(M_FILE_SAVE_EPS);
  static TString simdir(".");
  TGFileInfo file_info;
  file_info.fFileTypes = filetypes;
  file_info.fIniDir = StrDup(simdir);
  new TGFileDialog(fClient->GetRoot(), this, kFDOpen, &file_info);
  simdir = file_info.fIniDir;
  if (file_info.fFilename)
    {
      this->OpenSimFile(file_info.fFilename);
      this->SimFileShowInfo();
    }
  return;
}

void mkplmakeplots::MenuCloseExpFile()
{
  if (mkpl_file_menu->IsEntryEnabled(M_FILE_SAVE_HISTOGRAMS)) 
    mkpl_file_menu->DisableEntry(M_FILE_SAVE_HISTOGRAMS);
  if (mkpl_file_menu->IsEntryEnabled(M_FILE_SAVE_EPS)) 
    mkpl_file_menu->DisableEntry(M_FILE_SAVE_EPS);
  mkpl_graphics->GetCanvas()->Clear();
  mkpl_graphics->GetCanvas()->Modified();
  mkpl_graphics->GetCanvas()->Update();
  mkpl_info_exp_text->Clear();
  gSystem->ProcessEvents();
  if (mkpl_comparison_selector) delete mkpl_comparison_selector;
  mkpl_comparison_selector = 0;
  this->CloseExpFile();
  return;
}

void mkplmakeplots::MenuCloseSimFile()
{
  if (mkpl_file_menu->IsEntryEnabled(M_FILE_SAVE_HISTOGRAMS)) 
    mkpl_file_menu->DisableEntry(M_FILE_SAVE_HISTOGRAMS);
  if (mkpl_file_menu->IsEntryEnabled(M_FILE_SAVE_EPS)) 
    mkpl_file_menu->DisableEntry(M_FILE_SAVE_EPS);
  mkpl_graphics->GetCanvas()->Clear();
  mkpl_graphics->GetCanvas()->Modified();
  mkpl_graphics->GetCanvas()->Update();
  mkpl_info_sim_text->Clear();
  gSystem->ProcessEvents();
  if (mkpl_comparison_selector) delete mkpl_comparison_selector;
  mkpl_comparison_selector = 0;
  this->CloseSimFile();
  return;
}


void mkplmakeplots::ExpFileShowInfo()
{
  std::ostringstream summary;
  summary << mkpl_file_manager->GetExpFile()->GetBaseName() 
	  << "\n\n";
  mkpl_exp_data_manager->GetSummary(summary);
  string s(summary.str());
  mkpl_info_exp_text->LoadBuffer(s.c_str());;
  return;
}

void mkplmakeplots::SimFileShowInfo()
{
  std::ostringstream summary;
  summary << mkpl_file_manager->GetSimFile()->GetBaseName() 
	  << "\n\n";
  mkpl_sim_data_manager->GetSummary(summary);
  string s(summary.str());
  mkpl_info_sim_text->LoadBuffer(s.c_str());;
  return;
}


void mkplmakeplots::FillHistograms()
{
  if (mkpl_comparison_selector) delete mkpl_comparison_selector;
  mkpl_comparison_selector = 0;
  int retval;
  if (!mkpl_exp_data_manager)
    {
      new TGMsgBox(fClient->GetRoot(), this, "Warning","You should open a Experimental File",
		   kMBIconStop, kMBOk, &retval);
    }
  else if (!mkpl_sim_data_manager)
    {
      new TGMsgBox(fClient->GetRoot(), this, "Warning","You should open a Simulation File",
		   kMBIconStop, kMBOk, &retval);

    }
  else if (mkpl_sim_data_manager->GetTargetZ() != mkpl_exp_data_manager->GetTargetZ())
    {
      new TGMsgBox(fClient->GetRoot(), this, "Warning",
		   "Experimental and Simulation are not comparable",
		   kMBIconStop, kMBOk, &retval);
    }
  else
    {
      this->HideFrame(mkpl_tab);
      if (mkpl_tests_button->GetState() == kButtonDown)
	mkpl_sim_data_manager->PrepareTestHistograms();
      if (mkpl_comparisons_button->GetState() == kButtonDown)
	mkpl_sim_data_manager->PrepareComparisonHistograms(mkpl_exp_data_manager);
      mkpl_sim_data_manager->FillHistograms();
      this->ShowFrame(mkpl_tab);
    }
  mkpl_file_menu->EnableEntry(M_FILE_SAVE_HISTOGRAMS);
  mkpl_file_menu->EnableEntry(M_FILE_SAVE_EPS);
  return;
}

void mkplmakeplots::MenuSaveHistograms()
{
  static TString histosdir(".");
  TGFileInfo file_info;
  file_info.fFileTypes = filetypes;
  file_info.fIniDir = StrDup(histosdir);
  //  TString fname(mkpl_file_manager->GetSimFile()->GetBaseName().c_str());
  //  fname += "_histograms.root";
  new TGFileDialog(fClient->GetRoot(), this, kFDSave, &file_info);
  histosdir = file_info.fIniDir;
  if (file_info.fFilename)
    {

      TFile * hf = new TFile(file_info.fFilename,"RECREATE");

      vector<mkplcomparisonhistograms*> * ch = mkpl_sim_data_manager->GetComparisonHistograms();
      for (vector<mkplcomparisonhistograms*>::iterator it = ch->begin(); it != ch->end(); it++)
	{
	  (*it)->SaveHistograms(hf);
	}
     
      hf->Write();
      delete hf;
    }
  return;
}


void mkplmakeplots::MenuSaveEPS()
{
  TString fname(mkpl_file_manager->GetSimFile()->GetBaseName().c_str());
  bool comp(mkpl_components_button->GetState() == kButtonDown);
  //  Color_t canvas_color = gROOT->GetStyle("Plain")->GetCanvasColor();
  //  Color_t frame_color = gROOT->GetStyle("Plain")->GetFrameFillColor();
  Color_t canvas_color = gPad->GetCanvas()->GetFillColor();
  Color_t frame_color = gPad->GetFrame()->GetFillColor();
  gPad->GetCanvas()->SetFillColor(10);
  gPad->GetFrame()->SetFillColor(10);
  gPad->Modified();
  vector<mkplcomparisonhistograms*> * ch = mkpl_sim_data_manager->GetComparisonHistograms();
  for (vector<mkplcomparisonhistograms*>::iterator it = ch->begin(); it != ch->end(); it++)
    {
      (*it)->SaveEPS(fname,comp);
    }
  gPad->GetCanvas()->SetFillColor(canvas_color);
  gPad->GetFrame()->SetFillColor(frame_color);
  gPad->Modified();
  return;
}


void mkplmakeplots::PlotComparison(const Long_t p1, const Long_t p2)
{
  mkplcomparisonhistograms * ch = mkpl_sim_data_manager->GetComparisonHistograms()->operator[](p1);
  bool comp(mkpl_components_button->GetState() == kButtonDown);
  mkpl_graphics->GetCanvas()->Clear();
  if (p2 < 100)
    {
      if (p2 == 1) 
	{
	  mkpl_info_comp_text->Clear();
	  vector<mkplcomparisonagreement> agr = ch->PlotDE(comp);
	  for (vector<mkplcomparisonagreement>::iterator ds = agr.begin();
	       ds != agr.end(); ds++)
	    {
	      TString text("F = ");
	      text += ds->Fbar;
	      text += "\n<F> = ";
	      text += ds->Favg;
	      text += "\nFmin = ";
	      text += ds->Fmin;
	      text += "\nFmax = ";
	      text += ds->Fmax;
	      text += "\nNs = ";
	      text += ds->Ns;
	      text += "\n*************\n";
	      mkpl_info_comp_text->LoadBuffer(text);
	    }
	}
      else if (p2 == 2) 
	{
	  mkpl_info_comp_text->Clear();
	  vector<mkplcomparisonagreement> agr = ch->PlotDEcm(comp);
	  for (vector<mkplcomparisonagreement>::iterator ds = agr.begin();
	       ds != agr.end(); ds++)
	    {
	      TString text("F = ");
	      text += ds->Fbar;
	      text += "\n<F> = ";
	      text += ds->Favg;
	      text += "\nFmin = ";
	      text += ds->Fmin;
	      text += "\nFmax = ";
	      text += ds->Fmax;
	      text += "\nNs = ";
	      text += ds->Ns;
	      text += "\n*************\n";
	      mkpl_info_comp_text->LoadBuffer(text);
	    }
	}
      else if (p2 == 3) 
	{
	  mkpl_info_comp_text->Clear();
	  vector<mkplcomparisonagreement> agr = ch->PlotDA(comp);
	  for (vector<mkplcomparisonagreement>::iterator ds = agr.begin();
	       ds != agr.end(); ds++)
	    {
	      TString text("F = ");
	      text += ds->Fbar;
	      text += "\n<F> = ";
	      text += ds->Favg;
	      text += "\nFmin = ";
	      text += ds->Fmin;
	      text += "\nFmax = ";
	      text += ds->Fmax;
	      text += "\nNs = ";
	      text += ds->Ns;
	      text += "\n*************\n";
	      mkpl_info_comp_text->LoadBuffer(text);
	    }
	}
      else if (p2 == 4) 
	{
	  mkpl_info_comp_text->Clear();
	  vector<mkplcomparisonagreement> agr = ch->PlotDAcm(comp);
	  for (vector<mkplcomparisonagreement>::iterator ds = agr.begin();
	       ds != agr.end(); ds++)
	    {
	      TString text("F = ");
	      text += ds->Fbar;
	      text += "\n<F> = ";
	      text += ds->Favg;
	      text += "\nFmin = ";
	      text += ds->Fmin;
	      text += "\nFmax = ";
	      text += ds->Fmax;
	      text += "\nNs = ";
	      text += ds->Ns;
	      text += "\n*************\n";
	      mkpl_info_comp_text->LoadBuffer(text);
	    }
	}
    }
  else 
    {
      int c = p2/100;
      int w = p2 - 100*c - 1;
      if (c == 1) 
	{
	  mkpl_info_comp_text->Clear();
	  if (w == ch->GetNumOfAngles())
	    {
	      mkpl_graphics->GetCanvas()->SetLogy(1);
	    }
	  vector<mkplcomparisonagreement> agr = ch->PlotDD(w,comp); 
	  for (vector<mkplcomparisonagreement>::iterator ds = agr.begin();
	       ds != agr.end(); ds++)
	    {
	      TString text("F = ");
	      text += ds->Fbar;
	      text += "\n<F> = ";
	      text += ds->Favg;
	      text += "\nFmin = ";
	      text += ds->Fmin;
	      text += "\nFmax = ";
	      text += ds->Fmax;
	      text += "\nNs = ";
	      text += ds->Ns;
	      text += "\n*************\n";
	      mkpl_info_comp_text->LoadBuffer(text);
	    }
	}
      if (c == 2) 
	{
	  mkpl_info_comp_text->Clear();
	  vector<mkplcomparisonagreement> agr = ch->PlotDDcm(w,comp);
	  if (w == ch->GetNumOfAnglesCM())
	    {
	      mkpl_graphics->GetCanvas()->SetLogy(1);
	    }
	  for (vector<mkplcomparisonagreement>::iterator ds = agr.begin();
	       ds != agr.end(); ds++)
	    {
	      TString text("F = ");
	      text += ds->Fbar;
	      text += "\n<F> = ";
	      text += ds->Favg;
	      text += "\nFmin = ";
	      text += ds->Fmin;
	      text += "\nFmax = ";
	      text += ds->Fmax;
	      text += "\nNs = ";
	      text += ds->Ns;
	      text += "\n*************\n";
	      mkpl_info_comp_text->LoadBuffer(text);
	    }
	}
      if (c == 3) 
	{
	  mkpl_info_comp_text->Clear();
	  vector<mkplcomparisonagreement> agr = ch->PlotDDA(w,comp);
	  for (vector<mkplcomparisonagreement>::iterator ds = agr.begin();
	       ds != agr.end(); ds++)
	    {
	      TString text("F = ");
	      text += ds->Fbar;
	      text += "\n<F> = ";
	      text += ds->Favg;
	      text += "\nFmin = ";
	      text += ds->Fmin;
	      text += "\nFmax = ";
	      text += ds->Fmax;
	      text += "\nNs = ";
	      text += ds->Ns;
	      text += "\n*************\n";
	      mkpl_info_comp_text->LoadBuffer(text);
	    }
	}
      if (c == 4) 
	{
	  mkpl_info_comp_text->Clear();
	  vector<mkplcomparisonagreement> agr = ch->PlotDDAcm(w,comp);
	  for (vector<mkplcomparisonagreement>::iterator ds = agr.begin();
	       ds != agr.end(); ds++)
	    {
	      TString text("F = ");
	      text += ds->Fbar;
	      text += "\n<F> = ";
	      text += ds->Favg;
	      text += "\nFmin = ";
	      text += ds->Fmin;
	      text += "\nFmax = ";
	      text += ds->Fmax;
	      text += "\nNs = ";
	      text += ds->Ns;
	      text += "\n*************\n";
	      mkpl_info_comp_text->LoadBuffer(text);
	    }
	}

    }
  mkpl_graphics->GetCanvas()->Modified();
  mkpl_graphics->GetCanvas()->Update();
  gSystem->ProcessEvents(); 
  return;
} 


void mkplmakeplots::PlotTest(const Long_t p1)
{
  mkpltesthistograms * tst = mkpl_sim_data_manager->GetTestHistograms();
  mkpl_graphics->GetCanvas()->Clear();
  switch (p1)
    {
    case 1:
      tst->PlotBaryonNumber();
      break;
    case 2:
      tst->PlotCharge();
      break;
    case 3:
      tst->PlotEnergy();
      break;
    case 4:
      tst->PlotMomentum();
      break;
    case 5:
      tst->PlotPx();
      break;
    case 6:
      tst->PlotPy();
      break;
    case 7:
      tst->PlotPz();
      break;
    case 8:
      tst->PlotTheta();
      break;
    case 9:
      tst->PlotThetaPreeq();
      break;
    case 10:
      tst->PlotThetaEvap();
      break;
    case 11:
      tst->PlotPhi();
      break;
    case 12:
      tst->PlotPhiPreeq();
      break;
    case 13:
      tst->PlotPhiEvap();
      break;
    case 14:
      tst->PlotPhiNucleons();
      break;
    case 15:
      tst->PlotTypeOfFragments();
      break;
    default:
      break;
    }
  mkpl_graphics->GetCanvas()->Modified();
  mkpl_graphics->GetCanvas()->Update();
  gSystem->ProcessEvents();
  return;
}
