#include "mkplmainframe.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TMultiGraph.h"



#include "mkplmessages.h"


const char * mkplmainframe::filetypes[] = 
  {
    "Root files", "*.root",
    "All files", "*",
    0 , 0
  };


mkplmainframe::mkplmainframe(const TGWindow *p, UInt_t w, UInt_t h)
  : TGMainFrame(p,w,h), mkpl_exp_data_manager(0), mkpl_energy_selector(0), 
    mkpl_angle_selector(0), mkpl_range_selector(0)
{
  mkpl_file_manager = new mkplfilemanager();
  
  
  // File Menu
  mkpl_file_menu_layout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
  mkpl_file_menu = new TGPopupMenu(fClient->GetRoot());
  mkpl_file_menu->AddEntry("&Open", M_FILE_OPEN);
  mkpl_file_menu->AddEntry("&Close", M_FILE_CLOSE);
  mkpl_file_menu->AddSeparator();
  mkpl_file_menu->AddEntry("E&xit", M_FILE_EXIT);
  mkpl_file_menu->Associate(this);


  // Menu Bar
  mkpl_menu_bar_layout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 1, 1);
  mkpl_menu_bar = new TGMenuBar(this, 1, 1, kHorizontalFrame);
  mkpl_menu_bar->AddPopup("&File", mkpl_file_menu, mkpl_file_menu_layout);

  this->AddFrame(mkpl_menu_bar,mkpl_menu_bar_layout);
  this->MapSubwindows();

  // The TAB 
  mkpl_tab = new TGTab(this, 300, 300);
  mkpl_tab_layout =  new TGLayoutHints(kLHintsBottom | kLHintsExpandX | kLHintsExpandY, 2, 2, 5, 1);
  
  // Tab: The Summary
  TGCompositeFrame * compFrame = mkpl_tab->AddTab("Summary");
  mkpl_summary_layout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,
				       5, 5, 5, 5);
  mkpl_summary_text_viewer = new TGTextView(compFrame, 100, 100);
  compFrame->AddFrame(mkpl_summary_text_viewer, mkpl_summary_layout);

  // Tab: The Details
  compFrame = mkpl_tab->AddTab("Details");
  mkpl_details_layout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,
				       5, 5, 5, 5);
  mkpl_details_text_viewer = new TGTextView(compFrame, 100, 100);
  compFrame->AddFrame(mkpl_details_text_viewer, mkpl_details_layout);

  // Tab: The Graphics
  compFrame = mkpl_tab->AddTab("Graphics");
  mkpl_graphics_layout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,
					5, 5, 5, 5);
  mkpl_graphics_canvas = new TRootEmbeddedCanvas("GraphicsCanvas", compFrame, 200, 200);
  compFrame->AddFrame(mkpl_graphics_canvas, mkpl_graphics_layout);


  this->AddFrame(mkpl_tab, mkpl_tab_layout);
  this->MapSubwindows();

  this->SetWindowName("GSHOWDATA");
  TGDimension dim(GetSize());
  this->Resize(GetDefaultSize());
  this->Resize(dim);
  this->MapWindow();
}


mkplmainframe::~mkplmainframe()
{
  if (mkpl_exp_data_manager) delete mkpl_exp_data_manager;
  if (mkpl_file_manager) delete mkpl_file_manager;
  if (mkpl_energy_selector) delete mkpl_energy_selector;
  if (mkpl_angle_selector) delete mkpl_angle_selector;
  if (mkpl_range_selector) delete mkpl_range_selector;

  delete mkpl_graphics_canvas;
  delete mkpl_graphics_layout;
  delete mkpl_details_text_viewer;
  delete mkpl_details_layout;
  delete mkpl_summary_text_viewer;
  delete mkpl_summary_layout;
  delete mkpl_tab;
  delete mkpl_tab_layout;
  delete mkpl_file_menu;
  delete mkpl_file_menu_layout;
  delete mkpl_menu_bar;
  delete mkpl_menu_bar_layout;
}


Bool_t mkplmainframe::ProcessMessage(Long_t msg, Long_t param1, Long_t param2)
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
	    case M_FILE_OPEN:
	      this->MenuOpenFile();
	      break;
	    case M_FILE_CLOSE:
	      this->MenuCloseFile();
	      break;
	    case M_FILE_EXIT:
	      CloseWindow();
	      break;
	    default:
	      break;
	    }
	  break;
	default:
	  break;
	}
      break;
    case ENERGY_SELECTOR_MSG:
      switch (GET_SUBMSG(msg))
	{
	case ENERGY_SELECTOR_DELETED:
	  mkpl_energy_selector = 0;
	  this->MenuCloseFile();
	  break;
	case ENERGY_SELECTOR_DE:
	  this->ShowDEdata(param1);
	  break;
	case ENERGY_SELECTOR_DA:
	  this->ShowDAdata(param1);
	  break;
	case ENERGY_SELECTOR_DD:
	  this->ShowDDangles(param1);
	  break;
	case ENERGY_SELECTOR_DDA:
	  this->ShowDDAranges(param1);
	  break;
	default:
	  break;
	}
      break;
    case RANGE_SELECTOR_MSG:
      switch (GET_SUBMSG(msg))
	{
	case RANGE_SELECTOR_DELETED:
	  mkpl_range_selector = 0;
	  break;
	case RANGE_SELECTOR_ONE:
	  this->ShowDDAdata(param1,param2);
	  break;
	case RANGE_SELECTOR_CLOSE:
	  delete mkpl_range_selector;
	  mkpl_range_selector = 0;
	  break;
	default:
	  break;
	}
      break;
    case ANGLE_SELECTOR_MSG:
      switch (GET_SUBMSG(msg))
      {
      case ANGLE_SELECTOR_DELETED:
	mkpl_angle_selector = 0;
	break;
      case ANGLE_SELECTOR_ONE:
	this->ShowDDdata(param1,param2);
	break;
      case ANGLE_SELECTOR_ALL:
	this->ShowDDdataALL(param1);
	break;
      case ANGLE_SELECTOR_CLOSE:
	delete mkpl_angle_selector;
	mkpl_angle_selector = 0;
	break;
      default:
	break;
      }
      break;
    default:
      break;
    }
  return kTRUE;
}


void mkplmainframe::MenuOpenFile()
{
  static TString dir(".");
  TGFileInfo file_info;
  file_info.fFileTypes = filetypes;
  file_info.fIniDir = StrDup(dir);
  new TGFileDialog(fClient->GetRoot(), this, kFDOpen, &file_info);
  dir = file_info.fIniDir;
  if (file_info.fFilename)
    {
      this->OpenExpFile(file_info.fFilename);
      std::ostringstream summary;
      mkpl_exp_data_manager->GetSummary(summary);
      string s(summary.str());
      mkpl_summary_text_viewer->LoadBuffer(s.c_str());
      if (mkpl_range_selector) delete mkpl_range_selector;
      if (mkpl_angle_selector) delete mkpl_angle_selector;
      if (mkpl_energy_selector) delete mkpl_energy_selector;
      mkpl_energy_selector = new mkplenergyselector(fClient->GetRoot(), this, 400, 200, 
						    mkpl_exp_data_manager);
    }
  return;
}


void mkplmainframe::MenuCloseFile()
{
  mkpl_summary_text_viewer->Clear();
  mkpl_details_text_viewer->Clear();
  mkpl_graphics_canvas->GetCanvas()->Clear();
  mkpl_graphics_canvas->GetCanvas()->Modified();
  mkpl_graphics_canvas->GetCanvas()->Update();
  gSystem->ProcessEvents();
  if (mkpl_range_selector)
    {
      delete mkpl_range_selector;
      mkpl_range_selector = 0;
    }
  if (mkpl_angle_selector) 
    {
      delete mkpl_angle_selector;
      mkpl_angle_selector = 0;
    }
  if (mkpl_energy_selector) 
    {
      delete mkpl_energy_selector;
      mkpl_energy_selector = 0;
    }
  this->CloseExpFile();
  return;
}


void mkplmainframe::ShowDEdata(const int n)
{
  gROOT->SetStyle("Plain");
  mkpl_graphics_canvas->GetCanvas()->Clear();
  mkpl_graphics_canvas->GetCanvas()->SetGrid(1,1);
  TGraphErrors * graph = (TGraphErrors*)(mkpl_exp_data_manager->GetDE(n));
  if (graph)
    {
      std::ostringstream details;
      mkpl_exp_data_manager->GetDEdetails(n,details);
      string det(details.str());
      mkpl_details_text_viewer->LoadBuffer(det.c_str());
      graph->SetMarkerColor(2);
      graph->SetMarkerStyle(20);
      graph->SetMarkerSize(0.75);
      graph->Draw("ALP");
      graph->GetXaxis()->SetTitle("T (MeV)");
      graph->GetYaxis()->SetTitle("#frac{d#sigma}{dT} (mb/MeV)");
      graph->GetYaxis()->SetTitleOffset(1.1);
      graph->Draw("ALP");
      mkpl_graphics_canvas->GetCanvas()->Modified();
      mkpl_graphics_canvas->GetCanvas()->Update();
      gSystem->ProcessEvents();
    }
  return;
}

void mkplmainframe::ShowDAdata(const int n)
{
  gROOT->SetStyle("Plain");
  mkpl_graphics_canvas->GetCanvas()->Clear();
  mkpl_graphics_canvas->GetCanvas()->SetGrid(1,1);
  TGraphErrors * graph = (TGraphErrors*)(mkpl_exp_data_manager->GetDA(n));
  if (graph)
    {
      std::ostringstream details;
      mkpl_exp_data_manager->GetDAdetails(n,details);
      string det(details.str());
      mkpl_details_text_viewer->LoadBuffer(det.c_str());
      graph->SetMarkerColor(2);
      graph->SetMarkerStyle(20);
      graph->SetMarkerSize(0.75);
      graph->Draw("ALP");
      graph->GetXaxis()->SetTitle("#theta (deg)");
      graph->GetYaxis()->SetTitle("#frac{d#sigma}{d#theta} (mb/deg)");
      graph->GetYaxis()->SetTitleOffset(1.1);
      graph->Draw("ALP");
      mkpl_graphics_canvas->GetCanvas()->Modified();
      mkpl_graphics_canvas->GetCanvas()->Update();
      gSystem->ProcessEvents();
    }
  return;
}

void mkplmainframe::ShowDDangles(const int n)
{
  if (mkpl_angle_selector) delete mkpl_angle_selector;
  mkpl_angle_selector = new mkplangleselector(fClient->GetRoot(), this, 450, 220, 
					      mkpl_exp_data_manager,n);
  
  return;
}



void mkplmainframe::ShowDDdata(const int energy, const int angle)
{
  gROOT->SetStyle("Plain");
  mkpl_graphics_canvas->GetCanvas()->Clear();
  mkpl_graphics_canvas->GetCanvas()->SetGrid(1,1);
  TGraphErrors * graph = (TGraphErrors*)(mkpl_exp_data_manager->GetDD(energy,angle));
  if (graph)
    {
      std::ostringstream details;
      mkpl_exp_data_manager->GetDDdetails(energy,angle,details);
      string det(details.str());
      mkpl_details_text_viewer->LoadBuffer(det.c_str());
      graph->SetMarkerColor(2);
      graph->SetMarkerStyle(20);
      graph->SetMarkerSize(0.75);
      graph->Draw("ALP");
      graph->GetXaxis()->SetTitle("T (MeV)");
      graph->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{dT d#Omega} (mb/deg/Mev)");
      graph->GetYaxis()->SetTitleOffset(1.1);
      graph->Draw("ALP");
      mkpl_graphics_canvas->GetCanvas()->Modified();
      mkpl_graphics_canvas->GetCanvas()->Update();
      gSystem->ProcessEvents();
    }
  return;
}


void mkplmainframe::ShowDDdataALL(const int energy)
{
  gROOT->SetStyle("Plain");
  mkpl_graphics_canvas->GetCanvas()->Clear();
  mkpl_graphics_canvas->GetCanvas()->SetGrid(1,1);

  std::ostringstream details;
  mkpl_exp_data_manager->GetDDdetails(energy,details);
  string det(details.str());
  mkpl_details_text_viewer->LoadBuffer(det.c_str());

  TMultiGraph * mg = mkpl_exp_data_manager->GetDD(energy);
  if (mg)
  {
    mg->Draw("ALP");
    mg->GetXaxis()->SetTitle("T (MeV)");
    mg->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{dT d#Omega} (mb/deg/Mev)");
    mg->GetYaxis()->SetTitleOffset(1.1);
    mg->Draw("ALP");
  }
  mkpl_graphics_canvas->GetCanvas()->Modified();
  mkpl_graphics_canvas->GetCanvas()->Update();
  gSystem->ProcessEvents();
  
  return;
}

void mkplmainframe::ShowDDAranges(const int n)
{
  if (mkpl_range_selector) delete mkpl_range_selector;
  mkpl_range_selector = new mkplrangeselector(fClient->GetRoot(), this, 450, 220, 
					      mkpl_exp_data_manager,n);
  
  return;
}

void mkplmainframe::ShowDDAdata(const int energy, const int range)
{
  gROOT->SetStyle("Plain");
  mkpl_graphics_canvas->GetCanvas()->Clear();
  mkpl_graphics_canvas->GetCanvas()->SetGrid(1,1);
  TGraphErrors * graph = (TGraphErrors*)(mkpl_exp_data_manager->GetDDA(energy,range));
  if (graph)
    {
      std::ostringstream details;
      mkpl_exp_data_manager->GetDDAdetails(energy,range,details);
      string det(details.str());
      mkpl_details_text_viewer->LoadBuffer(det.c_str());
      graph->SetMarkerColor(2);
      graph->SetMarkerStyle(20);
      graph->SetMarkerSize(0.75);
      graph->Draw("ALP");
      graph->GetXaxis()->SetTitle("#theta (degrees)");
      graph->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} (mb/deg)");
      graph->GetYaxis()->SetTitleOffset(1.1);
      graph->Draw("ALP");
      mkpl_graphics_canvas->GetCanvas()->Modified();
      mkpl_graphics_canvas->GetCanvas()->Update();
      gSystem->ProcessEvents();
    }
  return;
}
