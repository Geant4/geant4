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
// $Id: G4PlotManager.cc 85314 2014-10-27 14:33:08Z ihrivnac $

// Author: Ivana Hrivnacova, 02/06/2015  (ivana@ipno.in2p3.fr)

#include "G4HnInformation.hh"
#include "G4PlotManager.hh"
#include "G4ios.hh"

#if defined(TOOLS_USE_FREETYPE)
#include <tools/sg/text_freetype>
#include <tools/xml/xml_style>
#include <tools/xml/wrap_viewplot_style> 
  //inlib/xml/viewplot.style file embeded in an inline function.

namespace {

// from g4tools/test/viewplot.cpp  
//_____________________________________________________________________________
void HD_style(tools::sg::plots& a_plots,float a_line_width) {
  std::vector<tools::sg::plotter*> plotters;
  a_plots.plotters(plotters);
  tools_vforcit(tools::sg::plotter*,plotters,it) {
    tools::sg::plotter* _plotter = *it;
    _plotter->bins_style(0).line_width = a_line_width;
    _plotter->inner_frame_style().line_width = a_line_width;
    _plotter->grid_style().line_width = a_line_width;
    _plotter->x_axis().line_style().width = a_line_width;
    _plotter->x_axis().ticks_style().width = a_line_width;
    _plotter->y_axis().line_style().width = a_line_width;
    _plotter->y_axis().ticks_style().width = a_line_width;
    _plotter->z_axis().line_style().width = a_line_width;
    _plotter->z_axis().ticks_style().width = a_line_width;

    // needed if font is hershey :
    _plotter->title_style().line_width = a_line_width;
    _plotter->infos_style().line_width = a_line_width;
    _plotter->title_box_style().line_width = a_line_width;

    _plotter->x_axis().labels_style().line_width = a_line_width;
    _plotter->x_axis().mag_style().line_width = a_line_width;
    _plotter->x_axis().title_style().line_width = a_line_width;

    _plotter->y_axis().labels_style().line_width = a_line_width;
    _plotter->y_axis().mag_style().line_width = a_line_width;
    _plotter->y_axis().title_style().line_width = a_line_width;

    _plotter->z_axis().labels_style().line_width = a_line_width;
    _plotter->z_axis().mag_style().line_width = a_line_width;
    _plotter->z_axis().title_style().line_width = a_line_width;
  }    
}

// from g4tools/test/viewplot.cpp  
//_____________________________________________________________________________
void regions_style(tools::sg::plots& a_plots,float a_plotter_scale = 1) {
  // Rescale some plotter parameters (for example margins) according to the number of regions.
  // We assume that these parameters had been set previously according to one plot per page.
  // Then this function must be applied after all the styles had been applied (because
  // a plotting style may set these parameters).

  float ww_wc = a_plots.width;
  float wh_wc = a_plots.height;
  float rw_wc = ww_wc/a_plots.cols;
  float rh_wc = wh_wc/a_plots.rows;

  float cooking = 1.2f; //if increased the data area is diminished.

  float wfac = (rw_wc/ww_wc)*cooking;
  float hfac = (rh_wc/wh_wc)*cooking;

  float label_cooking = 1.6f; //if increased the labels are bigger.

  if((a_plots.cols.value()>=4)&&(a_plots.cols.value()>a_plots.rows.value())) label_cooking = 0.9f;

  float title_cooking = 1.1f; //extra title cooking.

  a_plots.plotter_scale = a_plotter_scale;

  std::vector<tools::sg::plotter*> plotters;
  a_plots.plotters(plotters);
  tools_vforcit(tools::sg::plotter*,plotters,it) {
    tools::sg::plotter* _plotter = *it;

    _plotter->left_margin = _plotter->left_margin * wfac;
    _plotter->right_margin = _plotter->right_margin * wfac;
    _plotter->bottom_margin = _plotter->bottom_margin * hfac;
    _plotter->top_margin = _plotter->top_margin * hfac;

    _plotter->x_axis().tick_length = _plotter->x_axis().tick_length * wfac;
    _plotter->y_axis().tick_length = _plotter->y_axis().tick_length * hfac;

    _plotter->title_to_axis = _plotter->title_to_axis * hfac;
    _plotter->title_height = _plotter->title_height * hfac * title_cooking;

    _plotter->x_axis().label_height = _plotter->x_axis().label_height * hfac * label_cooking;
    _plotter->y_axis().label_height = _plotter->y_axis().label_height * hfac * label_cooking;

  }    
}

// from g4tools/test/viewplot.cpp  
//_____________________________________________________________________________
bool load_embeded_styles(tools::xml::styles& a_styles) {
  std::string ss;
  unsigned int linen;
  const char** lines = viewplot_style(linen);
  for(unsigned int index=0;index<linen;index++) {
    std::string s = lines[index];
    tools::replace(s,"@@double_quote@@","\"");
    tools::replace(s,"@@back_slash@@","\\");
    ss += s + "\n";
  }
  return tools::xml::load_style_string(a_styles,ss);
}

}
#endif

//
// static data
//

//_____________________________________________________________________________
G4PlotParameters G4PlotManager::fgPlotParameters;

//
// ctors, dtor
//

//_____________________________________________________________________________
G4PlotManager::G4PlotManager(const G4AnalysisManagerState& state)
 : fState(state),
   fViewer(nullptr),
   fFileName()
{
#if defined(TOOLS_USE_FREETYPE)
  //////////////////////////////////////////////////////////////////////////////
  /// plotting, high resolution with freetype fonts and by using styles : //////
  //////////////////////////////////////////////////////////////////////////////
#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    G4cout << "... using high resolution with Freetype fonts" << G4endl;
#endif
  //Have vertical A4 :
  // unsigned int ww = 2000; //to have better antialising on freetype fonts.
  // float A4 = 29.7f/21.0f;
  // unsigned int wh = (unsigned int)(float(ww)*A4*0.80);
  static tools::sg::text_freetype ttf;
  fViewer.reset(new tools::viewplot(G4cout, ttf,
                                    fgPlotParameters.GetColumns(),
                                    fgPlotParameters.GetRows(), 
                                    fgPlotParameters.GetWidth(), 
                                    fgPlotParameters.GetHeight()));
  fViewer->plots().view_border = false;
  load_embeded_styles(fViewer->styles());
  fViewer->styles().add_colormap("default",tools::sg::style_default_colormap());
  fViewer->styles().add_colormap("ROOT",tools::sg::style_ROOT_colormap());
#else
  // cretae a viewer with default parameters
#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    G4cout << "... using low resolution with Hershey fonts" << G4endl;
#endif
  fViewer.reset(new tools::viewplot(G4cout, 
                                    fgPlotParameters.GetColumns(), 
                                    fgPlotParameters.GetRows(), 
                                    fgPlotParameters.GetWidth(), 
                                    fgPlotParameters.GetHeight()));
  fViewer->plots().view_border = false;
#endif
}

//_____________________________________________________________________________
G4PlotManager::~G4PlotManager()
{}

// 
// private methods
//

//_____________________________________________________________________________
G4bool G4PlotManager::WritePage()
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("write a page in", "plot file", fFileName);
#endif

#if defined(TOOLS_USE_FREETYPE)
  HD_style(fViewer->plots(), 5);
  regions_style(fViewer->plots(), fgPlotParameters.GetScale());
#endif

  G4bool result = fViewer->write_page();
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot write a page in the plot file " << fFileName;
    G4Exception("G4PlotManager::WritePage()",
                "Analysis_W022", JustWarning, description);
  }

  // clear viewers plots
  fViewer->plots().init_sg(); 
    //it will recreate the sg::plotters and then reset the styles on new ones.

#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) 
    fState.GetVerboseL3()->Message("write a page in", "plot file", fFileName);
#endif

  return result;
}  

// 
// public methods
//

//_____________________________________________________________________________
G4bool G4PlotManager::OpenFile(const G4String& fileName)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("open", "plot file", fileName);
#endif

  // Keep filename for logging
  fFileName = fileName;

  G4bool result = fViewer->open_file(fileName);
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot open plot file " << fileName;
    G4Exception("G4PlotManager::OpenFile()",
                "Analysis_W001", JustWarning, description);
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()->Message("open", "plot file", fileName);
#endif

  return result;
}

//_____________________________________________________________________________
G4bool G4PlotManager::CloseFile()
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("close", "plot file", fFileName);
#endif

  G4bool result = fViewer->close_file();
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot close the plot file.";
    G4Exception("G4PlotManager::CloseFile()",
                "Analysis_W021", JustWarning, description);
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()->Message("close", "plot file", fFileName);
#endif

  return result;
}
