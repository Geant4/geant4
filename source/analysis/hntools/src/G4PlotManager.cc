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

// Author: Ivana Hrivnacova, 02/06/2015  (ivana@ipno.in2p3.fr)

#include "G4HnInformation.hh"
#include "G4PlotManager.hh"
#include "G4AnalysisUtilities.hh"
#include "G4ios.hh"

#if defined(TOOLS_USE_FREETYPE)
#include "toolx/sg/text_freetype"
#include "toolx/xml/xml_style"
#include "tools/xml/wrap_viewplot_fonts_google_style"
  //inlib/xml/viewplot.style file embeded in an inline function.
#include "tools/font/lato_regular_ttf"
#include "tools/font/roboto_bold_ttf"

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
  const char** lines = viewplot_fonts_google_style(linen);
  for(unsigned int index=0;index<linen;index++) {
    std::string s = lines[index];
    tools::replace(s,"@@double_quote@@","\"");
    tools::replace(s,"@@back_slash@@","\\");
    ss += s + "\n";
  }
  return toolx::xml::load_style_string(a_styles,ss);
}

}
#endif

using namespace G4Analysis;

//
// ctors, dtor
//

//_____________________________________________________________________________
G4PlotManager::G4PlotManager(const G4AnalysisManagerState& state)
 : fState(state)
{
#if defined(TOOLS_USE_FREETYPE)
  //////////////////////////////////////////////////////////////////////////////
  /// plotting, high resolution with freetype fonts and by using styles : //////
  //////////////////////////////////////////////////////////////////////////////
  fState.Message(kVL1,  "... using high resolution with Freetype fonts", "");
  //Have vertical A4 :
  // unsigned int ww = 2000; //to have better antialising on freetype fonts.
  // float A4 = 29.7f/21.0f;
  // unsigned int wh = (unsigned int)(float(ww)*A4*0.80);
  static toolx::sg::text_freetype ttf;
  ttf.add_embedded_font(tools::sg::font_lato_regular_ttf(),tools::font::lato_regular_ttf);
  ttf.add_embedded_font(tools::sg::font_roboto_bold_ttf(),tools::font::roboto_bold_ttf);
  fViewer = std::make_unique<tools::viewplot>(G4cout, ttf,
                                    fPlotParameters.GetColumns(),
                                    fPlotParameters.GetRows(),
                                    fPlotParameters.GetWidth(),
                                    fPlotParameters.GetHeight());
  fViewer->plots().view_border = false;
  load_embeded_styles(fViewer->styles());
  fViewer->styles().add_colormap("default",tools::sg::style_default_colormap());
  fViewer->styles().add_colormap("ROOT",tools::sg::style_ROOT_colormap());
#else
  // cretae a viewer with default parameters
  fState.Message(kVL1,  "... using low resolution with Hershey fonts", "");
  fViewer = std::make_unique<tools::viewplot>(G4cout,
                                    fPlotParameters.GetColumns(),
                                    fPlotParameters.GetRows(),
                                    fPlotParameters.GetWidth(),
                                    fPlotParameters.GetHeight());
  fViewer->plots().view_border = false;
#endif
}

//
// private methods
//

//_____________________________________________________________________________
G4bool G4PlotManager::WritePage()
{
  fState.Message(kVL4, "write a page in", "plot file", fFileName);

#if defined(TOOLS_USE_FREETYPE)
  HD_style(fViewer->plots(), 5);
  regions_style(fViewer->plots(), fPlotParameters.GetScale());
#endif

  G4bool result = fViewer->write_page();
  if ( ! result ) {
    Warn("Cannot write a page in the plot file " + fFileName,
      fkClass, "WritePage");
  }

  // clear viewers plots
  fViewer->plots().init_sg();
    //it will recreate the sg::plotters and then reset the styles on new ones.

  fState.Message(kVL3, "write a page in", "plot file", fFileName);

  return result;
}

//
// public methods
//

//_____________________________________________________________________________
G4bool G4PlotManager::OpenFile(const G4String& fileName)
{
  fState.Message(kVL4, "open", "plot file", fileName);

  // Keep filename for logging
  fFileName = fileName;

  G4bool result = fViewer->open_file(fileName);
  if ( ! result ) {
    Warn("Cannot open plot file " + fileName, fkClass, "OpenFile");
  }

  fState.Message(kVL1, "open", "plot file", fileName);

  return result;
}

//_____________________________________________________________________________
G4bool G4PlotManager::CloseFile()
{
  fState.Message(kVL4, "close", "plot file", fFileName);

  G4bool result = fViewer->close_file();
  if ( ! result ) {
    Warn("Cannot close the plot file", fkClass, "CloseFile");
  }

  fState.Message(kVL1, "close", "plot file", fFileName);

  return result;
}
