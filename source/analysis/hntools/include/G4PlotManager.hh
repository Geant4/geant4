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

// The manager class for batch plotting.

// Author: Ivana Hrivnacova, 02/06/2015  (ivana@ipno.in2p3.fr)

#ifndef G4PlotManager_h
#define G4PlotManager_h 1

#include "G4AnalysisManagerState.hh"
#include "G4PlotParameters.hh"
#include "G4HnInformation.hh"

#include "tools/viewplot"

#include <vector>
#include <memory>
#include <string_view>

class G4PlotManager
{
  public:
    explicit G4PlotManager(const G4AnalysisManagerState& state);
    ~G4PlotManager() = default;

    // deleted functions
    G4PlotManager() = delete;
    G4PlotManager(const G4PlotManager& rhs) = delete;
    G4PlotManager& operator=(const G4PlotManager& rhs) = delete;

  public:
    // Methods
    G4bool OpenFile(const G4String& fileName);
    template <typename HT>
    G4bool PlotAndWrite(const std::vector<std::pair<HT*, G4HnInformation*>>& hnVector);
    G4bool CloseFile();

  private:
    // Methods
    G4int  GetNofPlotsPerPage() const;
    G4bool WritePage();

    // Static data members
    static constexpr std::string_view fkClass { "G4PlotManager" };

    // Data members
    const G4AnalysisManagerState& fState;
    G4PlotParameters fPlotParameters;
    std::unique_ptr<tools::viewplot>  fViewer;
    G4String  fFileName;
};

// inline functions

//_____________________________________________________________________________
inline G4int  G4PlotManager::GetNofPlotsPerPage() const
{ return fPlotParameters.GetColumns()*fPlotParameters.GetRows(); }


//_____________________________________________________________________________
template <typename HT>
inline G4bool G4PlotManager::PlotAndWrite(
  const std::vector<std::pair<HT*, G4HnInformation*>>& hnVector)
{
  if ( ! hnVector.size() ) return true;

  fViewer->plots().init_sg();
    //it will recreate the sg::plotters and then reset the styles on new ones.
  fViewer->set_cols_rows(fPlotParameters.GetColumns(), fPlotParameters.GetRows());
  fViewer->plots().set_current_plotter(0);

  auto result = true;
  auto isWriteNeeded = false;

  for (const auto& [ht, info] : hnVector) {
    G4bool plotting = info->GetPlotting();
    G4bool activation = info->GetActivation();
    G4String name = info->GetName();
    // skip plotting if not selected for plotting or
    // if activation is enabled and HT is inactivated
    if ( ( ! plotting ) ||
         ( fState.GetIsActivation() && ( ! activation ) ) ) continue;

    // plot this object
    fViewer->plot(*ht);
    fViewer->set_current_plotter_style(fPlotParameters.GetStyle());

    // set color (only blue for the time being)
    tools::sg::plotter& plotter = fViewer->plots().current_plotter();
    // set plot properties (use info object to get these)
    plotter.bins_style(0).color = tools::colorf_blue();

    // get axis titles from base_histo (base of all T)
    G4String title;
    if ( ht->annotation(tools::histo::key_axis_x_title(), title) ) {
      plotter.x_axis().title = title;
    }
    if ( ht->annotation(tools::histo::key_axis_y_title(), title) ) {
      plotter.y_axis().title = title;
    }
    if ( ht->annotation(tools::histo::key_axis_z_title(), title) ) {
      plotter.z_axis().title = title;
    }

#ifndef TOOLS_USE_FREETYPE
    plotter.set_encoding_none();
#endif

    // get log axis parameters from G4HnInformation
    if ( info->GetIsLogAxis(G4Analysis::kX) ) {
      plotter.x_axis().labels_style().encoding = "PAW";
      plotter.x_axis_is_log = true;
    }
    if ( info->GetIsLogAxis(G4Analysis::kY) ) {
      plotter.y_axis().labels_style().encoding = "PAW";
      plotter.y_axis_is_log = true;
    }
    if ( info->GetIsLogAxis(G4Analysis::kZ) ) {
      plotter.z_axis().labels_style().encoding = "PAW";
      plotter.z_axis_is_log = true;
    }
    isWriteNeeded = true;

    fState.Message(G4Analysis::kVL3, "plotting", "hd|pd", name);

    // write a page if number of plots per page is achieved
    if ( G4int(fViewer->plots().current_index()) == (GetNofPlotsPerPage() - 1) ) {
      result &= WritePage();
      isWriteNeeded = false;
    }

    // Prepare for the next plot
    fViewer->plots().next();
  }

  // write a page if loop is finished and there are plots to be written
  if ( isWriteNeeded ) {
    result &= WritePage();
  }

  // add test of result
  return result;
}

#endif
