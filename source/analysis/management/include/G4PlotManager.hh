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
// $Id$

// The manager class for batch plotting.

// Author: Ivana Hrivnacova, 02/06/2015  (ivana@ipno.in2p3.fr)

#ifndef G4PlotManager_h
#define G4PlotManager_h 1

#include "G4AnalysisManagerState.hh"
#include "G4PlotParameters.hh"
#include "G4HnInformation.hh"

#include <tools/viewplot>

#include <vector>
#include <memory>

class G4PlotManager
{
  public:
    explicit G4PlotManager(const G4AnalysisManagerState& state);
    ~G4PlotManager();

    // deleted functions
    G4PlotManager(const G4PlotManager& rhs) = delete;
    G4PlotManager& operator=(const G4PlotManager& rhs) = delete;
    
  public:
    // methods
    G4bool OpenFile(const G4String& fileName);
    template <typename T>
    G4bool PlotAndWrite(const std::vector<T*>& htVector,
                        const std::vector<G4HnInformation*>& hnVector);
    G4bool CloseFile();

  private: 
    // methods
    G4int  GetNofPlotsPerPage() const;
    G4bool WritePage();

    // static data members
    static G4PlotParameters fgPlotParameters;

    // data members
    const G4AnalysisManagerState& fState;
    std::unique_ptr<tools::viewplot>  fViewer;
    G4String  fFileName;
};

// inline functions

//_____________________________________________________________________________
inline G4int  G4PlotManager::GetNofPlotsPerPage() const
{ return fgPlotParameters.GetColumns()*fgPlotParameters.GetRows(); }


//_____________________________________________________________________________
template <typename T>
inline G4bool G4PlotManager::PlotAndWrite(const std::vector<T*>& htVector,
                                          const std::vector<G4HnInformation*>& hnVector)
{
  if ( ! htVector.size() ) return true;

  fViewer->plots().init_sg(); 
    //it will recreate the sg::plotters and then reset the styles on new ones.
  fViewer->set_cols_rows(fgPlotParameters.GetColumns(), fgPlotParameters.GetRows());
  fViewer->plots().set_current_plotter(0);

  G4bool finalResult = true;
  G4bool isWriteNeeded = false;

  for ( G4int i=0; i<G4int(htVector.size()); ++i ) {
    G4HnInformation* info = hnVector[i];
    G4bool plotting = info->GetPlotting();
    G4bool activation = info->GetActivation();
    G4String name = info->GetName();
    // skip plotting if not selected for plotting or
    // if activation is enabled and HT is inactivated
    if ( ( ! plotting ) ||
         ( fState.GetIsActivation() && ( ! activation ) ) ) continue;

    T* ht = htVector[i];

    // plot this object
    fViewer->plot(*ht);
    fViewer->set_current_plotter_style(fgPlotParameters.GetStyle());

    // set color (only blue for the time being)
    tools::sg::plotter& plotter = fViewer->plots().current_plotter();
    // set plot properties (use info object to get these)
    plotter.bins_style(0).color = tools::colorf_blue();
    
    isWriteNeeded = true;

#ifdef G4VERBOSE
    if ( fState.GetVerboseL3() ) 
      fState.GetVerboseL3()->Message("plotting", "hd|pd", name);
#endif

    // write a page if number of plots per page is achieved
    if ( G4int(fViewer->plots().current_index()) == (GetNofPlotsPerPage() - 1) ) { 
      G4bool result = WritePage();
      finalResult = result && finalResult;
      isWriteNeeded = false;
    }

    // Prepare for the next plot
    fViewer->plots().next(); 
  }

  // write a page if loop is finished and there are plots to be written
  if ( isWriteNeeded ) {
    G4bool result = WritePage();
    finalResult = result && finalResult;
  }

  // add test of result
  return finalResult;
}

#endif
