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
///  \file G4H1Wrapper.cc
///  \brief This allows statistics compatible with FLUKA.
///         It is fully relying on G4VAnalysisManager's G4H1.
//
//  Author: G.Hugo, 08 December 2022
//
// ***************************************************************************
//
//      G4H1Wrapper
//
///  This allows event-level statistics.
///  This is necessary IN THIS SPECIFIC CASE ONLY, to be able to compare errors on results
///  with e.g. FLUKA, MCNP, PENELOPE, SERPENT, which all use an event-based error approach.
///
///  It is fully relying on G4VAnalysisManager's G4H1.
///  Collects info within an event, then flush it to G4VAnalysisManager's G4H1.
//
// ***************************************************************************

#include "G4H1Wrapper.hh"

#include "G4RootAnalysisManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4H1Wrapper::G4H1Wrapper(G4VAnalysisManager* const analysisManager,
                         const G4int histoIndex) :
  fAnalysisManager(analysisManager),
  fHistoIndex(histoIndex),
  fHisto(dynamic_cast<G4RootAnalysisManager*>(analysisManager)->GetH1(histoIndex)),
  //fAxis(fHisto->axis()),
  fNumBins(fHisto->axis().bins()),
  fIsSet(false),
  fUnderflow(0.),
  fOverflow(0.),
  fSumSquaredEventTotals(0.),
  fSumSquaredEventInRangeTotals(0.)
{
  fEventData.assign(fNumBins, 0.);
}


// ***************************************************************************
// Fill data WITHIN THE EVENT. Can be called an arbitrary number of times.
// ***************************************************************************
void G4H1Wrapper::Fill(const G4double abscissaValue, 
                       const G4double weight) {
  if (!fIsSet) { fIsSet = true; }

  const G4int binIndex = fHisto->axis().coord_to_index(abscissaValue);

  // underflow
  // See external/g4tools/include/tools/histo/axis.
  if (binIndex == tools::histo::axis_UNDERFLOW_BIN) {
    fUnderflow += weight;
  }
  // overflow
  // See external/g4tools/include/tools/histo/axis.
  else if (binIndex == tools::histo::axis_OVERFLOW_BIN) {
    fOverflow += weight;
  }
  // in range
  else {
    fEventData[binIndex] += weight;
  }
}


// ***************************************************************************
// End of event: flush all data collected during the event 
// to G4VAnalysisManager's G4H1, then reset event data.
// ***************************************************************************
void G4H1Wrapper::EndOfEvent() {

  // Only do the flush if there was at least one fill during the event.
  if (fIsSet) {

    G4double eventTotal = 0.;
    // FLUSH IN RANGE DATA.
    for (G4int binIndex = 0; binIndex < fNumBins; ++binIndex) {

      if (fEventData[binIndex] > 0.) {

        const G4double binKineticEnergy = fHisto->axis().bin_center(binIndex);
        const G4double binEventWeight = fEventData[binIndex];
        fAnalysisManager->FillH1(fHistoIndex, 
                                 binKineticEnergy, 
                                 binEventWeight);
        eventTotal += binEventWeight;

        fEventData[binIndex] = 0.;        
      }
    }
    // Also keep track of the event's total (INTEGRAL over abscissa range).
    const G4double squaredEventTotal = eventTotal * eventTotal;
    fSumSquaredEventInRangeTotals += squaredEventTotal;
    fSumSquaredEventTotals += squaredEventTotal;

    // FLUSH UNDERFLOW DATA.
    if (fUnderflow > 0.) {
      fAnalysisManager->FillH1(fHistoIndex, 
                               fHisto->axis().lower_edge() - 1., 
                               fUnderflow);
      fSumSquaredEventTotals += fUnderflow;
      fUnderflow = 0.;
    }

    // FLUSH OVERFLOW DATA.
    if (fOverflow > 0.) {
      fAnalysisManager->FillH1(fHistoIndex, 
                               fHisto->axis().upper_edge() + 1., 
                               fOverflow);
      fSumSquaredEventTotals += fOverflow;
      fOverflow = 0.;
    }

    fIsSet = false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
