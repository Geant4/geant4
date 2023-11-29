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
///  \file G4H1Wrapper.hh
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

#ifndef G4H1_WRAPPER_HH
#define G4H1_WRAPPER_HH

#include <vector>

#include "globals.hh"
#include "g4hntools_defs.hh"

class G4VAnalysisManager;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4H1Wrapper {

public:
  G4H1Wrapper(G4VAnalysisManager* const analysisManager,
              const G4int histoIndex);

  void Fill(const G4double abscissaValue, 
            const G4double weight);
  void EndOfEvent();

  G4double GetSumSquaredEventTotals() const { return fSumSquaredEventTotals; }
  G4double GetSumSquaredEventInRangeTotals() const { return fSumSquaredEventInRangeTotals; }
  G4H1* GetG4H1() const { return fHisto; } // RUN DATA (G4H1 from G4VAnalysisManager)


private:
  G4VAnalysisManager* fAnalysisManager = nullptr;
  G4int fHistoIndex;
  G4H1* fHisto = nullptr;
  G4int fNumBins;

  G4bool fIsSet;
  G4double fUnderflow;
  G4double fOverflow;
  G4double fSumSquaredEventTotals;
  G4double fSumSquaredEventInRangeTotals;
  std::vector<G4double> fEventData; // EVENT DATA (data stored on the fly WITHIN the event)
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif
