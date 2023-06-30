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
///  \file FinalStateHistoManager.hh
///  \brief Create a set of histos for final state study.
//
//  Author: G.Hugo, 08 December 2022
//
// ***************************************************************************
//
//      FinalStateHistoManager
//
///  Create a set of histos for final state study.
///  In practice, the interactions studied here are hadron nuclear inelastic interactions
///  (though the code is fully generic).
///
///  Energy spectra are plotted for all encountered secondaries 
///  (one histo per secondary).
///  In addition, the residual nuclei Z and A distributions are plotted.
///
///  All histograms are G4H1. 
///  They are created and filled via the G4VAnalysisManager.
///
///  The histograms can be dumped to all usual formats, including ROOT 
///  (via G4VAnalysisManager).
///  An interesting added feature here, is that the plots, while being allocated 
///  and filled via G4VAnalysisManager, are also dumped 
///  in a Flair-compatible format (via tools::histo::flair).
///
///  NB 1: Note that instead of a hardcoded number associated to a hardcoded set of particles,
///  particle PDG IDs are used to index the histos. 
///  This allows a dynamic storage of all particles encountered in the final states.
///
///  NB 2: tools::histo::flair code, which allows the dump of any G4H1 
///  into Flair-compatible format, is fully application-agnostic, 
///  and is placed in FlukaCern/utils. 
///  It could also be added as an extension of core G4 Analysis Manager.
//
// ***************************************************************************

#ifndef FINAL_STATE_HISTO_MANAGER_HH
#define FINAL_STATE_HISTO_MANAGER_HH

#include <memory>
#include <unordered_map>
#include <vector>

#include "globals.hh"

#include "G4SystemOfUnits.hh"

#include "G4H1Wrapper.hh"


class G4DynamicParticle;
class G4VAnalysisManager;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FinalStateHistoManager {

public:
  FinalStateHistoManager();

  void Book();
  void BeginOfEvent();
  void ScoreSecondary(const G4DynamicParticle* const secondary);
  void EndOfEvent();
  void EndOfRun() const;


private:
  void DumpAllG4H1IntoRootFile() const;
  void DumpAllG4H1IntoFlairFile(const std::map<G4String, 
                                const G4H1Wrapper*>& particlesHistos) const;
  
  G4String fOutputFileName = "all_secondaries";
  G4String fRootOutputFileName = fOutputFileName + ".root";
  G4String fFlairOutputFileName = fOutputFileName + ".hist";

  G4int fNumBins = 90;
  G4double fMinKineticEnergy = 10. * keV;
  G4double fMaxKineticEnergy = 10. * TeV;
  G4String fFunctionName = "none";
  G4String fBinSchemeName = "log";
  G4String fRootEnergyUnit = "MeV";

  G4int fNucleiZMax = 25;
  G4int fNucleiAMax = 50;

  G4int fNumEvents = 0;

  G4VAnalysisManager* fAnalysisManager = nullptr;

  // key is particle PDG ID:
  std::unordered_map<G4int, std::unique_ptr<G4H1Wrapper>> fParticleData;
  // key is nuclei Z or A score index:
  std::unordered_map<G4int, std::unique_ptr<G4H1Wrapper>> fNucleiData;
  G4int fNucleiZScoreIndex = 0;
  G4int fNucleiAScoreIndex = 1;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....


#endif
