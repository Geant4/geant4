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
///  They are created and filled solely via G4VAnalysisManager.
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

#include "FinalStateHistoManager.hh"

#include "G4RootAnalysisManager.hh"
//#include "G4AnalysisManager.hh"

#include "G4ParticleTable.hh"
#include "G4DynamicParticle.hh"

#include "G4ios.hh"
#include "G4Exception.hh"

#include "g4hntools_defs.hh"
#include "tools_histo_flair.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FinalStateHistoManager::FinalStateHistoManager() :
  fOutputFileName("all_secondaries"),
  fRootOutputFileName(fOutputFileName + ".root"),
  fFlairOutputFileName(fOutputFileName + ".hist"),
  fNumBins(90),
  fMinKineticEnergy(10. * keV),
  fMaxKineticEnergy(10. * TeV),
  fFunctionName("none"),
  fBinSchemeName("log"),
  fRootEnergyUnit("MeV"),
  fNucleiZMax(25),
  fNucleiAMax(50),
  fNumEvents(0),
  fAnalysisManager(G4RootAnalysisManager::Instance()),
  fNucleiZScoreIndex(0),
  fNucleiAScoreIndex(1)
{
  //fAnalysisManager = G4AnalysisManager::Instance();
  //fAnalysisManager->SetDefaultFileType("root");
  //fAnalysisManager->SetVerboseLevel(0);
  //fOutputFileName += fAnalysisManager->GetFileType();
}


// ***************************************************************************
// Open output file + create residual nuclei histograms considered for final state study.
// The histograms are G4H1, created via G4VAnalysisManager.
// ***************************************************************************
void FinalStateHistoManager::Book() {

  // Open file.
  if(!fAnalysisManager->OpenFile(fRootOutputFileName)) {

    G4ExceptionDescription msg;
    msg << "Booking histograms: cannot open file " 
        << fRootOutputFileName 
        << G4endl;
    G4Exception("FinalStateHistoManager::Book",
                "Cannot open file",
                FatalException,
                msg);
  }
  G4cout << "### FinalStateHistoManager::Book: Successfully opended file " 
         << fRootOutputFileName 
         << " for dumping histograms." 
         << G4endl;


  // Create the residual nuclei distributions (in Z and A).
  const G4int nucleiZHistoIndex = fAnalysisManager->CreateH1("nucleiZ", 
                                                           "Residual nuclei distribution in Z", 
                                                           fNucleiZMax, 
                                                           0.5,
                                                           fNucleiZMax + 0.5);
  auto nucleiZHistoWrapper = std::make_unique<G4H1Wrapper>(fAnalysisManager, 
                                                           nucleiZHistoIndex);
  fNucleiData.insert(std::make_pair(fNucleiZScoreIndex, std::move(nucleiZHistoWrapper)));
                                                           

  const G4int nucleiAHistoIndex = fAnalysisManager->CreateH1("nucleiA", 
                                                           "Residual nuclei distribution in A", 
                                                           fNucleiAMax, 
                                                           0.5,
                                                           fNucleiAMax + 0.5);
  auto nucleiAHistoWrapper = std::make_unique<G4H1Wrapper>(fAnalysisManager, 
                                                           nucleiAHistoIndex);
  fNucleiData.insert(std::make_pair(fNucleiAScoreIndex, std::move(nucleiAHistoWrapper)));
}


// ***************************************************************************
// Keep track of the total number of events (used later on for normalization).
// ***************************************************************************
void FinalStateHistoManager::BeginOfEvent() {
  fNumEvents++;
}


// ***************************************************************************
// Fill all plots (WITHIN event, ie the interaction).
// ***************************************************************************
void FinalStateHistoManager::ScoreSecondary(const G4DynamicParticle* const secondary) {

  // SELECT SPECIFIC SECONDARIES ONLY
  // Select by angle with beam direction
  /* if ( (std::pow(secondary->GetMomentumDirection().x(), 2.) 
     + std::pow(secondary->GetMomentumDirection().y(), 2.))
     <= 0.0001 ) {*/

  // Select by production tag
  /* if (secondary->GetProductionTag() == 6) {*/

  // Primary track
  /* if(track->GetParentID() == 0) {*/

  const auto& particle = secondary->GetDefinition();


  // SECONDARIES ENERGY SPECTRA
  // Dynamic creation of histos, so that all encountered particles have their own histos.

  // Check whether a particle has already been encountered.
  const auto found = fParticleData.find(secondary->GetPDGcode());

  G4H1Wrapper* particleHistoWrapper = nullptr;
  // If the particle has already been encountered, use the corresponding histos.
  if (found != fParticleData.end()) {
    particleHistoWrapper = found->second.get();
  }
  // Otherwise, create histos for that particle.
  else {
    const G4String& particleName = particle->GetParticleName();
    const G4int particlePDG = secondary->GetPDGcode();

    const G4String histoTitle = (particlePDG == 0 ? 
                                 "Particle pdg==0 spectrum"
                                 : G4String(particleName + " spectrum"));
    const G4int histoIndex = fAnalysisManager->CreateH1(particleName, 
                                                      histoTitle,
                                                      fNumBins, 
                                                      fMinKineticEnergy,
                                                      fMaxKineticEnergy,
                                                      fRootEnergyUnit,
                                                      fFunctionName,
                                                      fBinSchemeName);
    auto histoWrapper = std::make_unique<G4H1Wrapper>(fAnalysisManager, 
                                                      histoIndex);
    particleHistoWrapper = histoWrapper.get();
    fParticleData.insert(std::make_pair(particlePDG, std::move(histoWrapper)));
  }

  // Fill the G4H1Wrapper.
  const G4double kineticEnergy = secondary->GetKineticEnergy();
  particleHistoWrapper->Fill(kineticEnergy, 1.);


  // NUCLEI DISTRIBUTIONS IN Z AND A
  if (particle->GetParticleType() == "nucleus") {
    // Fill the G4H1Wrapper.
    const G4double Z = particle->GetPDGCharge() / eplus;
    fNucleiData[fNucleiZScoreIndex]->Fill(Z, 1.);

    // Fill the G4H1Wrapper.
    const G4double A = particle->GetBaryonNumber();
    fNucleiData[fNucleiAScoreIndex]->Fill(A, 1.);
  }

  //} // select secondaries
}


// ***************************************************************************
// End of event: all event-level G4H1 are flushed into the Analysis Manager G4H1.
// ***************************************************************************
void FinalStateHistoManager::EndOfEvent() {

  for (const auto& particleIt : fParticleData) {
    particleIt.second->EndOfEvent();
  }
  for (const auto& nucleiScoreIt : fNucleiData) {
    nucleiScoreIt.second->EndOfEvent();
  }
}


// ***************************************************************************
// Printout secondary counts + dump all plots into relevant formats.
// ***************************************************************************
void FinalStateHistoManager::EndOfRun() const {

  // PRINTOUT SECONDARYS COUNTS (FULL ENERGY RANGE).

  // Order the histos by particles names.
  std::map<G4String, const G4H1Wrapper*> particlesHistos;

  for (const auto& particleIt : fParticleData) {
    const G4int particlePdg = particleIt.first;
    const G4String particleName = G4ParticleTable::GetParticleTable()
      ->FindParticle(particlePdg)->GetParticleName();

    const G4H1Wrapper* const particleHisto = particleIt.second.get();
    particlesHistos.insert(std::make_pair(particleName, particleHisto));
  }
  
  // Printout secondarys counts (full energy range)
  // Values are averaged over the number of events.
  G4cout << "========================================================" << G4endl;
  G4cout << "Number of events                     " << fNumEvents << G4endl << G4endl;
  for (const auto& particleIt : particlesHistos) {

    // Note that the info is directly obtained from the histogram:
    // it is the integral over the full energy range.
    const G4int count = particleIt.second->GetG4H1()->sum_all_bin_heights();
   
    const G4double averageCount = static_cast<G4double>(count) / fNumEvents;

    G4cout << "Average (per event) number of " << particleIt.first
           << "              " << averageCount
           << G4endl;
  }
  G4cout << "========================================================" << G4endl;
  G4cout << G4endl;


  // DUMP G4H1 PLOTS INTO ROOT FILE
  DumpAllG4H1IntoRootFile();

  // DUMP G4H1 PLOTS INTO FLAIR FILE
  DumpAllG4H1IntoFlairFile(particlesHistos);


  // Close and clear fAnalysisManager.
  fAnalysisManager->CloseFile();
  fAnalysisManager->Clear();
}


// ***************************************************************************
// DUMP G4H1 PLOTS INTO ROOT FILE (via G4VAnalysisManager).
// ***************************************************************************
void FinalStateHistoManager::DumpAllG4H1IntoRootFile() const {

  if (!fAnalysisManager->Write()) {
    G4ExceptionDescription message;
    message << "Could not write ROOT file."; 
    G4Exception("FinalStateHistoManager::EndOfRun()",
                "I/O Error", 
                FatalException, 
                message);
  }
  G4cout << "### All histograms saved to " << fRootOutputFileName << G4endl;
}


// ***************************************************************************
// DUMP G4H1 PLOTS INTO FLAIR FILE (via tools::histo::flair).
// ***************************************************************************
void FinalStateHistoManager::DumpAllG4H1IntoFlairFile(
         const std::map<G4String, const G4H1Wrapper*>& particlesHistos) const {

  std::ofstream output;
  output.open(fFlairOutputFileName, std::ios_base::out);
  G4int indexInOutputFile = 1;

  // SECONDARIES ENERGY SPECTRA
  for (const auto& particleIt : particlesHistos) {

    const G4String& histoName = particleIt.first;
    const auto& histo = particleIt.second->GetG4H1();
    
    tools::histo::flair::dumpG4H1HistoInFlairFormat(output,
                                                    indexInOutputFile,
                                                    histoName,
                                                    histo,
                                                    tools::histo::flair::Abscissa::KineticEnergy,
                                                    fBinSchemeName,
                                                    fNumEvents,
                                                    particleIt.second
                                                    ->GetSumSquaredEventTotals(),
                                                    particleIt.second
                                                    ->GetSumSquaredEventInRangeTotals());
    ++indexInOutputFile;
  }

  // RESIDUAL NUCLEI DISTRIBUTIONS
  for (const auto& plotIt : fNucleiData) {
  
    const auto& histo = plotIt.second->GetG4H1();
    const G4String& histoName = (plotIt.first == fNucleiZScoreIndex ? 
                                 "nucleiZ" 
                                 : "nucleiA");
    const auto& abscissaKind = (plotIt.first == fNucleiZScoreIndex ? 
                                tools::histo::flair::Abscissa::Z 
                                : tools::histo::flair::Abscissa::A);
    
    tools::histo::flair::dumpG4H1HistoInFlairFormat(output,
                                                    indexInOutputFile,
                                                    histoName,
                                                    histo,
                                                    abscissaKind,
                                                    fBinSchemeName,
                                                    fNumEvents,   
                                                    plotIt.second
                                                    ->GetSumSquaredEventTotals(),
                                                    plotIt.second
                                                    ->GetSumSquaredEventInRangeTotals());
    ++indexInOutputFile;
  }

  output.close();
  G4cout << "### All histograms saved to " << fFlairOutputFileName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
