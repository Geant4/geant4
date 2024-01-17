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
///  \file XSHistoManager.cc
///  \brief Create a set of profiles for XS study.
//
//  Adapted from hadronic/Hadr00/src/HistoManager.cc
//  Author: G.Hugo, 06 January 2023
//
// ***************************************************************************
//
//      XSHistoManager
//
///  Create a set of profiles for XS study.
///
///  All profiles are G4H1. 
///  They are created and filled via G4VAnalysisManager.
///
///  The profiles can be dumped to all usual formats, including ROOT 
///  (via G4VAnalysisManager).
///  They are also dumped in a format compatible with Flair
///  (via tools::histo::flair).
//
// ***************************************************************************

#include "XSHistoManager.hh"

#include "G4ios.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4NistManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"

#include "G4HadronicProcessStore.hh"

//#include "G4AnalysisManager.hh"
#include "G4RootAnalysisManager.hh"

#include "tools_histo_flair.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XSHistoManager::XSHistoManager() :
  fMessenger(new XSHistoManagerMessenger(this)),
  fOutputFileName("all_XS"),
  fRootOutputFileName("all_XS.root"),
  fFlairOutputFileName("all_XS.hist"),
  fParticle(nullptr),
  fElement(nullptr),
  fMaterial(nullptr),
  fNumBins(10000),
  fMinKineticEnergy(1.*keV),
  fMaxKineticEnergy(10.*TeV),
  fFunctionName("none"),
  fBinSchemeName("log"),
  fRootEnergyUnit("MeV"),
  fAnalysisManager(G4RootAnalysisManager::Instance()),
  fElasticXSIndex(0),
  fInelasticXSIndex(1),
  fCaptureXSIndex(2),
  fFissionXSIndex(3),
  fChargeExchangeXSIndex(4),
  fTotalXSIndex(5),
  fElasticPerVolumeXSIndex(6),
  fInelasticPerVolumeXSIndex(7)
{
  //G4NistManager::Instance()->ListMaterials("all");
}


// ***************************************************************************
// Set output files names: 2 formats supported, ROOT and Flair.
// ***************************************************************************
void XSHistoManager::SetOutputFileName(const G4String& outputFileName) { 
  fOutputFileName = outputFileName;
  fRootOutputFileName = outputFileName + ".root";
  fFlairOutputFileName = outputFileName + ".hist";
}


// ***************************************************************************
// Set the particle considered for XS study.
// ***************************************************************************
void XSHistoManager::SetParticle(const G4String& particleName) { 
  fParticle = G4ParticleTable::GetParticleTable()->FindParticle(particleName);
}


// ***************************************************************************
// Set the target element considered for XS study.
// ***************************************************************************
void XSHistoManager::SetElement(const G4String& elementName) {
  fElement = G4NistManager::Instance()->FindOrBuildElement(elementName);
  // Also needs to set material!
  SetMaterial(elementName);
}


// ***************************************************************************
// Set the target material considered for XS study.
// ***************************************************************************
void XSHistoManager::SetMaterial(const G4String& materialName) {

  // Check that material is not set already.
  if (fMaterial) {
    G4ExceptionDescription msg;
    msg << "Please use UI command /allXS/elementName"
        << " OR UI command /allXS/nonElementaryMaterialName,"
        << " BUT NOT BOTH!"
        << G4endl;
    G4Exception("XSHistoManager::SetMaterial",
                "A target material is already defined.",
                FatalException,
                msg);
  }

  fMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_" + materialName);
}


// ***************************************************************************
// Open output file + create all profiles considered for XS study.
// All profiles are G4H1, created via G4VAnalysisManager.
// ***************************************************************************
void XSHistoManager::Book() {

  // Check all XSHistoManager data is set properly.
  CheckInput();

  // Open file.
  if (!fAnalysisManager->OpenFile(fRootOutputFileName)) {

    G4ExceptionDescription msg;
    msg << "Booking profiles: cannot open file " << fRootOutputFileName 
        << G4endl;
    G4Exception("XSHistoManager::Book",
                "Cannot open file",
                FatalException,
                msg);
  }
  G4cout << "### XSHistoManager::Book: Successfully opened file "
         << fRootOutputFileName 
         << " for dumping profiles." 
         << G4endl;

  // Create all G4H1, and keep track of each histo index in fXSProfileIndex.
  const G4int elasticXSProfileIndex = fAnalysisManager->CreateH1("ElasticXS", 
                                                             "Elastic XS", 
                                                             fNumBins, 
                                                             fMinKineticEnergy,
                                                             fMaxKineticEnergy,
                                                             fRootEnergyUnit,
                                                             fFunctionName,
                                                             fBinSchemeName);
  fXSProfileIndex.insert(std::make_pair(fElasticXSIndex, 
                                        elasticXSProfileIndex));

  const G4int inelasticXSProfileIndex = fAnalysisManager->CreateH1("InelasticXS", 
                                                               "Inelastic XS", 
                                                               fNumBins, 
                                                               fMinKineticEnergy,
                                                               fMaxKineticEnergy,
                                                               fRootEnergyUnit,
                                                               fFunctionName,
                                                               fBinSchemeName);
  fXSProfileIndex.insert(std::make_pair(fInelasticXSIndex, 
                                        inelasticXSProfileIndex));

  const G4int captureXSProfileIndex = fAnalysisManager->CreateH1("CaptureXS", 
                                                             "Capture XS", 
                                                             fNumBins, 
                                                             fMinKineticEnergy,
                                                             fMaxKineticEnergy,
                                                             fRootEnergyUnit,
                                                             fFunctionName,
                                                             fBinSchemeName);
  fXSProfileIndex.insert(std::make_pair(fCaptureXSIndex, 
                                        captureXSProfileIndex));

  const G4int fissionXSProfileIndex = fAnalysisManager->CreateH1("FissionXS", 
                                                             "Fission XS", 
                                                             fNumBins, 
                                                             fMinKineticEnergy,
                                                             fMaxKineticEnergy,
                                                             fRootEnergyUnit,
                                                             fFunctionName,
                                                             fBinSchemeName);
  fXSProfileIndex.insert(std::make_pair(fFissionXSIndex, 
                                        fissionXSProfileIndex));

  const G4int chargeExchangeXSProfileIndex = fAnalysisManager->CreateH1("ChargeExchangeXS", 
                                                                    "Charge exchange XS", 
                                                                    fNumBins, 
                                                                    fMinKineticEnergy,
                                                                    fMaxKineticEnergy,
                                                                    fRootEnergyUnit,
                                                                    fFunctionName,
                                                                    fBinSchemeName);
  fXSProfileIndex.insert(std::make_pair(fChargeExchangeXSIndex, 
                                        chargeExchangeXSProfileIndex));

  const G4int totalXSProfileIndex = fAnalysisManager->CreateH1("TotalXS", 
                                                           "Total XS", 
                                                           fNumBins, 
                                                           fMinKineticEnergy,
                                                           fMaxKineticEnergy,
                                                           fRootEnergyUnit,
                                                           fFunctionName,
                                                           fBinSchemeName);
  fXSProfileIndex.insert(std::make_pair(fTotalXSIndex,
                                        totalXSProfileIndex));

  const G4int elasticPerVolumeXSProfileIndex = fAnalysisManager->CreateH1("ElasticPerVolumeXS", 
                                                                      "Elastic XS per volume", 
                                                                      fNumBins, 
                                                                      fMinKineticEnergy,
                                                                      fMaxKineticEnergy,
                                                                      fRootEnergyUnit,
                                                                      fFunctionName,
                                                                      fBinSchemeName);
  fXSProfileIndex.insert(std::make_pair(fElasticPerVolumeXSIndex, 
                                        elasticPerVolumeXSProfileIndex));

  const G4int inelasticPerVolumeXSProfileIndex =
    fAnalysisManager->CreateH1("InelasticPerVolumeXS",
                               "Inelastic XS per volume",
                               fNumBins,
                               fMinKineticEnergy,
                               fMaxKineticEnergy,
                               fRootEnergyUnit,
                               fFunctionName,
                               fBinSchemeName);
  fXSProfileIndex.insert(std::make_pair(fInelasticPerVolumeXSIndex, 
                                        inelasticPerVolumeXSProfileIndex));
}


// ***************************************************************************
// Fill all plots, then dump them into relevant formats.
// ***************************************************************************
void XSHistoManager::EndOfRun() {
        G4cout << "### XSHistoManager::EndOfRun: Compute & fill XS for " 
               << fParticle->GetParticleName() 
               << " in " << (fElement ? fElement->GetName() : fMaterial->GetName())
               << G4endl;

        G4HadronicProcessStore* const store = G4HadronicProcessStore::Instance();

        // Fill XS profiles.
        const G4double logMinKineticEnergy = std::log10(fMinKineticEnergy);
        const G4double logMaxKineticEnergy = std::log10(fMaxKineticEnergy);
        const G4double deltaLogKineticEnergy =
          (logMaxKineticEnergy - logMinKineticEnergy) / fNumBins;

        G4double logKineticEnergy = logMinKineticEnergy - deltaLogKineticEnergy/2.;
        
        // Loop on all kinetic energies of interest.
        for (G4int binIndex = 0; binIndex < fNumBins; ++binIndex) {

          logKineticEnergy += deltaLogKineticEnergy;
          const G4double kineticEnergy = std::pow(10., logKineticEnergy) * MeV;

          G4double totalXS = 0.;
          if (fElement) {
            // ELASTIC (ELEMENTARY MATERIAL)
            const G4double elasticXS = store->GetElasticCrossSectionPerAtom(fParticle,
                                                                          kineticEnergy,
                                                                          fElement,
                                                                          fMaterial);         
            fAnalysisManager->FillH1(fXSProfileIndex[fElasticXSIndex], 
                                     kineticEnergy, 
                                     elasticXS/barn);
            totalXS += elasticXS;
 
            // INELASTIC (ELEMENTARY MATERIAL)
            const G4double inelasticXS = store->GetInelasticCrossSectionPerAtom(fParticle,
                                                                              kineticEnergy,
                                                                              fElement,
                                                                              fMaterial);
            fAnalysisManager->FillH1(fXSProfileIndex[fInelasticXSIndex], 
                                     kineticEnergy, 
                                     inelasticXS/barn);
            totalXS += inelasticXS;

            if (fParticle == G4Neutron::Definition()) {
              // NEUTRON CAPTURE (ELEMENTARY MATERIAL)
              const G4double captureXS = store->GetCaptureCrossSectionPerAtom(fParticle,
                                                                            kineticEnergy,
                                                                            fElement,
                                                                            fMaterial);
              fAnalysisManager->FillH1(fXSProfileIndex[fCaptureXSIndex], 
                                       kineticEnergy, 
                                       captureXS/barn);
              totalXS += captureXS;
   
              // FISSION (ELEMENTARY MATERIAL)
              const G4double fissionXS = store->GetFissionCrossSectionPerAtom(fParticle,
                                                                            kineticEnergy,
                                                                            fElement,
                                                                            fMaterial);
              totalXS += fissionXS;
              fAnalysisManager->FillH1(fXSProfileIndex[fFissionXSIndex], 
                                       kineticEnergy, 
                                       fissionXS/barn);
            }

            // CHARGE EXCHANGE (ELEMENTARY MATERIAL)
            const G4double chargeExchangeXS = 
              store->GetChargeExchangeCrossSectionPerAtom(fParticle,
                                                          kineticEnergy,
                                                          fElement,
                                                          fMaterial);
            fAnalysisManager->FillH1(fXSProfileIndex[fChargeExchangeXSIndex], 
                                     kineticEnergy, 
                                     chargeExchangeXS/barn);
            totalXS += chargeExchangeXS;
   
            // TOTAL (ELEMENTARY MATERIAL)
            fAnalysisManager->FillH1(fXSProfileIndex[fTotalXSIndex], 
                                     kineticEnergy, 
                                     totalXS/barn);
          }

          if (fMaterial) {
            const G4double materialSurfacicDensity = (fMaterial ? 
                                                    fMaterial->GetDensity() / (g/cm2) 
                                                    : 1.);

            // ELASTIC
            const G4double elasticPerVolumeXS = 
              store->GetElasticCrossSectionPerVolume(fParticle,
                                                     kineticEnergy,
                                                     fMaterial);
            fAnalysisManager->FillH1(fXSProfileIndex[fElasticPerVolumeXSIndex], 
                                     kineticEnergy, 
                                     elasticPerVolumeXS/materialSurfacicDensity);

            // INELASTIC
            const G4double inelasticPerVolumeXS = 
              store->GetInelasticCrossSectionPerVolume(fParticle,
                                                       kineticEnergy,
                                                       fMaterial);
            fAnalysisManager->FillH1(fXSProfileIndex[fInelasticPerVolumeXSIndex], 
                                     kineticEnergy, 
                                     inelasticPerVolumeXS/materialSurfacicDensity);
          }

        }


        // DUMP G4H1 PLOTS INTO ROOT FILE
        DumpAllG4H1IntoRootFile();

        // DUMP G4H1 PLOTS INTO FLAIR FILE
        DumpAllG4H1IntoFlairFile();


        // Close and clear fAnalysisManager.
        fAnalysisManager->CloseFile();
        fAnalysisManager->Clear();
}


// ***************************************************************************
// Checks that particle and material are set 
// (all others have relevant default values).
// ***************************************************************************
void XSHistoManager::CheckInput() {

  if (!fParticle) {
    G4ExceptionDescription msg;
    msg << "Please add a particle to study XS: UI command /allXS/particleName"
        << G4endl;
    G4Exception("XSHistoManager::CheckInput()",
                "Print XS: no input particle defined.",
                FatalException,
                msg);
  }

  if (!fMaterial) {
    G4ExceptionDescription msg;
    msg << "Please add a material to study XS:"
        << " UI command /allXS/elementName for an elementary material,"
        << " or UI command /allXS/nonElementaryMaterialName for a compound/mixture material."
        << G4endl;
    G4Exception("XSHistoManager::CheckInput()",
                "Print XS: no target material defined.",
                FatalException,
                msg);
  }
}


// ***************************************************************************
// DUMP G4H1 PLOTS INTO ROOT FILE (via G4VAnalysisManager).
// ***************************************************************************
void XSHistoManager::DumpAllG4H1IntoRootFile() const {

        if (!fAnalysisManager->Write()) {
          G4ExceptionDescription message;
          message << "Could not write ROOT file."; 
          G4Exception("XSHistoManager::EndOfRun()",
                      "I/O Error", 
                      FatalException, 
                      message);
        }
        G4cout << "### All profiles saved to " << fRootOutputFileName << G4endl;
}


// ***************************************************************************
// DUMP G4H1 PLOTS INTO FLAIR FILE (via tools::histo::flair).
// ***************************************************************************
void XSHistoManager::DumpAllG4H1IntoFlairFile() const {

  std::ofstream output;
  output.open(fFlairOutputFileName, std::ios_base::out);
  auto const rootAnalysisManager = dynamic_cast<G4RootAnalysisManager*>(fAnalysisManager);

  G4int indexInOutputFile = 1;
  for (G4int xsIndex = fElasticXSIndex; xsIndex <= fInelasticPerVolumeXSIndex; ++xsIndex) {

    const G4int histoIndex = fXSProfileIndex.at(xsIndex);
    const G4String& histoName = fAnalysisManager->GetH1Name(histoIndex);
    const auto& histo = rootAnalysisManager->GetH1(histoIndex);
    
    tools::histo::flair::dumpG4H1ProfileInFlairFormat(output,
                                                      indexInOutputFile,
                                                      histoName,
                                                      histo,
                                                      tools::histo::flair::Abscissa::KineticEnergy,
                                                      fBinSchemeName);
    ++indexInOutputFile;
  }
  output.close();
  G4cout << "### All profiles saved to " << fFlairOutputFileName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
