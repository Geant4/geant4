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
/// \file optical/OpNovice2/src/Run.cc
/// \brief Implementation of the Run class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"

#include "DetectorConstruction.hh"
#include "HistoManager.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include <numeric>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Run::Run()
  : G4Run()
{
  fParticle     = nullptr;
  fEkin         = -1.;
  fPolarized    = false;
  fPolarization = 0.;

  fCerenkovEnergy       = 0.0;
  fScintEnergy          = 0.0;
  fWLSAbsorptionEnergy  = 0.0;
  fWLSEmissionEnergy    = 0.0;
  fWLS2AbsorptionEnergy = 0.0;
  fWLS2EmissionEnergy   = 0.0;

  fCerenkovCount       = 0;
  fScintCount          = 0;
  fWLSAbsorptionCount  = 0;
  fWLSEmissionCount    = 0;
  fWLS2AbsorptionCount = 0;
  fWLS2EmissionCount   = 0;
  fRayleighCount       = 0;

  fOpAbsorption      = 0;
  fOpAbsorptionPrior = 0;

  fTotalSurface = 0;

  fBoundaryProcs.assign(43, 0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Run::~Run() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy,
                     G4bool polarized, G4double polarization)
{
  fParticle     = particle;
  fEkin         = energy;
  fPolarized    = polarized;
  fPolarization = polarization;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);

  // pass information about primary particle
  fParticle     = localRun->fParticle;
  fEkin         = localRun->fEkin;
  fPolarized    = localRun->fPolarized;
  fPolarization = localRun->fPolarization;

  fCerenkovEnergy += localRun->fCerenkovEnergy;
  fScintEnergy += localRun->fScintEnergy;
  fWLSAbsorptionEnergy += localRun->fWLSAbsorptionEnergy;
  fWLSEmissionEnergy += localRun->fWLSEmissionEnergy;
  fWLS2AbsorptionEnergy += localRun->fWLS2AbsorptionEnergy;
  fWLS2EmissionEnergy += localRun->fWLS2EmissionEnergy;

  fCerenkovCount += localRun->fCerenkovCount;
  fScintCount += localRun->fScintCount;
  fWLSAbsorptionCount += localRun->fWLSAbsorptionCount;
  fWLSEmissionCount += localRun->fWLSEmissionCount;
  fWLS2AbsorptionCount += localRun->fWLS2AbsorptionCount;
  fWLS2EmissionCount += localRun->fWLS2EmissionCount;
  fRayleighCount += localRun->fRayleighCount;

  fTotalSurface += localRun->fTotalSurface;

  fOpAbsorption += localRun->fOpAbsorption;
  fOpAbsorptionPrior += localRun->fOpAbsorptionPrior;

  for(size_t i = 0; i < fBoundaryProcs.size(); ++i)
  {
    fBoundaryProcs[i] += localRun->fBoundaryProcs[i];
  }

  G4Run::Merge(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Run::EndOfRun()
{
  if(numberOfEvent == 0)
    return;
  G4double TotNbofEvents = (G4double) numberOfEvent;

  G4AnalysisManager* analysisMan = G4AnalysisManager::Instance();
  G4int id                       = analysisMan->GetH1Id("Cerenkov spectrum");
  analysisMan->SetH1XAxisTitle(id, "Energy [eV]");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("Scintillation spectrum");
  analysisMan->SetH1XAxisTitle(id, "Energy [eV]");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("Scintillation time");
  analysisMan->SetH1XAxisTitle(id, "Creation time [ns]");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("WLS abs");
  analysisMan->SetH1XAxisTitle(id, "Energy [eV]");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("WLS em");
  analysisMan->SetH1XAxisTitle(id, "Energy [eV]");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("WLS time");
  analysisMan->SetH1XAxisTitle(id, "Creation time [ns]");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("WLS2 abs");
  analysisMan->SetH1XAxisTitle(id, "Energy [eV]");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("WLS2 em");
  analysisMan->SetH1XAxisTitle(id, "Energy [eV]");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("WLS2 time");
  analysisMan->SetH1XAxisTitle(id, "Creation time [ns]");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("bdry status");
  analysisMan->SetH1XAxisTitle(id, "Status code");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("x_backward");
  analysisMan->SetH1XAxisTitle(id, "Direction cosine");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("y_backward");
  analysisMan->SetH1XAxisTitle(id, "Direction cosine");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("z_backward");
  analysisMan->SetH1XAxisTitle(id, "Direction cosine");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("x_forward");
  analysisMan->SetH1XAxisTitle(id, "Direction cosine");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("y_forward");
  analysisMan->SetH1XAxisTitle(id, "Direction cosine");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("z_forward");
  analysisMan->SetH1XAxisTitle(id, "Direction cosine");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("x_fresnel");
  analysisMan->SetH1XAxisTitle(id, "Direction cosine");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("y_fresnel");
  analysisMan->SetH1XAxisTitle(id, "Direction cosine");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("z_fresnel");
  analysisMan->SetH1XAxisTitle(id, "Direction cosine");
  analysisMan->SetH1YAxisTitle(id, "Number of photons");

  id = analysisMan->GetH1Id("Fresnel reflection");
  analysisMan->SetH1XAxisTitle(id, "Angle [deg]");
  analysisMan->SetH1YAxisTitle(id, "Fraction of photons");

  id = analysisMan->GetH1Id("Fresnel refraction");
  analysisMan->SetH1XAxisTitle(id, "Angle [deg]");
  analysisMan->SetH1YAxisTitle(id, "Fraction of photons");

  id = analysisMan->GetH1Id("Total internal reflection");
  analysisMan->SetH1XAxisTitle(id, "Angle [deg]");
  analysisMan->SetH1YAxisTitle(id, "Fraction of photons");

  id = analysisMan->GetH1Id("Fresnel reflection plus TIR");
  analysisMan->SetH1XAxisTitle(id, "Angle [deg]");
  analysisMan->SetH1YAxisTitle(id, "Fraction of photons");

  id = analysisMan->GetH1Id("Absorption");
  analysisMan->SetH1XAxisTitle(id, "Angle [deg]");
  analysisMan->SetH1YAxisTitle(id, "Fraction of photons");

  id = analysisMan->GetH1Id("Transmitted");
  analysisMan->SetH1XAxisTitle(id, "Angle [deg]");
  analysisMan->SetH1YAxisTitle(id, "Fraction of photons");

  id = analysisMan->GetH1Id("Spike reflection");
  analysisMan->SetH1XAxisTitle(id, "Angle [deg]");
  analysisMan->SetH1YAxisTitle(id, "Fraction of photons");

  const DetectorConstruction* det =
    (const DetectorConstruction*) (G4RunManager::GetRunManager()
                                     ->GetUserDetectorConstruction());

  std::ios::fmtflags mode = G4cout.flags();
  G4int prec              = G4cout.precision(2);

  G4cout << "\n    Run Summary\n";
  G4cout << "---------------------------------\n";
  G4cout << "Primary particle was: " << fParticle->GetParticleName()
         << " with energy " << G4BestUnit(fEkin, "Energy") << "." << G4endl;
  G4cout << "Number of events: " << numberOfEvent << G4endl;

  G4cout << "Material of world: " << det->GetWorldMaterial()->GetName()
         << G4endl;
  G4cout << "Material of tank:  " << det->GetTankMaterial()->GetName() << G4endl
         << G4endl;

  if(fParticle->GetParticleName() != "opticalphoton")
  {
    G4cout << "Average energy of Cerenkov photons created per event: "
           << (fCerenkovEnergy / eV) / TotNbofEvents << " eV." << G4endl;
    G4cout << "Average number of Cerenkov photons created per event: "
           << fCerenkovCount / TotNbofEvents << G4endl;
    if(fCerenkovCount > 0)
    {
      G4cout << " Average energy per photon: "
             << (fCerenkovEnergy / eV) / fCerenkovCount << " eV." << G4endl;
    }
    G4cout << "Average energy of scintillation photons created per event: "
           << (fScintEnergy / eV) / TotNbofEvents << " eV." << G4endl;
    G4cout << "Average number of scintillation photons created per event: "
           << fScintCount / TotNbofEvents << G4endl;
    if(fScintCount > 0)
    {
      G4cout << " Average energy per photon: "
             << (fScintEnergy / eV) / fScintCount << " eV." << G4endl;
    }
  }

  G4cout << "Average number of photons absorbed by WLS per event: "
         << fWLSAbsorptionCount / G4double(TotNbofEvents) << " " << G4endl;
  if(fWLSAbsorptionCount > 0)
  {
    G4cout << " Average energy per photon: "
           << (fWLSAbsorptionEnergy / eV) / fWLSAbsorptionCount << " eV."
           << G4endl;
  }
  G4cout << "Average number of photons created by WLS per event: "
         << fWLSEmissionCount / TotNbofEvents << G4endl;
  if(fWLSEmissionCount > 0)
  {
    G4cout << " Average energy per photon: "
           << (fWLSEmissionEnergy / eV) / fWLSEmissionCount << " eV." << G4endl;
  }
  G4cout << "Average energy of WLS photons created per event: "
         << (fWLSEmissionEnergy / eV) / TotNbofEvents << " eV." << G4endl;

  G4cout << "Average number of photons absorbed by WLS2 per event: "
         << fWLS2AbsorptionCount / G4double(TotNbofEvents) << " " << G4endl;
  if(fWLS2AbsorptionCount > 0)
  {
    G4cout << " Average energy per photon: "
           << (fWLS2AbsorptionEnergy / eV) / fWLS2AbsorptionCount << " eV."
           << G4endl;
  }
  G4cout << "Average number of photons created by WLS2 per event: "
         << fWLS2EmissionCount / TotNbofEvents << G4endl;
  if(fWLS2EmissionCount > 0)
  {
    G4cout << " Average energy per photon: "
           << (fWLS2EmissionEnergy / eV) / fWLS2EmissionCount << " eV."
           << G4endl;
  }
  G4cout << "Average energy of WLS2 photons created per event: "
         << (fWLS2EmissionEnergy / eV) / TotNbofEvents << " eV." << G4endl;

  G4cout << "Average number of OpRayleigh per event:   "
         << fRayleighCount / TotNbofEvents << G4endl;
  G4cout << "Average number of OpAbsorption per event: "
         << fOpAbsorption / TotNbofEvents << G4endl;
  G4cout << "\nSurface events (on +X surface, maximum one per photon) this run:"
         << G4endl;
  G4cout << "# of primary particles:      " << std::setw(8) << TotNbofEvents
         << G4endl;
  G4cout << "OpAbsorption before surface: " << std::setw(8)
         << fOpAbsorptionPrior << G4endl;
  G4cout << "Total # of surface events:   " << std::setw(8) << fTotalSurface
         << G4endl;
  if(fParticle->GetParticleName() == "opticalphoton")
  {
    G4cout << "Unaccounted for:             " << std::setw(8)
           << fTotalSurface + fOpAbsorptionPrior - TotNbofEvents << G4endl;
  }
  G4cout << "\nSurface events by process:" << G4endl;
  if(fBoundaryProcs[Transmission] > 0)
  {
    G4cout << "  Transmission:              " << std::setw(8)
           << fBoundaryProcs[Transmission] << G4endl;
  }
  if(fBoundaryProcs[FresnelRefraction] > 0)
  {
    G4cout << "  Fresnel refraction:        " << std::setw(8)
           << fBoundaryProcs[FresnelRefraction] << G4endl;
  }
  if(fBoundaryProcs[FresnelReflection] > 0)
  {
    G4cout << "  Fresnel reflection:        " << std::setw(8)
           << fBoundaryProcs[FresnelReflection] << G4endl;
  }
  if(fBoundaryProcs[TotalInternalReflection] > 0)
  {
    G4cout << "  Total internal reflection: " << std::setw(8)
           << fBoundaryProcs[TotalInternalReflection] << G4endl;
  }
  if(fBoundaryProcs[LambertianReflection] > 0)
  {
    G4cout << "  Lambertian reflection:     " << std::setw(8)
           << fBoundaryProcs[LambertianReflection] << G4endl;
  }
  if(fBoundaryProcs[LobeReflection] > 0)
  {
    G4cout << "  Lobe reflection:           " << std::setw(8)
           << fBoundaryProcs[LobeReflection] << G4endl;
  }
  if(fBoundaryProcs[SpikeReflection] > 0)
  {
    G4cout << "  Spike reflection:          " << std::setw(8)
           << fBoundaryProcs[SpikeReflection] << G4endl;
  }
  if(fBoundaryProcs[BackScattering] > 0)
  {
    G4cout << "  Backscattering:            " << std::setw(8)
           << fBoundaryProcs[BackScattering] << G4endl;
  }
  if(fBoundaryProcs[Absorption] > 0)
  {
    G4cout << "  Absorption:                " << std::setw(8)
           << fBoundaryProcs[Absorption] << G4endl;
  }
  if(fBoundaryProcs[Detection] > 0)
  {
    G4cout << "  Detection:                 " << std::setw(8)
           << fBoundaryProcs[Detection] << G4endl;
  }
  if(fBoundaryProcs[NotAtBoundary] > 0)
  {
    G4cout << "  Not at boundary:           " << std::setw(8)
           << fBoundaryProcs[NotAtBoundary] << G4endl;
  }
  if(fBoundaryProcs[SameMaterial] > 0)
  {
    G4cout << "  Same material:             " << std::setw(8)
           << fBoundaryProcs[SameMaterial] << G4endl;
  }
  if(fBoundaryProcs[StepTooSmall] > 0)
  {
    G4cout << "  Step too small:            " << std::setw(8)
           << fBoundaryProcs[StepTooSmall] << G4endl;
  }
  if(fBoundaryProcs[NoRINDEX] > 0)
  {
    G4cout << "  No RINDEX:                 " << std::setw(8)
           << fBoundaryProcs[NoRINDEX] << G4endl;
  }
  // LBNL polished
  if(fBoundaryProcs[PolishedLumirrorAirReflection] > 0)
  {
    G4cout << "  Polished Lumirror Air reflection: " << std::setw(8)
           << fBoundaryProcs[PolishedLumirrorAirReflection] << G4endl;
  }
  if(fBoundaryProcs[PolishedLumirrorGlueReflection] > 0)
  {
    G4cout << "  Polished Lumirror Glue reflection: " << std::setw(8)
           << fBoundaryProcs[PolishedLumirrorGlueReflection] << G4endl;
  }
  if(fBoundaryProcs[PolishedAirReflection] > 0)
  {
    G4cout << "  Polished Air reflection: " << std::setw(8)
           << fBoundaryProcs[PolishedAirReflection] << G4endl;
  }
  if(fBoundaryProcs[PolishedTeflonAirReflection] > 0)
  {
    G4cout << "  Polished Teflon Air reflection: " << std::setw(8)
           << fBoundaryProcs[PolishedTeflonAirReflection] << G4endl;
  }
  if(fBoundaryProcs[PolishedTiOAirReflection] > 0)
  {
    G4cout << "  Polished TiO Air reflection: " << std::setw(8)
           << fBoundaryProcs[PolishedTiOAirReflection] << G4endl;
  }
  if(fBoundaryProcs[PolishedTyvekAirReflection] > 0)
  {
    G4cout << "  Polished Tyvek Air reflection: " << std::setw(8)
           << fBoundaryProcs[PolishedTyvekAirReflection] << G4endl;
  }
  if(fBoundaryProcs[PolishedVM2000AirReflection] > 0)
  {
    G4cout << "  Polished VM2000 Air reflection: " << std::setw(8)
           << fBoundaryProcs[PolishedVM2000AirReflection] << G4endl;
  }
  if(fBoundaryProcs[PolishedVM2000GlueReflection] > 0)
  {
    G4cout << "  Polished VM2000 Glue reflection: " << std::setw(8)
           << fBoundaryProcs[PolishedVM2000GlueReflection] << G4endl;
  }
  // LBNL etched
  if(fBoundaryProcs[EtchedLumirrorAirReflection] > 0)
  {
    G4cout << "  Etched Lumirror Air reflection: " << std::setw(8)
           << fBoundaryProcs[EtchedLumirrorAirReflection] << G4endl;
  }
  if(fBoundaryProcs[EtchedLumirrorGlueReflection] > 0)
  {
    G4cout << "  Etched Lumirror Glue reflection: " << std::setw(8)
           << fBoundaryProcs[EtchedLumirrorGlueReflection] << G4endl;
  }
  if(fBoundaryProcs[EtchedAirReflection] > 0)
  {
    G4cout << "  Etched Air reflection: " << std::setw(8)
           << fBoundaryProcs[EtchedAirReflection] << G4endl;
  }
  if(fBoundaryProcs[EtchedTeflonAirReflection] > 0)
  {
    G4cout << "  Etched Teflon Air reflection: " << std::setw(8)
           << fBoundaryProcs[EtchedTeflonAirReflection] << G4endl;
  }
  if(fBoundaryProcs[EtchedTiOAirReflection] > 0)
  {
    G4cout << "  Etched TiO Air reflection: " << std::setw(8)
           << fBoundaryProcs[EtchedTiOAirReflection] << G4endl;
  }
  if(fBoundaryProcs[EtchedTyvekAirReflection] > 0)
  {
    G4cout << "  Etched Tyvek Air reflection: " << std::setw(8)
           << fBoundaryProcs[EtchedTyvekAirReflection] << G4endl;
  }
  if(fBoundaryProcs[EtchedVM2000AirReflection] > 0)
  {
    G4cout << "  Etched VM2000 Air reflection: " << std::setw(8)
           << fBoundaryProcs[EtchedVM2000AirReflection] << G4endl;
  }
  if(fBoundaryProcs[EtchedVM2000GlueReflection] > 0)
  {
    G4cout << "  Etched VM2000 Glue reflection: " << std::setw(8)
           << fBoundaryProcs[EtchedVM2000GlueReflection] << G4endl;
  }
  // LBNL ground
  if(fBoundaryProcs[GroundLumirrorAirReflection] > 0)
  {
    G4cout << "  Ground Lumirror Air reflection: " << std::setw(8)
           << fBoundaryProcs[GroundLumirrorAirReflection] << G4endl;
  }
  if(fBoundaryProcs[GroundLumirrorGlueReflection] > 0)
  {
    G4cout << "  Ground Lumirror Glue reflection: " << std::setw(8)
           << fBoundaryProcs[GroundLumirrorGlueReflection] << G4endl;
  }
  if(fBoundaryProcs[GroundAirReflection] > 0)
  {
    G4cout << "  Ground Air reflection: " << std::setw(8)
           << fBoundaryProcs[GroundAirReflection] << G4endl;
  }
  if(fBoundaryProcs[GroundTeflonAirReflection] > 0)
  {
    G4cout << "  Ground Teflon Air reflection: " << std::setw(8)
           << fBoundaryProcs[GroundTeflonAirReflection] << G4endl;
  }
  if(fBoundaryProcs[GroundTiOAirReflection] > 0)
  {
    G4cout << "  Ground TiO Air reflection: " << std::setw(8)
           << fBoundaryProcs[GroundTiOAirReflection] << G4endl;
  }
  if(fBoundaryProcs[GroundTyvekAirReflection] > 0)
  {
    G4cout << "  Ground Tyvek Air reflection: " << std::setw(8)
           << fBoundaryProcs[GroundTyvekAirReflection] << G4endl;
  }
  if(fBoundaryProcs[GroundVM2000AirReflection] > 0)
  {
    G4cout << "  Ground VM2000 Air reflection: " << std::setw(8)
           << fBoundaryProcs[GroundVM2000AirReflection] << G4endl;
  }
  if(fBoundaryProcs[GroundVM2000GlueReflection] > 0)
  {
    G4cout << "  Ground VM2000 Glue reflection: " << std::setw(8)
           << fBoundaryProcs[GroundVM2000GlueReflection] << G4endl;
  }
  if(fBoundaryProcs[CoatedDielectricRefraction] > 0)
  {
    G4cout << "  CoatedDielectricRefraction: " << std::setw(8)
           << fBoundaryProcs[CoatedDielectricRefraction] << G4endl;
  }
  if(fBoundaryProcs[CoatedDielectricReflection] > 0)
  {
    G4cout << "  CoatedDielectricReflection: " << std::setw(8)
           << fBoundaryProcs[CoatedDielectricReflection] << G4endl;
  }
  if(fBoundaryProcs[CoatedDielectricFrustratedTransmission] > 0)
  {
    G4cout << "  CoatedDielectricFrustratedTransmission: " << std::setw(8)
           << fBoundaryProcs[CoatedDielectricFrustratedTransmission] << G4endl;
  }

  G4int sum = std::accumulate(fBoundaryProcs.begin(), fBoundaryProcs.end(), 0);
  G4cout << " Sum:                        " << std::setw(8) << sum << G4endl;
  G4cout << " Unaccounted for:            " << std::setw(8)
         << fTotalSurface - sum << G4endl;

  G4cout << "---------------------------------\n";
  G4cout.setf(mode, std::ios::floatfield);
  G4cout.precision(prec);

  G4int histo_id_refract = analysisMan->GetH1Id("Fresnel refraction");
  G4int histo_id_reflect = analysisMan->GetH1Id("Fresnel reflection plus TIR");
  G4int histo_id_spike   = analysisMan->GetH1Id("Spike reflection");
  G4int histo_id_absorption = analysisMan->GetH1Id("Absorption");

  if(analysisMan->GetH1Activation(histo_id_refract) &&
     analysisMan->GetH1Activation(histo_id_reflect))
  {
    G4double rindex1 = det->GetTankMaterial()
                         ->GetMaterialPropertiesTable()
                         ->GetProperty(kRINDEX)
                         ->Value(fEkin);
    G4double rindex2 = det->GetWorldMaterial()
                         ->GetMaterialPropertiesTable()
                         ->GetProperty(kRINDEX)
                         ->Value(fEkin);

    auto histo_refract = analysisMan->GetH1(histo_id_refract);
    auto histo_reflect = analysisMan->GetH1(histo_id_reflect);
    // std::vector<G4double> refract;
    std::vector<G4double> reflect;
    // std::vector<G4double> tir;
    std::vector<G4double> tot;
    for(size_t i = 0; i < histo_refract->axis().bins(); ++i)
    {
      // refract.push_back(histo_refract->bin_height(i));
      reflect.push_back(histo_reflect->bin_height(i));
      // tir.push_back(histo_TIR->bin_height(i));
      tot.push_back(histo_refract->bin_height(i) +
                    histo_reflect->bin_height(i));
    }

    // find Brewster angle: Rp = 0
    //  need enough statistics for this method to work
    G4double min_angle = -1.;
    G4double min_val   = DBL_MAX;
    G4double bin_width = 0.;
    for(size_t i = 0; i < reflect.size(); ++i)
    {
      if(reflect[i] < min_val)
      {
        min_val   = reflect[i];
        min_angle = histo_reflect->axis().bin_lower_edge(i);
        bin_width = histo_reflect->axis().bin_upper_edge(i) -
                    histo_reflect->axis().bin_lower_edge(i);
        min_angle += bin_width / 2.;
      }
    }
    G4cout << "Polarization of primary optical photons: "
           << fPolarization / deg << " deg." << G4endl;
    if(fPolarized && fPolarization == 0.0)
    {
      G4cout << "Reflectance shows a minimum at: " << min_angle << " +/- "
             << bin_width / 2;
      G4cout << " deg. Expected Brewster angle: "
             << (360. / CLHEP::twopi) * std::atan(rindex2 / rindex1)
             << " deg. " << G4endl;
    }

    // find angle of total internal reflection:  T -> 0
    //   last bin for T > 0
    min_angle = -1.;
    min_val   = DBL_MAX;
    for(size_t i = 0; i < histo_refract->axis().bins() - 1; ++i)
    {
      if(histo_refract->bin_height(i) > 0. &&
         histo_refract->bin_height(i + 1) == 0.)
      {
        min_angle = histo_refract->axis().bin_lower_edge(i);
        bin_width = histo_reflect->axis().bin_upper_edge(i) -
                    histo_reflect->axis().bin_lower_edge(i);
        min_angle += bin_width / 2.;
        break;
      }
    }
    if(fPolarized)
    {
      G4cout << "Fresnel transmission goes to 0 at: " << min_angle << " +/- "
             << bin_width / 2. << " deg."
             << " Expected: "
             << (360. / CLHEP::twopi) * std::asin(rindex2 / rindex1)
             << " deg." << G4endl;
    }

    // Normalize the transmission/reflection histos so that max is 1.
    // Only if x values are the same
    if((analysisMan->GetH1Nbins(histo_id_refract) ==
        analysisMan->GetH1Nbins(histo_id_reflect)) &&
       (analysisMan->GetH1Xmin(histo_id_refract) ==
        analysisMan->GetH1Xmin(histo_id_reflect)) &&
       (analysisMan->GetH1Xmax(histo_id_refract) ==
        analysisMan->GetH1Xmax(histo_id_reflect)))
    {
      unsigned int ent;
      G4double sw;
      G4double sw2;
      G4double sx2;
      G4double sx2w;
      for(size_t bin = 0; bin < histo_refract->axis().bins(); ++bin)
      {
        // "bin+1" below because bin 0 is underflow bin
        // NB. We are ignoring underflow/overflow bins
        histo_refract->get_bin_content(bin + 1, ent, sw, sw2, sx2, sx2w);
        if(tot[bin] > 0)
        {
          sw /= tot[bin];
          // bin error is sqrt(sw2)
          sw2 /= (tot[bin] * tot[bin]);
          sx2 /= (tot[bin] * tot[bin]);
          sx2w /= (tot[bin] * tot[bin]);
          histo_refract->set_bin_content(bin + 1, ent, sw, sw2, sx2, sx2w);
        }

        histo_reflect->get_bin_content(bin + 1, ent, sw, sw2, sx2, sx2w);
        if(tot[bin] > 0)
        {
          sw /= tot[bin];
          // bin error is sqrt(sw2)
          sw2 /= (tot[bin] * tot[bin]);
          sx2 /= (tot[bin] * tot[bin]);
          sx2w /= (tot[bin] * tot[bin]);
          histo_reflect->set_bin_content(bin + 1, ent, sw, sw2, sx2, sx2w);
        }

        G4int histo_id_fresnelrefl =
          analysisMan->GetH1Id("Fresnel reflection");
        auto histo_fresnelreflect = analysisMan->GetH1(histo_id_fresnelrefl);
        histo_fresnelreflect->get_bin_content(bin + 1, ent, sw, sw2, sx2,
                                              sx2w);
        if(tot[bin] > 0)
        {
          sw /= tot[bin];
          // bin error is sqrt(sw2)
          sw2 /= (tot[bin] * tot[bin]);
          sx2 /= (tot[bin] * tot[bin]);
          sx2w /= (tot[bin] * tot[bin]);
          histo_fresnelreflect->set_bin_content(bin + 1, ent, sw, sw2, sx2,
                                                sx2w);
        }

        G4int histo_id_TIR =
          analysisMan->GetH1Id("Total internal reflection");
        auto histo_TIR = analysisMan->GetH1(histo_id_TIR);
        if(analysisMan->GetH1Activation(histo_id_TIR))
        {
          histo_TIR->get_bin_content(bin + 1, ent, sw, sw2, sx2, sx2w);
          if(tot[bin] > 0)
          {
            sw /= tot[bin];
            // bin error is sqrt(sw2)
            sw2 /= (tot[bin] * tot[bin]);
            sx2 /= (tot[bin] * tot[bin]);
            sx2w /= (tot[bin] * tot[bin]);
            histo_TIR->set_bin_content(bin + 1, ent, sw, sw2, sx2, sx2w);
          }
        }
      }
    }
    else
    {
      G4cout << "Not going to normalize refraction and reflection "
             << "histograms because bins are not the same." << G4endl;
    }
  }

  // complex index of refraction; have spike reflection and absorption
  // Only works for polished surfaces. Ground surfaces neglected.
  else if(analysisMan->GetH1Activation(histo_id_absorption) &&
          analysisMan->GetH1Activation(histo_id_spike))
  {
    auto histo_spike      = analysisMan->GetH1(histo_id_spike);
    auto histo_absorption = analysisMan->GetH1(histo_id_absorption);

    std::vector<G4double> tot;
    for(size_t i = 0; i < histo_absorption->axis().bins(); ++i)
    {
      tot.push_back(histo_absorption->bin_height(i) +
                    histo_spike->bin_height(i));
    }

    if((analysisMan->GetH1Nbins(histo_id_absorption) ==
        analysisMan->GetH1Nbins(histo_id_spike)) &&
       (analysisMan->GetH1Xmin(histo_id_absorption) ==
        analysisMan->GetH1Xmin(histo_id_spike)) &&
       (analysisMan->GetH1Xmax(histo_id_absorption) ==
        analysisMan->GetH1Xmax(histo_id_spike)))
    {
      unsigned int ent;
      G4double sw;
      G4double sw2;
      G4double sx2;
      G4double sx2w;
      for(size_t bin = 0; bin < histo_absorption->axis().bins(); ++bin)
      {
        // "bin+1" below because bin 0 is underflow bin
        // NB. We are ignoring underflow/overflow bins
        histo_absorption->get_bin_content(bin + 1, ent, sw, sw2, sx2, sx2w);
        if(tot[bin] > 0)
        {
          sw /= tot[bin];
          // bin error is sqrt(sw2)
          sw2 /= (tot[bin] * tot[bin]);
          sx2 /= (tot[bin] * tot[bin]);
          sx2w /= (tot[bin] * tot[bin]);
          histo_absorption->set_bin_content(bin + 1, ent, sw, sw2, sx2, sx2w);
        }

        histo_spike->get_bin_content(bin + 1, ent, sw, sw2, sx2, sx2w);
        if(tot[bin] > 0)
        {
          sw /= tot[bin];
          // bin error is sqrt(sw2)
          sw2 /= (tot[bin] * tot[bin]);
          sx2 /= (tot[bin] * tot[bin]);
          sx2w /= (tot[bin] * tot[bin]);
          histo_spike->set_bin_content(bin + 1, ent, sw, sw2, sx2, sx2w);
        }
      }
    }
    else
    {
      G4cout << "Not going to normalize spike reflection and absorption "
             << "histograms because bins are not the same." << G4endl;
    }
  }
}
