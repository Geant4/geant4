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
/// \file optical/OpNovice2/include/RunAction.hh
/// \brief Definition of the RunAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "G4OpBoundaryProcess.hh"
#include "G4Run.hh"

class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class Run : public G4Run
{
 public:
  Run();
  ~Run() override = default;

  void SetPrimary(G4ParticleDefinition* particle, G4double energy,
                  G4bool polarized, G4double polarization);

  //  particle energy
  void AddCerenkovEnergy(G4double en) { fCerenkovEnergy += en; }
  void AddScintillationEnergy(G4double en) { fScintEnergy += en; }
  void AddWLSAbsorptionEnergy(G4double en) { fWLSAbsorptionEnergy += en; }
  void AddWLSEmissionEnergy(G4double en) { fWLSEmissionEnergy += en; }
  void AddWLS2AbsorptionEnergy(G4double en) { fWLS2AbsorptionEnergy += en; }
  void AddWLS2EmissionEnergy(G4double en) { fWLS2EmissionEnergy += en; }

  // number of particles
  void AddCerenkov() { fCerenkovCount += 1; }
  void AddScintillation() { fScintCount += 1; }
  void AddRayleigh() { fRayleighCount += 1; }
  void AddWLSAbsorption() { fWLSAbsorptionCount += 1; }
  void AddWLSEmission() { fWLSEmissionCount += 1; }
  void AddWLS2Absorption() { fWLS2AbsorptionCount += 1; }
  void AddWLS2Emission() { fWLS2EmissionCount += 1; }

  void AddOpAbsorption() { fOpAbsorption += 1; }
  void AddOpAbsorptionPrior() { fOpAbsorptionPrior += 1; }

  void AddFresnelRefraction() { fBoundaryProcs[FresnelRefraction] += 1; }
  void AddFresnelReflection() { fBoundaryProcs[FresnelReflection] += 1; }
  void AddTransmission() { fBoundaryProcs[Transmission] += 1; }
  void AddTotalInternalReflection()
  {
    fBoundaryProcs[TotalInternalReflection] += 1;
  }
  void AddLambertianReflection()
  {
    fBoundaryProcs[LambertianReflection] += 1;
  }
  void AddLobeReflection() { fBoundaryProcs[LobeReflection] += 1; }
  void AddSpikeReflection() { fBoundaryProcs[SpikeReflection] += 1; }
  void AddBackScattering() { fBoundaryProcs[BackScattering] += 1; }
  void AddAbsorption() { fBoundaryProcs[Absorption] += 1; }
  void AddDetection() { fBoundaryProcs[Detection] += 1; }
  void AddNotAtBoundary() { fBoundaryProcs[NotAtBoundary] += 1; }
  void AddSameMaterial() { fBoundaryProcs[SameMaterial] += 1; }
  void AddStepTooSmall() { fBoundaryProcs[StepTooSmall] += 1; }
  void AddNoRINDEX() { fBoundaryProcs[NoRINDEX] += 1; }

  void AddTotalSurface() { fTotalSurface += 1; }
  void AddPolishedLumirrorAirReflection()
  {
    fBoundaryProcs[PolishedLumirrorAirReflection] += 1;
  }
  void AddPolishedLumirrorGlueReflection()
  {
    fBoundaryProcs[PolishedLumirrorGlueReflection] += 1;
  }
  void AddPolishedAirReflection()
  {
    fBoundaryProcs[PolishedAirReflection] += 1;
  }
  void AddPolishedTeflonAirReflection()
  {
    fBoundaryProcs[PolishedTeflonAirReflection] += 1;
  }
  void AddPolishedTiOAirReflection()
  {
    fBoundaryProcs[PolishedTiOAirReflection] += 1;
  }
  void AddPolishedTyvekAirReflection()
  {
    fBoundaryProcs[PolishedTyvekAirReflection] += 1;
  }
  void AddPolishedVM2000AirReflection()
  {
    fBoundaryProcs[PolishedVM2000AirReflection] += 1;
  }
  void AddPolishedVM2000GlueReflection()
  {
    fBoundaryProcs[PolishedVM2000GlueReflection] += 1;
  }

  void AddEtchedLumirrorAirReflection()
  {
    fBoundaryProcs[EtchedLumirrorAirReflection] += 1;
  }
  void AddEtchedLumirrorGlueReflection()
  {
    fBoundaryProcs[EtchedLumirrorGlueReflection] += 1;
  }
  void AddEtchedAirReflection()
  {
    fBoundaryProcs[EtchedAirReflection] += 1;
  }
  void AddEtchedTeflonAirReflection()
  {
    fBoundaryProcs[EtchedTeflonAirReflection] += 1;
  }
  void AddEtchedTiOAirReflection()
  {
    fBoundaryProcs[EtchedTiOAirReflection] += 1;
  }
  void AddEtchedTyvekAirReflection()
  {
    fBoundaryProcs[EtchedTyvekAirReflection] += 1;
  }
  void AddEtchedVM2000AirReflection()
  {
    fBoundaryProcs[EtchedVM2000AirReflection] += 1;
  }
  void AddEtchedVM2000GlueReflection()
  {
    fBoundaryProcs[EtchedVM2000GlueReflection] += 1;
  }

  void AddGroundLumirrorAirReflection()
  {
    fBoundaryProcs[GroundLumirrorAirReflection] += 1;
  }
  void AddGroundLumirrorGlueReflection()
  {
    fBoundaryProcs[GroundLumirrorGlueReflection] += 1;
  }
  void AddGroundAirReflection()
  {
    fBoundaryProcs[GroundAirReflection] += 1;
  }
  void AddGroundTeflonAirReflection()
  {
    fBoundaryProcs[GroundTeflonAirReflection] += 1;
  }
  void AddGroundTiOAirReflection()
  {
    fBoundaryProcs[GroundTiOAirReflection] += 1;
  }
  void AddGroundTyvekAirReflection()
  {
    fBoundaryProcs[GroundTyvekAirReflection] += 1;
  }
  void AddGroundVM2000AirReflection()
  {
    fBoundaryProcs[GroundVM2000AirReflection] += 1;
  }
  void AddGroundVM2000GlueReflection()
  {
    fBoundaryProcs[GroundVM2000GlueReflection] += 1;
  }

  void AddDichroic() { fBoundaryProcs[Dichroic] += 1; }
  void AddCoatedDielectricRefraction()
  {
    fBoundaryProcs[CoatedDielectricRefraction] += 1;
  }
  void AddCoatedDielectricReflection()
  {
    fBoundaryProcs[CoatedDielectricReflection] += 1;
  }
  void AddCoatedDielectricFrustratedTransmission()
  {
    fBoundaryProcs[CoatedDielectricFrustratedTransmission] += 1;
  }

  void Merge(const G4Run*) override;

  void EndOfRun();

 private:
  // primary particle
  G4ParticleDefinition* fParticle = nullptr;
  G4double fEkin = -1.;
  G4bool fPolarized = false;
  G4double fPolarization = 0.;

  G4double fCerenkovEnergy = 0.;
  G4double fScintEnergy = 0.;
  G4double fWLSAbsorptionEnergy = 0.;
  G4double fWLSEmissionEnergy = 0.;
  G4double fWLS2AbsorptionEnergy = 0.;
  G4double fWLS2EmissionEnergy = 0.;

  // number of particles
  G4int fCerenkovCount = 0;
  G4int fScintCount = 0;
  G4int fWLSAbsorptionCount = 0;
  G4int fWLSEmissionCount = 0;
  G4int fWLS2AbsorptionCount = 0;
  G4int fWLS2EmissionCount = 0;
  // number of events
  G4int fRayleighCount = 0;

  // non-boundary processes
  G4int fOpAbsorption = 0;

  // prior to boundary:
  G4int fOpAbsorptionPrior = 0;

  // boundary proc
  std::vector<G4int> fBoundaryProcs;

  G4int fTotalSurface = 0;
};

#endif /* Run_h */
