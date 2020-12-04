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
//
//---------------------------------------------------------------------------
//
// ClassName:   G4OpticalPhysics
//
// Author:      P.Gumplinger 30.09.2009
//
// Modified:    P.Gumplinger 29.09.2011
//              (based on code from I. Hrivnacova)
//
//----------------------------------------------------------------------------
//

#include "G4OpticalPhysics.hh"
#include "G4OpticalParameters.hh"

#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpWLS.hh"
#include "G4OpWLS2.hh"
#include "G4Scintillation.hh"
#include "G4Cerenkov.hh"

#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
G4_DECLARE_PHYSCONSTR_FACTORY(G4OpticalPhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4OpticalPhysics::G4OpticalPhysics(G4int verbose, const G4String& name)
  : G4VPhysicsConstructor(name)
{
  verboseLevel = verbose;
  G4OpticalParameters::Instance()->SetVerboseLevel(verboseLevel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4OpticalPhysics::~G4OpticalPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4OpticalPhysics::PrintStatistics() const
{
  G4OpticalParameters::Instance()->Dump();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4OpticalPhysics::ConstructParticle()
{
  G4OpticalPhoton::OpticalPhotonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4OpticalPhysics::ConstructProcess()
{
  if(verboseLevel>0)
         G4cout <<"G4OpticalPhysics:: Add Optical Physics Processes"<< G4endl;

  auto params = G4OpticalParameters::Instance();

  // Add Optical Processes

  G4ProcessManager* pManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
  if (!pManager) {
     G4ExceptionDescription ed;
     ed << "Optical Photon without a Process Manager";
     G4Exception("G4OpticalPhysics::ConstructProcess()","",
                  FatalException,ed);
     return;
  }

  G4OpAbsorption* absorption  = new G4OpAbsorption();
  if (params->GetProcessActivation("OpAbsorption")) pManager->AddDiscreteProcess(absorption);

  G4OpRayleigh* rayleigh = new G4OpRayleigh();
  if (params->GetProcessActivation("OpRayleigh")) pManager->AddDiscreteProcess(rayleigh);

  G4OpMieHG* mie = new G4OpMieHG();
  if (params->GetProcessActivation("OpMieHG")) pManager->AddDiscreteProcess(mie);

  G4OpBoundaryProcess* boundary = new G4OpBoundaryProcess();
  if (params->GetProcessActivation("OpBoundary")) pManager->AddDiscreteProcess(boundary);

  G4OpWLS* wls = new G4OpWLS();
  if (params->GetProcessActivation("OpWLS")) pManager->AddDiscreteProcess(wls);

  G4OpWLS2* wls2 = new G4OpWLS2();
  if (params->GetProcessActivation("OpWLS2")) pManager->AddDiscreteProcess(wls2);

  G4Scintillation* scint = new G4Scintillation();
  G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
  scint->AddSaturation(emSaturation);

  G4Cerenkov* cerenkov = new G4Cerenkov();

  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();

  while( (*myParticleIterator)() ){

    G4ParticleDefinition* particle = myParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    pManager = particle->GetProcessManager();
    if (!pManager) {
       G4ExceptionDescription ed;
       ed << "Particle " << particleName << "without a Process Manager";
       G4Exception("G4OpticalPhysics::ConstructProcess()","",
                    FatalException, ed);
       return;                 // else coverity complains for pManager use below
    }

    if (cerenkov->IsApplicable(*particle) && params->GetProcessActivation("Cerenkov")) {
          pManager->AddProcess(cerenkov);
          pManager->SetProcessOrdering(cerenkov,idxPostStep);
    }
    if (scint->IsApplicable(*particle) && params->GetProcessActivation("Scintillation")) {
          pManager->AddProcess(scint);
          pManager->SetProcessOrderingToLast(scint,idxAtRest);
          pManager->SetProcessOrderingToLast(scint,idxPostStep);
    }
    if (boundary->IsApplicable(*particle) && params->GetProcessActivation("OpBoundary")) {
          pManager->SetProcessOrderingToLast(boundary,idxPostStep);
    }
  }


  if (verboseLevel > 1) PrintStatistics();
  if (verboseLevel > 0)
    G4cout << "### " << namePhysics << " physics constructed." << G4endl;
}

// DEPRECATED
// the methods below are kept for backwards compatibility, and are to be
// removed in future. Please use the methods in
// G4OpticalParameters instead.

void G4OpticalPhysics::Configure(G4OpticalProcessIndex index, G4bool val) {
	G4OpticalParameters* params = G4OpticalParameters::Instance();
  if      (index == kCerenkov)      params->SetProcessActivation("Cerenkov", val);
  else if (index == kScintillation) params->SetProcessActivation("Scintillation", val);
  else if (index == kAbsorption)    params->SetProcessActivation("Absorption", val);
  else if (index == kRayleigh)      params->SetProcessActivation("Rayleigh", val);
  else if (index == kMieHG)         params->SetProcessActivation("MieHG", val);
  else if (index == kBoundary)      params->SetProcessActivation("Boundary", val);
  else if (index == kWLS)           params->SetProcessActivation("WLS", val);
  else if (index == kWLS2)          params->SetProcessActivation("WLS2", val);
}

void G4OpticalPhysics::SetTrackSecondariesFirst(G4OpticalProcessIndex index, G4bool val) {
  if (index == kCerenkov)
    G4OpticalParameters::Instance()->SetCerenkovTrackSecondariesFirst(val);
  else if (index == kScintillation)
    G4OpticalParameters::Instance()->SetScintTrackSecondariesFirst(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetTrackSecondariesFirst is deprecated." << G4endl
     << "Use G4OpticalParameters::Set[Cerenkov/Scint]TrackSecondariesFirst(G4bool) instead.";
  PrintWarning(ed);
}

// Cerenkov
void G4OpticalPhysics::SetMaxNumPhotonsPerStep(G4int val) {
  G4OpticalParameters::Instance()->SetCerenkovMaxPhotonsPerStep(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetMaxNumPhotonsPerStep is deprecated." << G4endl
     << "Use G4OpticalParameters::SetCerenkovMaxPhotonsPerStep(G4int) instead.";
  PrintWarning(ed);
}

void G4OpticalPhysics::SetMaxBetaChangePerStep(G4double val) {
  G4OpticalParameters::Instance()->SetCerenkovMaxBetaChange(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetMaxBetaChangePerStep is deprecated." << G4endl
     << "Use G4OpticalParameters::SetCerenkovMaxBetaChange(G4double) instead.";
  PrintWarning(ed);
}

void G4OpticalPhysics::SetCerenkovStackPhotons(G4bool val) {
  G4OpticalParameters::Instance()->SetCerenkovStackPhotons(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetCerenkovStackPhotons is deprecated." << G4endl
     << "Use G4OpticalParameters::SetCerenkovStackPhotons(G4int) "
     << "instead.";
  PrintWarning(ed);
}

void G4OpticalPhysics::SetCerenkovTrackSecondariesFirst(G4bool val) {
  G4OpticalParameters::Instance()->SetCerenkovTrackSecondariesFirst(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetCerenkovTrackSecondariesFirst is deprecated." << G4endl
     << "Use G4OpticalParameters::SetCerenkovTrackSecondariesFirst(G4int) instead.";
  PrintWarning(ed);
}

void G4OpticalPhysics::SetCerenkovVerbosity(G4int val) {
  G4OpticalParameters::Instance()->SetCerenkovVerboseLevel(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetCerenkovVerbosity is deprecated." << G4endl
     << "Use G4OpticalParameters::SetCerenkovVerbosity(G4int) instead.";
  PrintWarning(ed);
}

// Scintillation
void G4OpticalPhysics::SetScintillationYieldFactor(G4double val) {
  G4OpticalParameters::Instance()->SetScintYieldFactor(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetScintillationYieldFactor is deprecated." << G4endl
     << "Use G4OpticalParameters::SetScintYieldFactor(G4double) instead.";
  PrintWarning(ed);
}

void G4OpticalPhysics::SetScintillationExcitationRatio(G4double val) {
  G4OpticalParameters::Instance()->SetScintExcitationRatio(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetScintillationExcitationRatio is deprecated." << G4endl
     << "Use G4OpticalParameters::SetScintExcitationRatio(G4double) instead.";
  PrintWarning(ed);
}

void G4OpticalPhysics::SetScintillationByParticleType(G4bool val) {
  G4OpticalParameters::Instance()->SetScintByParticleType(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetScintillationByParticleType is deprecated." << G4endl
     << "Use G4OpticalParameters::SetScintByParticleType(G4bool) instead.";
  PrintWarning(ed);
}

void G4OpticalPhysics::SetScintillationTrackInfo(G4bool val) {
  G4OpticalParameters::Instance()->SetScintTrackInfo(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetScintillationTrackInfo is deprecated." << G4endl
     << "Use G4OpticalParameters::SetScintTrackInfo(G4bool) instead.";
  PrintWarning(ed);
}

void G4OpticalPhysics::SetScintillationTrackSecondariesFirst(G4bool val) {
  G4OpticalParameters::Instance()->SetScintTrackSecondariesFirst(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetScintillationTrackSecondariesFirst is deprecated." << G4endl
     << "Use G4OpticalParameters::SetScintTrackSecondariesFirst(G4bool) instead.";
  PrintWarning(ed);
}

void G4OpticalPhysics::SetFiniteRiseTime(G4bool val) {
  G4OpticalParameters::Instance()->SetScintFiniteRiseTime(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetFiniteRiseTime is deprecated." << G4endl
     << "Use G4OpticalParameters::SetScintFiniteRiseTime(G4bool) instead.";
  PrintWarning(ed);
}

void G4OpticalPhysics::SetScintillationStackPhotons(G4bool val) {
  G4OpticalParameters::Instance()->SetScintStackPhotons(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetScintillationStackPhotons is deprecated." << G4endl
     << "Use G4OpticalParameters::SetScintStackPhotons(G4bool) instead.";
  PrintWarning(ed);
}

void G4OpticalPhysics::SetScintillationVerbosity(G4int val) {
  G4OpticalParameters::Instance()->SetScintVerboseLevel(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetScintillationVerbosity is deprecated." << G4endl
     << "Use G4OpticalParameters::SetScintVerboseLevel(G4int) instead.";
  PrintWarning(ed);
}

void G4OpticalPhysics::SetScintillationEnhancedTimeConstants(G4bool val) {
  G4OpticalParameters::Instance()->SetScintEnhancedTimeConstants(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetScintillationEnhanceTimeConstants is deprecated." << G4endl
     << "Use G4OpticalParameters::SetScintEnhancedTimeConstants(G4bool) instead.";
  PrintWarning(ed);
}

//void AddScintillationSaturation(G4EmSaturation* );

// WLS
void G4OpticalPhysics::SetWLSTimeProfile(G4String val) {
  G4OpticalParameters::Instance()->SetWLSTimeProfile(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetWLSTimeProfile is deprecated." << G4endl
     << "Use G4OpticalParameters::SetWLSTimeProfile(G4String) instead.";
  PrintWarning(ed);
}

void G4OpticalPhysics::SetWLSVerbosity(G4int val) {
  G4OpticalParameters::Instance()->SetWLSVerboseLevel(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetWLSVerbosity is deprecated." << G4endl
     << "Use G4OpticalParameters::SetWLSVerboseLevel(G4int) instead.";
  PrintWarning(ed);
}

//boundary
void G4OpticalPhysics::SetBoundaryVerbosity(G4int val) {
  G4OpticalParameters::Instance()->SetBoundaryVerboseLevel(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetBoundaryVerbosity is deprecated." << G4endl
     << "Use G4OpticalParameters::SetBoundaryVerboseLevel(G4int) instead.";
  PrintWarning(ed);
}

void G4OpticalPhysics::SetInvokeSD(G4bool val) {
  G4OpticalParameters::Instance()->SetBoundaryInvokeSD(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetInvokeSD is deprecated." << G4endl
     << "Use G4OpticalParameters::SetBoundaryInvokeSD(G4bool) instead.";
  PrintWarning(ed);
}

void G4OpticalPhysics::SetAbsorptionVerbosity(G4int val) {
  G4OpticalParameters::Instance()->SetAbsorptionVerboseLevel(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetAbsorptionVerbosity is deprecated." << G4endl
     << "Use G4OpticalParameters::SetAbsorptionVerboseLevel(G4int) instead.";
  PrintWarning(ed);
}

void G4OpticalPhysics::SetRayleighVerbosity(G4int val) {
  G4OpticalParameters::Instance()->SetRayleighVerboseLevel(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetRayleighVerbosity is deprecated." << G4endl
     << "Use G4OpticalParameters::SetRayleighVerboseLevel(G4int) instead.";
  PrintWarning(ed);
}

void G4OpticalPhysics::SetMieVerbosity(G4int val) {
  G4OpticalParameters::Instance()->SetMieVerboseLevel(val);
  G4ExceptionDescription ed;
  ed << "Method G4OpticalPhysics::SetMieVerbosity is deprecated." << G4endl
     << "Use G4OpticalParameters::SetMieVerboseLevel(G4int) instead.";
  PrintWarning(ed);
}

void G4OpticalPhysics::PrintWarning(G4ExceptionDescription& ed) const {
    G4Exception("G4OpticalPhysics", "Optical0021", JustWarning, ed);
}

