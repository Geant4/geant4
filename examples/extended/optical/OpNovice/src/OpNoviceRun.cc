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
/// \file OpNovice/src/OpNoviceRun.cc
/// \brief Implementation of the OpNoviceRun class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "OpNoviceRun.hh"

#include "G4ParticleDefinition.hh"
#include "G4Run.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceRun::OpNoviceRun()
  : G4Run()
{
  fParticle             = nullptr;
  fEnergy               = -1.;
  fCerenkovCounter      = 0.;
  fCerenkov2            = 0.;
  fScintillationCounter = 0.;
  fScintillation2       = 0.;
  fRayleighCounter      = 0.;
  fRayleigh2            = 0.;
  fAbsorptionCounter    = 0.;
  fAbsorption2          = 0.;
  fMieCounter           = 0.;
  fMie2                 = 0.;
  fBoundaryCounter      = 0.;
  fBoundary2            = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceRun::~OpNoviceRun() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceRun::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{
  fParticle = particle;
  fEnergy   = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceRun::Merge(const G4Run* run)
{
  const OpNoviceRun* localRun = static_cast<const OpNoviceRun*>(run);

  fParticle = localRun->fParticle;
  fEnergy   = localRun->fEnergy;

  fCerenkovCounter += localRun->fCerenkovCounter;
  fCerenkov2 += localRun->fCerenkov2;
  fScintillationCounter += localRun->fScintillationCounter;
  fScintillation2 += localRun->fScintillation2;

  fRayleighCounter += localRun->fRayleighCounter;
  fRayleigh2 += localRun->fRayleigh2;
  fAbsorptionCounter += localRun->fAbsorptionCounter;
  fAbsorption2 += localRun->fAbsorption2;
  fMieCounter += localRun->fMieCounter;
  fMie2 += localRun->fMie2;
  fBoundaryCounter += localRun->fBoundaryCounter;
  fBoundary2 += localRun->fBoundary2;

  G4Run::Merge(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void OpNoviceRun::EndOfRun()
{
  if(numberOfEvent == 0)
    return;
  G4double TotNbofEvents = G4double(numberOfEvent);

  fCerenkovCounter /= TotNbofEvents;
  fCerenkov2 /= TotNbofEvents;
  G4double rmsCerenkov = fCerenkov2 - fCerenkovCounter * fCerenkovCounter;
  if(rmsCerenkov > 0.)
    rmsCerenkov = std::sqrt(rmsCerenkov);
  else
    rmsCerenkov = 0.;

  fScintillationCounter /= TotNbofEvents;
  fScintillation2 /= TotNbofEvents;
  G4double rmsScint =
    fScintillation2 - fScintillationCounter * fScintillationCounter;
  if(rmsScint > 0.)
    rmsScint = std::sqrt(rmsScint);
  else
    rmsScint = 0.;

  fRayleighCounter /= TotNbofEvents;
  fRayleigh2 /= TotNbofEvents;
  G4double rmsRayleigh = fRayleigh2 - fRayleighCounter * fRayleighCounter;
  if(rmsRayleigh > 0.)
    rmsRayleigh = std::sqrt(rmsRayleigh);
  else
    rmsRayleigh = 0.;

  fAbsorptionCounter /= TotNbofEvents;
  fAbsorption2 /= TotNbofEvents;
  G4double rmsAbsorption =
    fAbsorption2 - fAbsorptionCounter * fAbsorptionCounter;
  if(rmsAbsorption > 0.)
    rmsAbsorption = std::sqrt(rmsAbsorption);
  else
    rmsAbsorption = 0.;

  fMieCounter /= TotNbofEvents;
  fMie2 /= TotNbofEvents;
  G4double rmsMie = fMie2 - fMieCounter * fMieCounter;
  if(rmsMie > 0.)
    rmsMie = std::sqrt(rmsMie);
  else
    rmsMie = 0.;

  fBoundaryCounter /= TotNbofEvents;
  fBoundary2 /= TotNbofEvents;
  G4double rmsBoundary = fBoundary2 - fBoundaryCounter * fBoundaryCounter;
  if(rmsBoundary > 0.)
    rmsBoundary = std::sqrt(rmsBoundary);
  else
    rmsBoundary = 0.;

  G4int prec = G4cout.precision(3);
  G4cout << "\n ======================== run summary ======================\n";

  G4cout << "Primary particle was: " << fParticle->GetParticleName()
         << " with energy " << G4BestUnit(fEnergy, "Energy") << "." << G4endl;
  G4cout << "Number of events: " << numberOfEvent << G4endl;

  G4cout << "Average number of Cerenkov photons created per event: "
         << fCerenkovCounter << " +- " << rmsCerenkov << G4endl;
  G4cout << "Average number of scintillation photons created per event: "
         << fScintillationCounter << " +- " << rmsScint << G4endl;

  G4cout << "Average number of optical Rayleigh interactions per event: "
         << fRayleighCounter << " +- " << rmsRayleigh << G4endl;
  G4cout << "Average number of optical absorption interactions per event: "
         << fAbsorptionCounter << " +- " << rmsAbsorption << G4endl;
  G4cout << "Average number of optical Mie interactions per event: "
         << fMieCounter << " +- " << rmsMie << G4endl;
  G4cout << "Average number of optical boundary interactions per event: "
         << fBoundaryCounter << " +- " << rmsBoundary << G4endl;

  G4cout << G4endl;
  G4cout.precision(prec);
}
