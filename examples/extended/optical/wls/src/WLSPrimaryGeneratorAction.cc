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
/// \file optical/wls/src/WLSPrimaryGeneratorAction.cc
/// \brief Implementation of the WLSPrimaryGeneratorAction class
//
//

#include "WLSPrimaryGeneratorAction.hh"

#include "WLSDetectorConstruction.hh"
#include "WLSPrimaryGeneratorMessenger.hh"

#include "G4AutoLock.hh"
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalPhoton.hh"
#include "G4PhysicsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

namespace
{
  G4Mutex gen_mutex = G4MUTEX_INITIALIZER;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool WLSPrimaryGeneratorAction::fFirst = false;

WLSPrimaryGeneratorAction::WLSPrimaryGeneratorAction(
  WLSDetectorConstruction* dc)
{
  fDetector      = dc;
  fIntegralTable = nullptr;

  fParticleGun = new G4GeneralParticleSource();

  fGunMessenger = new WLSPrimaryGeneratorMessenger(this);

  fTimeConstant     = 0.;
  fUseSampledEnergy = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSPrimaryGeneratorAction::~WLSPrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;
  if(fIntegralTable)
  {
    fIntegralTable->clearAndDestroy();
    delete fIntegralTable;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPrimaryGeneratorAction::SetDecayTimeConstant(G4double time)
{
  fTimeConstant = time;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPrimaryGeneratorAction::BuildEmissionSpectrum()
{
  if(fIntegralTable)
    return;

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  G4int numOfMaterials                    = G4Material::GetNumberOfMaterials();

  if(!fIntegralTable)
    fIntegralTable = new G4PhysicsTable(numOfMaterials);

  for(G4int i = 0; i < numOfMaterials; ++i)
  {
    G4PhysicsFreeVector* vec = new G4PhysicsFreeVector();

    G4MaterialPropertiesTable* MPT =
      (*theMaterialTable)[i]->GetMaterialPropertiesTable();

    if(MPT)
    {
      G4MaterialPropertyVector* theWLSVector = MPT->GetProperty("WLSCOMPONENT");

      if(theWLSVector)
      {
        G4double currentIN = (*theWLSVector)[0];
        if(currentIN >= 0.0)
        {
          G4double currentPM  = theWLSVector->Energy(0);
          G4double currentCII = 0.0;
          vec->InsertValues(currentPM, currentCII);
          G4double prevPM  = currentPM;
          G4double prevCII = currentCII;
          G4double prevIN  = currentIN;

          for(size_t j = 1; j < theWLSVector->GetVectorLength(); ++j)
          {
            currentPM  = theWLSVector->Energy(j);
            currentIN  = (*theWLSVector)[j];
            currentCII = 0.5 * (prevIN + currentIN);
            currentCII = prevCII + (currentPM - prevPM) * currentCII;
            vec->InsertValues(currentPM, currentCII);
            prevPM  = currentPM;
            prevCII = currentCII;
            prevIN  = currentIN;
          }
        }
      }
    }
    fIntegralTable->insertAt(i, vec);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if(!fFirst)
  {
    fFirst = true;
    BuildEmissionSpectrum();
  }

  if(fUseSampledEnergy)
  {
    const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

    G4double sampledEnergy = 3. * eV;

    for(size_t j = 0; j < theMaterialTable->size(); ++j)
    {
      G4Material* fMaterial = (*theMaterialTable)[j];
      if(fMaterial->GetName() == "PMMA")
      {
        auto WLSIntensity =
          fMaterial->GetMaterialPropertiesTable()->GetProperty("WLSCOMPONENT");

        if(WLSIntensity)
        {
          G4PhysicsFreeVector* WLSIntegral =
            (G4PhysicsFreeVector*) ((*fIntegralTable)(
              fMaterial->GetIndex()));

          G4double CIImax   = WLSIntegral->GetMaxValue();
          G4double CIIvalue = G4UniformRand() * CIImax;

          sampledEnergy = WLSIntegral->GetEnergy(CIIvalue);
        }
      }
    }

    // this does not work.
    G4String cmd = "/gun/energy " + G4UIcommand::ConvertToString(sampledEnergy / eV) + " eV";
    G4UImanager::GetUIpointer()->ApplyCommand(cmd);
  }

  // The code behind this line is not thread safe because polarization
  // and time are randomly selected and GPS properties are global

  G4AutoLock l(&gen_mutex);
  if(fParticleGun->GetParticleDefinition() == G4OpticalPhoton::Definition())
  {
    SetOptPhotonPolar();
    SetOptPhotonTime();
  }

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPrimaryGeneratorAction::SetOptPhotonPolar()
{
  G4double angle = G4UniformRand() * 360.0 * deg;
  SetOptPhotonPolar(angle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
  if(fParticleGun->GetParticleDefinition()->GetParticleName() !=
     "opticalphoton")
  {
    G4cout << "-> warning from WLSPrimaryGeneratorAction::SetOptPhotonPolar()"
           << ":  the ParticleGun is not an opticalphoton" << G4endl;
    return;
  }

  G4ThreeVector normal(1., 0., 0.);
  G4ThreeVector kphoton = fParticleGun->GetParticleMomentumDirection();
  G4ThreeVector product = normal.cross(kphoton);
  G4double modul2       = product * product;

  G4ThreeVector e_perpend(0., 0., 1.);
  if(modul2 > 0.)
    e_perpend = (1. / std::sqrt(modul2)) * product;
  G4ThreeVector e_paralle = e_perpend.cross(kphoton);

  G4ThreeVector polar =
    std::cos(angle) * e_paralle + std::sin(angle) * e_perpend;
  fParticleGun->SetParticlePolarization(polar);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPrimaryGeneratorAction::SetOptPhotonTime()
{
  G4double time = -std::log(G4UniformRand()) * fTimeConstant;
  fParticleGun->SetParticleTime(time);
}
