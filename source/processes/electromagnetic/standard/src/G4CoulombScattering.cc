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
// $Id: G4CoulombScattering.cc,v 1.28 2010-05-25 18:41:12 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4CoulombScattering
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 22.08.2004
//
// Modifications:
// 01.08.06 V.Ivanchenko add choice between G4eCoulombScatteringModel and
//          G4CoulombScatteringModel
//

//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4CoulombScattering.hh"
#include "G4CoulombScatteringModel.hh"
#include "G4eCoulombScatteringModel.hh"
//#include "G4hCoulombScatteringModel.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4LossTableManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4CoulombScattering::G4CoulombScattering(const G4String& name)
  : G4VEmProcess(name),thetaMin(0.0),thetaMax(pi),q2Max(TeV*TeV),
    isInitialised(false)
{
  //  G4cout << "G4CoulombScattering constructor "<< G4endl;
  SetBuildTableFlag(true);
  SetStartFromNullFlag(false);
  SetIntegral(true);
  thEnergy = PeV;
  thEnergyElec = PeV;
  if(name == "CoulombScat") {
    thEnergy = 10.*MeV;
    thEnergyElec = 10.*GeV;
  }
  SetSecondaryParticle(G4Electron::Electron());
  SetProcessSubType(fCoulombScattering);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CoulombScattering::~G4CoulombScattering()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4CoulombScattering::IsApplicable(const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && !p.IsShortLived());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CoulombScattering::InitialiseProcess(const G4ParticleDefinition* p)
{
  //G4cout << "### G4CoulombScattering::InitialiseProcess : "
  //	 << p->GetParticleName() << G4endl;
  G4double a = 
    G4LossTableManager::Instance()->FactorForAngleLimit()*CLHEP::hbarc/CLHEP::fermi;
  q2Max = 0.5*a*a;
 
  // second initialisation
  if(isInitialised) {
    G4VEmModel* mod = GetModelByIndex(0);
    mod->SetPolarAngleLimit(PolarAngleLimit());
    mod = GetModelByIndex(1);
    if(mod) { mod->SetPolarAngleLimit(PolarAngleLimit()); }

    // first initialisation
  } else {
    isInitialised = true;
    aParticle = p;
    G4double mass = p->GetPDGMass();
    G4String name = p->GetParticleName();
    //G4cout << name << "  type: " << p->GetParticleType() 
    //<< " mass= " << mass << G4endl;
    if (mass > GeV || p->GetParticleType() == "nucleus") {
      SetBuildTableFlag(false);
      if(name != "GenericIon") { SetVerboseLevel(0); }
    } else {
      if(name != "e-" && name != "e+" &&
         name != "mu+" && name != "mu-" && name != "pi+" && 
	 name != "kaon+" && name != "proton" ) { SetVerboseLevel(0); }
    }

    G4double emin = MinKinEnergy();
    G4double emax = MaxKinEnergy();
    //    G4gCoulombScatteringModel* model = new G4hCoulombScatteringModel();
    G4eCoulombScatteringModel* model = new G4eCoulombScatteringModel();
    model->SetPolarAngleLimit(PolarAngleLimit());
    model->SetLowEnergyLimit(emin);
    model->SetHighEnergyLimit(emax);
    AddEmModel(1, model);
    /*

    G4double eth  = thEnergy;
    if(mass < MeV) eth  = thEnergyElec; 
    if(eth > emin) {
      G4eCoulombScatteringModel* model = new G4eCoulombScatteringModel();
      model->SetPolarAngleLimit(PolarAngleLimit());
      model->SetLowEnergyLimit(emin);
      model->SetHighEnergyLimit(std::min(eth,emax));
      AddEmModel(1, model);
    }
    if(eth < emax) {
      G4CoulombScatteringModel* model = new G4CoulombScatteringModel();
      model->SetPolarAngleLimit(PolarAngleLimit());
      model->SetLowEnergyLimit(eth);
      model->SetHighEnergyLimit(emax);
      AddEmModel(2, model);
    }
    */
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CoulombScattering::PrintInfo()
{
  G4cout << "      " << PolarAngleLimit()/degree
	 << " < Theta(degree) < 180"
	 << ", Eth(MeV)= ";

  if(aParticle->GetPDGMass() < MeV) G4cout << thEnergyElec;
  else                              G4cout << thEnergy;

  if(q2Max < DBL_MAX) { G4cout << "; pLimit(GeV^1)= " << sqrt(q2Max)/GeV; }
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
