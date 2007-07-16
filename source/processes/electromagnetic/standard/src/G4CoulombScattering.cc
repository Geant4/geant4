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
// $Id: G4CoulombScattering.cc,v 1.8 2007-07-16 08:45:34 vnivanch Exp $
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
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4CoulombScattering::G4CoulombScattering(const G4String& name)
  : G4VEmProcess(name),thetaMin(0.0),thetaMax(pi),q2Max(DBL_MAX),
    isInitialised(false)
{
  G4VEmProcess::SetBuildTableFlag(true);
  SetStartFromNullFlag(false);
  SetIntegral(true);
  SetMinKinEnergy(keV);
  SetMaxKinEnergy(PeV);
  SetLambdaBinning(120);
  SetSecondaryParticle(G4Electron::Electron());
  buildElmTableFlag = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CoulombScattering::~G4CoulombScattering()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CoulombScattering::InitialiseProcess(const G4ParticleDefinition* p)
{
  if(!isInitialised) {
    isInitialised = true;
    aParticle = p;
    if (p->GetParticleType() == "nucleus") {
      buildElmTableFlag = false;
      verboseLevel = 0;
    } else {
      G4String name = p->GetParticleName();
      if(name != "e-" && name != "e+" &&
         name != "mu+" && name != "mu-" && name != "pi+" && 
	 name != "kaon+" && name != "proton" ) verboseLevel = 0;
    }

    G4double emin = MinKinEnergy();
    G4double emax = MaxKinEnergy();
    if(GetProcessName() == "eCoulombScat") {
      G4eCoulombScatteringModel* model = 
	new G4eCoulombScatteringModel(thetaMin,thetaMax,buildElmTableFlag,q2Max);
      model->SetLowEnergyLimit(emin);
      model->SetHighEnergyLimit(emax);
      AddEmModel(1, model);
    } else {
      G4CoulombScatteringModel* model = 
	new G4CoulombScatteringModel(thetaMin,thetaMax,buildElmTableFlag,q2Max);
      model->SetLowEnergyLimit(emin);
      model->SetHighEnergyLimit(emax);
      AddEmModel(1, model);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CoulombScattering::PrintInfo()
{
  G4cout << " Scattering of " << aParticle->GetParticleName()
	 << " with   " << thetaMin/degree
	 << " < Theta(degree) < " << thetaMax/degree;
  if(q2Max < DBL_MAX) G4cout << "; q2Max(GeV^2)= " << q2Max/(GeV*GeV);
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
