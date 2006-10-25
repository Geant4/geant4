//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4mplIonisation.cc,v 1.1 2006-10-25 17:37:44 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4mplIonisation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 25.08.2005
//
// Modifications:
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4mplIonisation.hh"
#include "G4Electron.hh"
#include "G4Monopole.hh"
#include "G4mplIonisationModel.hh"
#include "G4BohrFluctuations.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4mplIonisation::G4mplIonisation(const G4String& name)
  : G4VEnergyLossProcess(name),
    isInitialised(false)
{
  SetDEDXBinning(120);
  SetLambdaBinning(120);
  SetMinKinEnergy(0.1*keV);
  SetMaxKinEnergy(100.0*TeV);
  SetVerboseLevel(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4mplIonisation::~G4mplIonisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4mplIonisation::InitialiseEnergyLossProcess(const G4ParticleDefinition*,
						  const G4ParticleDefinition*)
{
  if(isInitialised) return;

  SetBaseParticle(0);
  SetSecondaryParticle(G4Electron::Electron());

  G4BohrFluctuations* fluct  = new G4BohrFluctuations();   
  G4mplIonisationModel* ion  = new G4mplIonisationModel();
  ion->SetLowEnergyLimit(0.1*keV);
  ion->SetHighEnergyLimit(100.*TeV);
  AddEmModel(0,ion,fluct);

  SetStepLimits(0.2, 0.1*mm);

  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4mplIonisation::PrintInfo()
{
  G4cout << "      No delta-electron production, only dE/dx"
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
