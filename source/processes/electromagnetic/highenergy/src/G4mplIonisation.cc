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
// $Id: G4mplIonisation.cc 106715 2017-10-20 09:39:06Z gcosmo $
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
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4mplIonisationModel.hh"
#include "G4mplIonisationWithDeltaModel.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4mplIonisation::G4mplIonisation(G4double mCharge, const G4String& name)
  : G4VEnergyLossProcess(name),
    magneticCharge(mCharge),
    isInitialised(false)
{
  // By default classical magnetic charge is used
  if(magneticCharge == 0.0) { magneticCharge = eplus*0.5/fine_structure_const; }

  SetVerboseLevel(0);
  SetProcessSubType(fIonisation);
  SetStepFunction(0.2, 1*mm);
  SetSecondaryParticle(G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4mplIonisation::~G4mplIonisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4mplIonisation::IsApplicable(const G4ParticleDefinition&)
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4mplIonisation::InitialiseEnergyLossProcess(const G4ParticleDefinition* p,
						  const G4ParticleDefinition*)
{
  if(isInitialised) { return; }

  SetBaseParticle(0);

  // monopole model is responsible both for energy loss and fluctuations
  G4mplIonisationWithDeltaModel* ion =
    new G4mplIonisationWithDeltaModel(magneticCharge,"PAI");
  ion->SetParticle(p);

  // define size of dedx and range tables
  G4EmParameters* param = G4EmParameters::Instance();
  G4double emin  = std::min(param->MinKinEnergy(),ion->LowEnergyLimit());
  G4double emax  = std::max(param->MaxKinEnergy(),ion->HighEnergyLimit());
  G4int bin = G4lrint(param->NumberOfBinsPerDecade()*std::log10(emax/emin));
  ion->SetLowEnergyLimit(emin);
  ion->SetHighEnergyLimit(emax);
  SetMinKinEnergy(emin);
  SetMaxKinEnergy(emax);
  SetDEDXBinning(bin);

  SetEmModel(ion);
  AddEmModel(1,ion,ion);

  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4mplIonisation::PrintInfo()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4mplIonisation::ProcessDescription(std::ostream& out) const
{
  out << "No description available.";
  out << "<br>\n";
  G4VEnergyLossProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

