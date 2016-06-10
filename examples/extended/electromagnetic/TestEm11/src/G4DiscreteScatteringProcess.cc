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
#include "G4DiscreteScatteringProcess.hh"
#include "G4DiscreteScatteringModel.hh"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DiscreteScatteringProcess::G4DiscreteScatteringProcess(G4int iNumAngles)
: G4VEmProcess("DiscreteScattering", fElectromagnetic), 
  fIsInitialised(false), fNumAngles(iNumAngles)
{
  SetBuildTableFlag(true);
  SetStartFromNullFlag(false);
  SetMinKinEnergy(2*CLHEP::keV);
  SetMaxKinEnergy(100*CLHEP::MeV);
  SetIntegral(false);
  SetSplineFlag(false);
  SetProcessSubType(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DiscreteScatteringProcess::~G4DiscreteScatteringProcess(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool G4DiscreteScatteringProcess::IsApplicable(const G4ParticleDefinition& p){
  return (&p == G4Electron::Definition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DiscreteScatteringProcess::PrintInfo(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4DiscreteScatteringProcess::InitialiseProcess(const G4ParticleDefinition*)
{
  if(fIsInitialised) { return; }
  fIsInitialised = true;
  G4VEmModel* model = new G4DiscreteScatteringModel(fNumAngles);
  AddEmModel(1, model);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
