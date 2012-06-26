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
// G4MuElecInelastic.cc, 2011/08/29 A.Valentin, M. Raine
//
// Based on the following publications
//
//          - Inelastic cross-sections of low energy electrons in silicon
//	    for the simulation of heavy ion tracks with theGeant4-DNA toolkit,
//	    NSS Conf. Record 2010, p80-85
//	    - Geant4 physics processes for microdosimetry simulation:
//	    very low energy electromagnetic models for electrons in Si,
//	    to be published in TNS
//	    - Geant4 physics processes for microdosimetry simulation:
//	    very low energy electromagnetic models for protons and
//	    heavy ions in Si, to be published in NIMB
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 


#include "G4MuElecInelastic.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4MuElecInelastic::G4MuElecInelastic(const G4String& processName,
  G4ProcessType type):G4VEmProcess (processName, type),
    isInitialised(false)
{
  SetProcessSubType(53);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4MuElecInelastic::~G4MuElecInelastic()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4MuElecInelastic::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Electron::Electron() ||
	&p == G4Proton::Proton()  || 
	(p.GetPDGCharge() != 0.0 && !p.IsShortLived() && p.GetParticleType() == "nucleus"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuElecInelastic::InitialiseProcess(const G4ParticleDefinition* p)
{
  if(!isInitialised) 
  {
    isInitialised = true;
    SetBuildTableFlag(false);
    G4String name = p->GetParticleName();

    if(name == "e-")
    {
      if(!Model()) SetModel(new G4MuElecInelasticModel);
      Model()->SetLowEnergyLimit(16.7*eV);
      Model()->SetHighEnergyLimit(50*keV);

      AddEmModel(1, Model());   
    }

    else if(name == "proton")
    {
      if(!Model()) SetModel(new G4MuElecInelasticModel);
      Model()->SetLowEnergyLimit(50.*keV);
      Model()->SetHighEnergyLimit(50*MeV);

      AddEmModel(1, Model());   
    }

    else
    {
      if(!Model()) SetModel(new G4MuElecInelasticModel);
      Model()->SetLowEnergyLimit(50.*keV);
      Model()->SetHighEnergyLimit(100.*GeV);

      AddEmModel(1, Model());   
    }
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuElecInelastic::PrintInfo()
{
  if (Model(2))
  {
    G4cout
      << " Total cross sections computed from " 
      << Model(1)->GetName() 
      << " and "
      << Model(2)->GetName() 
      << " models"
      << G4endl;
  } 
  else
  {
    G4cout
      << " Total cross sections computed from " 
      << Model()->GetName() 
      << G4endl;
  }
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
