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
// $Id: G4DNAExcitation.cc,v 1.7 2010-10-08 08:53:17 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4DNAExcitation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4DNAExcitation::G4DNAExcitation(const G4String& processName,
  G4ProcessType type):G4VEmProcess (processName, type),
    isInitialised(false)
{
  SetProcessSubType(51);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4DNAExcitation::~G4DNAExcitation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4DNAExcitation::IsApplicable(const G4ParticleDefinition& p)
{

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  return 
    (
       &p == G4Electron::Electron() 
    || &p == G4Proton::ProtonDefinition()
    || &p == instance->GetIon("hydrogen")
    || &p == instance->GetIon("alpha++")
    || &p == instance->GetIon("alpha+")
    || &p == instance->GetIon("helium")
    );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAExcitation::InitialiseProcess(const G4ParticleDefinition* p)
{
  if(!isInitialised) 
  {
    isInitialised = true;
    SetBuildTableFlag(false);
    
    G4String name = p->GetParticleName();

    if(name == "e-")
    {

      // Emfietzoglou model
/*
      if(!Model()) SetModel(new G4DNAEmfietzoglouExcitationModel);
      Model()->SetLowEnergyLimit(8.23*eV);
      Model()->SetHighEnergyLimit(10*MeV);
*/
      // Born model

      if(!Model()) SetModel(new G4DNABornExcitationModel);
      Model()->SetLowEnergyLimit(9*eV);
      Model()->SetHighEnergyLimit(1*MeV);

      AddEmModel(1, Model());   
    }
    
    if(name == "proton")
    {
      if(!Model(1)) SetModel(new G4DNAMillerGreenExcitationModel,1);
      Model(1)->SetLowEnergyLimit(10*eV);
      Model(1)->SetHighEnergyLimit(500*keV);

      if(!Model(2)) SetModel(new G4DNABornExcitationModel,2);
      Model(2)->SetLowEnergyLimit(500*keV);
      Model(2)->SetHighEnergyLimit(100*MeV);
    
      AddEmModel(1, Model(1));   
      AddEmModel(2, Model(2));   
    }

    if(name == "hydrogen")
    {
      if(!Model()) SetModel(new G4DNAMillerGreenExcitationModel);
      Model()->SetLowEnergyLimit(10*eV);
      Model()->SetHighEnergyLimit(500*keV);
   
      AddEmModel(1, Model());   
    }


    if( name == "alpha" || name == "alpha+" || name == "helium" )
    {
      if(!Model()) SetModel(new G4DNAMillerGreenExcitationModel);
      Model()->SetLowEnergyLimit(1*keV);
      Model()->SetHighEnergyLimit(400*MeV);

      AddEmModel(1, Model());   
    }

  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DNAExcitation::PrintInfo()
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
