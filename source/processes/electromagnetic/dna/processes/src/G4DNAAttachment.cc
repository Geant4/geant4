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
// $Id: G4DNAAttachment.cc 91992 2015-08-13 07:20:24Z gcosmo $
// GEANT4 tag $Name:  $

#include "G4DNAAttachment.hh"
#include "G4LEPTSAttachmentModel.hh"

#include "G4SystemOfUnits.hh"
#include "G4Positron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4DNAAttachment::G4DNAAttachment(const G4String& processName,
  G4ProcessType type):G4VEmProcess (processName, type),
    isInitialised(false)
{
  SetProcessSubType(55);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4DNAAttachment::~G4DNAAttachment()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4DNAAttachment::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Electron::Electron() || &p == G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAAttachment::InitialiseProcess(const G4ParticleDefinition* p)
{
  if(!isInitialised) 
  {
    isInitialised = true;
    SetBuildTableFlag(false);

    G4String name = p->GetParticleName();

    if(name == "e-" )
    {
      if(!EmModel())
      {
        SetEmModel(new G4DNAMeltonAttachmentModel);
        EmModel()->SetLowEnergyLimit(4.*eV);
        EmModel()->SetHighEnergyLimit(13.*eV);
      }
      AddEmModel(1, EmModel());   
    }
    else if(name == "e+" )
    {
      if(!EmModel())
      {
        SetEmModel(new G4LEPTSAttachmentModel);
        EmModel()->SetLowEnergyLimit(1.*eV);
        EmModel()->SetHighEnergyLimit(1.*MeV);
      }
      AddEmModel(1, EmModel());   
    }
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DNAAttachment::PrintInfo()
{
  G4cout << " Total cross sections computed from "
         << EmModel()->GetName()
         << G4endl;
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
