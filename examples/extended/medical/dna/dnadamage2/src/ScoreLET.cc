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
// This example is provided by the Geant4-DNA collaboration
// chem6 example is derived from chem4 and chem5 examples
//
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// J. Appl. Phys. 125 (2019) 104301
// Med. Phys. 45 (2018) e722-e739
// J. Comput. Phys. 274 (2014) 841-882
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157-178
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: W. G. Shin and S. Incerti (CENBG, France)
//
// $Id$
//
/// \file ScoreLET.cc
/// \brief Implementation of the ScoreLET class

#include "ScoreLET.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ScoreLET::ScoreLET(G4String name)
  :G4VPrimitiveScorer(name),G4UImessenger(),fEvtMap(0)
{
  fpLETDir = new G4UIdirectory("/scorer/LET/");
  fpLETDir->SetGuidance("LET scorer commands");

  fpCutoff = new G4UIcmdWithADoubleAndUnit("/scorer/LET/cutoff", this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ScoreLET::~ScoreLET()
{
  delete fpLETDir;
  delete fpCutoff;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreLET::SetNewValue(G4UIcommand* command, G4String newValue){
  if(command == fpCutoff) fCutoff = atof(newValue);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool ScoreLET::ProcessHits(G4Step* aStep,G4TouchableHistory* /*TH*/)
{
  // In order to follow the primary track
  // regardless charge increasing or decreasing
  if(aStep->GetTrack()->GetTrackID() != 1 &&
     aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding() != 11){
    G4int subType = aStep->GetTrack()->GetCreatorProcess()
                    ->GetProcessSubType();
    if(subType == 56 || subType == 57){
       fTrackID = aStep->GetTrack()->GetTrackID();
    }
  }

  // Ignore the step if it is not primary.
  if(aStep->GetTrack()->GetTrackID() != fTrackID) return false;
  else{
    fStepL += aStep->GetStepLength()/um;
    fEdep += aStep->GetTotalEnergyDeposit()/keV;

    G4int subType = aStep->GetPostStepPoint()->
                    GetProcessDefinedStep()->GetProcessSubType();

    // Don't add the kinetic energy of primary particle
    if(subType == 56 || subType == 57) return false;

    const std::vector<const G4Track*>* secondary =
                                       aStep->GetSecondaryInCurrentStep();

    size_t nbtrk = (*secondary).size();

    if(nbtrk){
      for(size_t lp=0;lp<nbtrk;lp++){
        // Store the kinetic energy of secondaries
        // which less than cutoff energy.
        if((*secondary)[lp]->GetKineticEnergy()/eV<fCutoff){
          fEdep += (*secondary)[lp]->GetKineticEnergy()/keV;
        }
      }
    }
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreLET::Initialize(G4HCofThisEvent* HCE)
{
  fEvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(),
                                     GetName());

  fEvtMap->clear();
  fNEvent = 0;
static G4int HCID = -1;
  if(HCID < 0) {HCID = GetCollectionID(0);}
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)fEvtMap);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreLET::EndOfEvent(G4HCofThisEvent*)
{
  if (fStepL > 0) {
    fLET = fEdep / fStepL;
    if(!G4RunManager::GetRunManager()->GetCurrentEvent()->IsAborted())
    {
      fEvtMap->add(fNEvent,fLET);
      fNEvent++;
    }
  }
  fTrackID = 1;
  fLET = 0;
  fEdep = 0;
  fStepL = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4int ScoreLET::GetIndex(G4Step* /*aStep*/)
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
