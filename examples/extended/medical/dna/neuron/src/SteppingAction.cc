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
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// and papers
// M. Batmunkh et al. J Radiat Res Appl Sci 8 (2015) 498-507
// O. Belov et al. Physica Medica 32 (2016) 1510-1520
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// -------------------------------------------------------------------
// November 2016
// -------------------------------------------------------------------
//
/// \file SteppingAction.cc 
/// \brief Implementation of the SteppingAction class

#include "Analysis.hh"
#include "SteppingAction.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
//
//#include "NeuronLoadDataFile.hh"
#include "G4ios.hh"  
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "Run.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(RunAction* run)  
: G4UserSteppingAction(), fRunAction(run)
{
 //fNeuronLoadParamz = new NeuronLoadDataFile(); 
 //fpNeuronHits = new NeuronHitCompartments();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{ 

  if (step->GetPostStepPoint()->GetProcessDefinedStep()->
      GetProcessName()!="Transportation")
  {
   G4double edepStep = step->GetTotalEnergyDeposit() ;
   G4VPhysicalVolume* volumeStep = step->GetPreStepPoint()->
                                   GetPhysicalVolume();
   G4TouchableHandle touchStep = step->GetPreStepPoint()->GetTouchableHandle();
   G4VPhysicalVolume* volumeMedium = G4PhysicalVolumeStore::GetInstance()
                                   ->GetVolume("Medium");
   G4VPhysicalVolume* volumeSlice = G4PhysicalVolumeStore::GetInstance()
                                   ->GetVolume("BoundingSlice");
 
   // count processes
   // 
   const G4StepPoint* endPoint = step->GetPostStepPoint();
   const G4VProcess* process   = endPoint->GetProcessDefinedStep();
   Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
   run->CountProcesses(process);

 // Edep in all volume
 if (edepStep > 0. ) 
 {   
   fRunAction->AddEdepALL(edepStep);
   //fpNeuronHits->RegisterEdep(volumeStep, edepStep);

   // Edep in a target volume (Bounding slice)
   if (volumeStep == volumeSlice) 
   {
     fRunAction->AddEdepSlice(edepStep);

   }

   // Edep in Soma
   if (volumeStep->GetName() == "Soma")
   {
     fRunAction->AddEdepSoma(edepStep);  
  // fRunAction->AddSomaCompart(touchStep->GetCopyNumber(), edepStep);
     run->AddSomaCompart(touchStep->GetCopyNumber(), edepStep);    
   }

   // Edep in Dendrites
   if (volumeStep->GetName() == "Dendrites")
   {
     fRunAction->AddEdepDend(edepStep);
     run->AddDendCompart(touchStep->GetCopyNumber(), edepStep);
   }

   // Edep in Axon
   if (volumeStep->GetName() == "Axon")
   {   
     fRunAction->AddEdepAxon(edepStep); 
     run->AddAxonCompart(touchStep->GetCopyNumber(), edepStep);
   } 
  
   // Edep in whole Neuron
   if (volumeStep != volumeMedium && volumeStep != volumeSlice)
   {
     fRunAction->AddEdepNeuron(edepStep); 
   }

   // Edep outside bounding scice
   if (volumeStep == volumeMedium)
   {
     fRunAction->AddEdepMedium(edepStep);
   } 
 
 } // end edep>0  
  }

}    
