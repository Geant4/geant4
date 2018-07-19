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
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id$
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "Analysis.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Threading.hh"


//#ifdef G4MULTITHREADED
//  #include "G4MTRunManager.hh"
//#else
  #include "G4RunManager.hh"
//#endif



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::SteppingAction(EventAction* event)
: G4UserSteppingAction(),
  fpEventaction(event)
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::~SteppingAction()
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::UserSteppingAction(const G4Step* step)
{
    //#ifdef G4MULTITHREADED
    //    G4MTRunManager *rm = G4MTRunManager::GetRunManager();
    //#else
    G4RunManager *rm = G4RunManager::GetRunManager();
    //#endif

    G4int eventID= rm->GetCurrentEvent()->GetEventID();

    G4double flagProcess (0.);
    G4double x,y,z;
    G4double dE;
    G4String process_name;
    process_name=step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();


    // only ionizations in the target will be recorded
    if(step->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Target")
    {
        if (process_name=="e-_G4DNAElastic")
            flagProcess=1;
        else if(process_name=="e-_G4DNAExcitation")
            flagProcess=2;
        else if(process_name=="e-_G4DNAIonisation")
            flagProcess=3;

        // Energy deposited and position of the interaction are saved
        dE= step->GetTotalEnergyDeposit()/eV;

        if (dE!=0)
        {
            x=step->GetPostStepPoint()->GetPosition().x()/nanometer;
            y=step->GetPostStepPoint()->GetPosition().y()/nanometer;
            z=step->GetPostStepPoint()->GetPosition().z()/nanometer;

            // get analysis manager
            G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

            // fill ntuple
            analysisManager->FillNtupleDColumn(2,0, eventID);
            analysisManager->FillNtupleDColumn(2,1, flagProcess);
            analysisManager->FillNtupleDColumn(2,2, x);
            analysisManager->FillNtupleDColumn(2,3, y);
            analysisManager->FillNtupleDColumn(2,4, z);
            analysisManager->FillNtupleDColumn(2,5, dE);
            analysisManager->AddNtupleRow(2);

            // histogram
            if (flagProcess==3)
                fpEventaction-> AddEventIn(1);
        }
    }
}    
