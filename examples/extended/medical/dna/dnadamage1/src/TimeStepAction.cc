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
/// \file TimeStepAction.hh
/// \brief Implementation of the TimeStepAction class

#include "TimeStepAction.hh"
#include <G4Scheduler.hh>
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Molecule.hh"
#include "g4root.hh"
#include "G4DNAMolecule.hh"
#include "G4MoleculeTable.hh"
#include "G4OH.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction::TimeStepAction() 
    : G4UserTimeStepAction()
{
    AddTimeStep(1*picosecond, 0.35*picosecond);
    AddTimeStep(10*picosecond, 1*picosecond);
    AddTimeStep(100*picosecond, 3*picosecond);
    AddTimeStep(1000*picosecond, 10*picosecond);
    AddTimeStep(10000*picosecond, 100*picosecond);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction::~TimeStepAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction::TimeStepAction(const TimeStepAction& other) 
    : G4UserTimeStepAction(other)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction&
TimeStepAction::operator=(const TimeStepAction& rhs)
{
    if (this == &rhs) return *this;
    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::UserReactionAction(const G4Track& trackA,
                                        const G4Track& trackB,
                                        const std::vector<G4Track*>* 
                                        pProducts)
{
    if(pProducts && (GetMolecule((*pProducts)[0])->
    GetDefinition() != G4DamagedDeoxyribose::Definition()))
    {
        return;
    }
    
    //check for the case "no product"
    if(!pProducts)
    {
        return;
    }
    
    auto phyEventId = G4EventManager::GetEventManager()->
    GetConstCurrentEvent()->GetEventID();
    
    const G4Track* DNAElement = nullptr;
    const G4Track* radical    = nullptr;
    if(GetMolecule(&trackA)->
    GetDefinition() == G4Deoxyribose::Definition())
    {
        DNAElement = &trackA;
        radical    = &trackB;
    }
    else 
    {
        DNAElement = &trackB;
        radical    = &trackA;
    }
    
    if(GetMolecule(radical)->GetDefinition() != G4OH::Definition())
    {
        return;
    }
    
    G4AnalysisManager* analysisManager = 
    G4AnalysisManager::Instance();
    analysisManager->FillNtupleDColumn(2, 0, 
    DNAElement->GetPosition().getX()/nm);
    analysisManager->FillNtupleDColumn(2, 1, 
    DNAElement->GetPosition().getY()/nm);
    analysisManager->FillNtupleDColumn(2, 2, 
    DNAElement->GetPosition().getZ()/nm);
    analysisManager->FillNtupleSColumn(2, 3, 
    GetMolecule(radical)->GetName());
    analysisManager->FillNtupleIColumn(2, 4, 
    G4int(phyEventId));
    analysisManager->AddNtupleRow(2);
}

