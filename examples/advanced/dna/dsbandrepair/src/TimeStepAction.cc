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
/// \file TimeStepAction.cc
/// \brief Implementation of the TimeStepAction class

#include "TimeStepAction.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4DNAChemistryManager.hh"
#include "G4MoleculeCounter.hh"
#include "G4MoleculeTable.hh"
#include "G4Scheduler.hh"
#include "UserMolecule.hh"
#include "G4ChemTimeStepModel.hh"
#include "G4EmParameters.hh"
#include "Analysis.hh"
#include "G4DNAMolecule.hh"
#include "G4H3O.hh"
#include "G4OH.hh"
#include "G4H2O2.hh"
#include "G4H2.hh"
#include "G4Hydrogen.hh"
#include "G4Electron_aq.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction::TimeStepAction()
: G4UserTimeStepAction()
{
    auto param = G4EmParameters::Instance();
    if (param->GetTimeStepModel() == G4ChemTimeStepModel::SBS) {
        AddTimeStep(1*picosecond,0.1*picosecond);
        AddTimeStep(10*picosecond,1*picosecond);
        AddTimeStep(100*picosecond,3*picosecond);
        AddTimeStep(1000*picosecond,10*picosecond);
        AddTimeStep(10000*picosecond,100*picosecond);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


TimeStepAction::TimeStepAction(const TimeStepAction& other) 
: G4UserTimeStepAction(other)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction& TimeStepAction::operator=(const TimeStepAction& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::StartProcessing()
{

    G4cout<<"Start chemistry "<<G4Scheduler::Instance()->GetGlobalTime()/s<<" s"<<G4endl;
    G4cout<<"Chemical endTime "<<G4Scheduler::Instance()->GetEndTime()/nanosecond<<" ns"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::UserReactionAction(const G4Track& a, const G4Track& b,
    const std::vector<G4Track*>* products)
{

    // ***************************************
    // Flag Part
    // ***************************************
    //
    // set the flags

    fReactif1 = 0;
    fReactif2 = 0;
    fProduct1 = 0;
    fProduct2 = 0;

    fReactif1 = SetParticleFlag(a.GetDynamicParticle()->GetDefinition());
    fReactif2 = SetParticleFlag(b.GetDynamicParticle()->GetDefinition());
    if(products)
    {
        for(unsigned int i=0;i<products->size();++i)
        {
            switch(i)
            {
                case 0:
                fProduct1 = SetParticleFlag(((*products)[i])->GetDynamicParticle()->GetDefinition());
                break;
                case 1:
                fProduct2 = SetParticleFlag(((*products)[i])->GetDynamicParticle()->GetDefinition());;
                break;
                default:
                break;
            }
        }
    }

    // ***************************************
    // Save Part
    // ***************************************
    // For DBScan
    G4int base;

    if( (fReactif1==3 || fReactif2 ==3) // OH
        && (
            (fReactif1==8 || fReactif2 ==8) // Desoxyribose
                || (fReactif1==9 || fReactif2 ==9) // Phosphate
		|| (fReactif1==10|| fReactif2 ==10) // Adenine
	    	|| (fReactif1==11|| fReactif2 ==11) // Thymine
	    	|| (fReactif1==12|| fReactif2 ==12) // Guanine
	    	|| (fReactif1==13|| fReactif2 ==13) // Cytosine
                ) )
    {
        // Retrieve the dna molecule
        const G4Track* dnaMolecule = nullptr;
        const G4Track* radical = nullptr;
        //
        if (a.GetDynamicParticle()->GetDefinition() == G4OH::Definition()) {
            dnaMolecule = &b;
            radical = &a;
        }
        else {
            dnaMolecule = &a;
            radical = &b;
        }

        auto _dnaMolecule = dynamic_cast<UserMolecule*> (dnaMolecule->GetUserInformation());
        if (_dnaMolecule) {
            G4int copyNo = _dnaMolecule->GetCopyNumber();
            G4int strand = _dnaMolecule->GetStrand();
            
            if((fReactif1==8 || fReactif2 ==8) // Desoxyribose
                || (fReactif1==9 || fReactif2 ==9)) // Phosphate
                {
                    base=0;
                }
            else 
                {
                    base=1;
                }

            // Retrieve the molecule coordinates
            auto analysisManager = Analysis::GetAnalysis()->GetAnalysisManager();
            analysisManager->FillNtupleIColumn(1, 0, strand);
            analysisManager->FillNtupleIColumn(1, 1, copyNo);
            analysisManager->FillNtupleDColumn(1, 2, radical->GetPosition().getX()/nm);
            analysisManager->FillNtupleDColumn(1, 3, radical->GetPosition().getY()/nm);
            analysisManager->FillNtupleDColumn(1, 4, radical->GetPosition().getZ()/nm);
            analysisManager->FillNtupleDColumn(1, 5, G4Scheduler::Instance()->GetGlobalTime()/ns );
            analysisManager->FillNtupleIColumn(1, 6, base);
            analysisManager->AddNtupleRow(1);
        } else {
            G4String msg = "Error in Downcasting";
            G4Exception("TimeStepAction::UserReactionAction", "", FatalException, msg);
        }
        
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int TimeStepAction::SetParticleFlag(const G4ParticleDefinition* partDef)
{
    G4int flag (0);

    if(partDef == G4H3O::Definition()) flag = 1;
    else if(partDef == G4OH::Definition()) flag = 3;
    else if(partDef == G4Electron_aq::Definition()) flag = 4;
    else if(partDef == G4Hydrogen::Definition()) flag = 5;
    else if(partDef == G4H2::Definition()) flag = 6;
    else if(partDef == G4H2O2::Definition()) flag = 7;
    else if(partDef == G4Deoxyribose::Definition()) flag = 8;
    else if(partDef == G4Phosphate::Definition()) flag = 9;
    else if(partDef == G4Adenine::Definition()) flag = 10;
    else if(partDef == G4Thymine::Definition()) flag = 11;
    else if(partDef == G4Guanine::Definition()) flag = 12;
    else if(partDef == G4Cytosine::Definition()) flag = 13;
    else if(partDef == G4Histone::Definition()) flag = 14;
    else if(partDef == G4DamagedDeoxyribose::Definition()) flag = 15;
    else if(partDef == G4DamagedAdenine::Definition()) flag = 16;
    else if(partDef == G4DamagedThymine::Definition()) flag = 17;
    else if(partDef == G4DamagedCytosine::Definition()) flag = 18;
    else if(partDef == G4DamagedGuanine::Definition()) flag = 19;
    else {
        G4ParticleDefinition* OHm =G4ParticleTable::GetParticleTable()->FindParticle("OHm");
        if (OHm && partDef == OHm) {
            flag = 2;
        } 
    }
    return flag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
