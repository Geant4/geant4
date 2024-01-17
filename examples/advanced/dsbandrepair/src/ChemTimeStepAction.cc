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
/// \file ChemTimeStepAction.cc
/// \brief Implementation of the ChemTimeStepAction class

#include "ChemTimeStepAction.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4DNAChemistryManager.hh"
#include "G4MoleculeCounter.hh"
#include "G4MoleculeTable.hh"
#include "G4Scheduler.hh"
#include "UserMolecule.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChemTimeStepAction::ChemTimeStepAction(ChemNtupleManager* histo, TimeStepModel timeStepModel) 
: G4UserTimeStepAction(),fpHisto(histo)
{
    if (timeStepModel == fSBS) {
        AddTimeStep(1*picosecond,0.1*picosecond);
        AddTimeStep(10*picosecond,1*picosecond);
        AddTimeStep(100*picosecond,3*picosecond);
        AddTimeStep(1000*picosecond,10*picosecond);
        AddTimeStep(10000*picosecond,100*picosecond);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


ChemTimeStepAction::ChemTimeStepAction(const ChemTimeStepAction& other) 
: G4UserTimeStepAction(other)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChemTimeStepAction& ChemTimeStepAction::operator=(const ChemTimeStepAction& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemTimeStepAction::StartProcessing()
{

    G4cout<<"Start chemistry "<<G4Scheduler::Instance()->GetGlobalTime()/s<<" s"<<G4endl;
    G4cout<<"Chemical endTime "<<G4Scheduler::Instance()->GetEndTime()/nanosecond<<" ns"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemTimeStepAction::UserReactionAction(const G4Track& a, const G4Track& b,
    const std::vector<G4Track*>* products)
{

    // ***************************************
    // Flag Part
    // ***************************************
    //
    // set the flags

    fReactif1 = 0.;
    fReactif2 = 0.;
    fProduct1 = 0.;
    fProduct2 = 0.;

    G4String aName = a.GetDynamicParticle()->GetDefinition()->GetParticleName();
    fReactif1 = SetFlag(aName);

    G4String bName = b.GetDynamicParticle()->GetDefinition()->GetParticleName();
    fReactif2 = SetFlag(bName);

    if(products)
    {
        for(unsigned int i=0;i<products->size();++i)
        {
            switch(i)
            {
                case 0:
                fProduct1 = SetFlag(((*products)[i])->GetDynamicParticle()->GetDefinition()->GetParticleName());
                break;
                case 1:
                fProduct2 = SetFlag(((*products)[i])->GetDynamicParticle()->GetDefinition()->GetParticleName());;
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
        const G4Track* dnaMolecule;
        const G4Track* radical;
        //
        if (aName == "OH") {
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
            fpHisto->FillNtupleDColumn(1, 0, strand);
            fpHisto->FillNtupleDColumn(1, 1, copyNo);
            fpHisto->FillNtupleDColumn(1, 2, radical->GetPosition().getX()/nm);
            fpHisto->FillNtupleDColumn(1, 3, radical->GetPosition().getY()/nm);
            fpHisto->FillNtupleDColumn(1, 4, radical->GetPosition().getZ()/nm);
            fpHisto->FillNtupleDColumn(1, 5, G4Scheduler::Instance()->GetGlobalTime()/ns );
            fpHisto->FillNtupleDColumn(1, 6, base);
            fpHisto->AddNtupleRow(1);
        } else {
            G4String msg = "Error in Downcasting";
            G4Exception("ChemTimeStepAction::UserReactionAction", "", FatalException, msg);
        }
        
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ChemTimeStepAction::SetFlag(const G4String& val)
{
    G4double flag (0);

    if(val == "H3O") flag = 1;
    else if(val == "OHm") flag = 2;
    else if(val == "OH") flag = 3;
    else if(val == "e_aq") flag = 4;
    else if(val == "H") flag = 5;
    else if(val == "H_2") flag = 6;
    else if(val == "H2O2") flag = 7;
    else if(val == "Deoxyribose") flag = 8;
    else if(val == "Phosphate") flag = 9;
    else if(val == "Adenine") flag = 10;
    else if(val == "Thymine") flag = 11;
    else if(val == "Guanine") flag = 12;
    else if(val == "Cytosine") flag = 13;
    else if(val == "Histone") flag = 14;
    else if(val == "Damaged_Deoxyribose") flag = 15;
    else if(val == "Damaged_Adenine") flag = 16;
    else if(val == "Damaged_Thymine") flag = 17;
    else if(val == "Damaged_Cytosine") flag = 18;
    else if(val == "Damaged_Guanine") flag = 19;

    return flag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
