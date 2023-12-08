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
/// \file ChemPhysicsList.cc
/// \brief Implementation of the ChemPhysicsList class

#include "ChemPhysicsList.hh"
#include "G4EmDNAChemistry_option3_Extended.hh"

#include "G4SystemOfUnits.hh"
#include "G4EmDNAPhysics.hh"
#include "G4EmDNAPhysics_option1.hh"
#include "G4EmDNAPhysics_option2.hh"
#include "G4EmDNAPhysics_option3.hh"
#include "G4EmDNAPhysics_option4.hh"
#include "G4EmDNAPhysics_option5.hh"
#include "G4EmDNAPhysics_option6.hh"
#include "G4EmDNAPhysics_option7.hh"
#include "G4EmDNAPhysics_option8.hh"
#include "G4EmDNAChemistry.hh"
#include "G4EmDNAChemistry_option1.hh"
#include "G4EmDNAChemistry_option2.hh"
#include "G4EmDNAChemistry_option3.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4PhysicsConstructorRegistry.hh"
#include "G4LeptonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4EmParameters.hh"
#include "G4ChemicalMoleculeFinder.hh"
#include "G4Filesystem.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChemPhysicsList::ChemPhysicsList(G4String chemListName)
    : G4VModularPhysicsList()
{
    fmsg = new ChemPhysicsMessenger(this);
    SetDefaultCutValue(1.0*nanometer);
    SetVerboseLevel(1);
    RegisterPhysListConstructor(ExtractPhysDNAName());
    RegisterChemListConstructor(chemListName);
    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV, 1*GeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemPhysicsList::ConstructParticle()
{
    if(fDNAPhysicsList != nullptr)    
    {
        fDNAPhysicsList->ConstructParticle(); 
    }
    
    if(fChemistryList_option2 != nullptr)
    {
        fChemistryList_option2->ConstructParticle();
    }
    if(fChemistryList_option3_mod != nullptr)
    {
        fChemistryList_option3_mod->ConstructParticle();
    }
    
    // construct following particles to get rid of warning
    G4LeptonConstructor lConstructor;
    lConstructor.ConstructParticle();
    G4BosonConstructor  pBosonConstructor;
    pBosonConstructor.ConstructParticle();
    G4MesonConstructor pMesonConstructor;
    pMesonConstructor.ConstructParticle();
    G4BaryonConstructor pBaryonConstructor;
    pBaryonConstructor.ConstructParticle();
    G4ShortLivedConstructor pShortLivedConstructor;
    pShortLivedConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemPhysicsList::ConstructProcess()
{  
    AddTransportation();

    if(fDNAPhysicsList != nullptr)    
    {
        fDNAPhysicsList->ConstructProcess(); 
    }
    
    if(fChemistryList_option2 != nullptr)
    {
        fChemistryList_option2->ConstructProcess();
    }
    
    if(fChemistryList_option3_mod != nullptr)
    {
        fChemistryList_option3_mod->ConstructProcess();
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemPhysicsList::RegisterPhysListConstructor(const G4String& name)
{
    if(name == fPhysDNAName ) 
    {
        return; 
    }
    
    if(verboseLevel > 0) 
    {
        if (fPhysDNAName == "") G4cout << "===== Register constructor ==== " << name << G4endl; 
        else G4cout << "===== Change Physics constructor from "<<fPhysDNAName<<" to "<<name<<" ==== " << G4endl; 
    }
    
    if(name == "G4EmDNAPhysics") 
    {
        fDNAPhysicsList.reset(new G4EmDNAPhysics(verboseLevel));
        fPhysDNAName = name;
    } 
    else if(name == "G4EmDNAPhysics_option1") 
    {
        fDNAPhysicsList.reset(new G4EmDNAPhysics_option1(verboseLevel));
        fPhysDNAName = name;
    } 
    else if(name == "G4EmDNAPhysics_option2") 
    {
        fDNAPhysicsList.reset(new G4EmDNAPhysics_option2(verboseLevel));
        fPhysDNAName = name;
    } 
    else if(name == "G4EmDNAPhysics_option3") 
    {
        fDNAPhysicsList.reset(new G4EmDNAPhysics_option3(verboseLevel));
        fPhysDNAName = name;
    } 
    else if(name == "G4EmDNAPhysics_option4") 
    {
        fDNAPhysicsList.reset(new G4EmDNAPhysics_option4(verboseLevel));
        fPhysDNAName = name;
    } 
    else if(name == "G4EmDNAPhysics_option5") 
    {
        fDNAPhysicsList.reset(new G4EmDNAPhysics_option5(verboseLevel));
        fPhysDNAName = name;
    } 
    else if(name == "G4EmDNAPhysics_option6") 
    {
        fDNAPhysicsList.reset(new G4EmDNAPhysics_option6(verboseLevel));
        fPhysDNAName = name;
    } 
    else if(name == "G4EmDNAPhysics_option7") 
    {
        fDNAPhysicsList.reset(new G4EmDNAPhysics_option7(verboseLevel));
        fPhysDNAName = name;
    } 
    else if(name == "G4EmDNAPhysics_option8") 
    {
        fDNAPhysicsList.reset(new G4EmDNAPhysics_option8(verboseLevel));
        fPhysDNAName = name;
    } 
    else 
    {
        G4cout << "PhysicsList::RegisterPhysListConstructor: <" << name << ">"
               << " fails - name is not defined"
               << G4endl;    
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemPhysicsList::RegisterChemListConstructor(const G4String& name)
{
    if(name == fChemListName ) 
    {
        return; 
    }
    
    if(verboseLevel > 0) 
    {
        G4cout << "===== Register constructor for Chemitry ==== " << name << G4endl; 
    }
  
    if(name == "G4EmDNAChemistry_option2")
    {
        if ( fChemistryList_option2 != nullptr) return;
        fChemistryList_option2.reset(new G4EmDNAChemistry_option2());
        fChemistryList_option2->SetVerboseLevel(verboseLevel);
        fTimeStepModel = fSBS;
        fChemListName = name;
    }
    else if(name == "G4EmDNAChemistry_option3") 
    {
        if ( fChemistryList_option3_mod != nullptr) return;
        fChemistryList_option3_mod.reset(new G4EmDNAChemistry_option3_Extended());
        fChemistryList_option3_mod->SetVerboseLevel(verboseLevel);
        fTimeStepModel = fIRT_syn;
        fChemistryList_option3_mod->SetTimeStepModel(fTimeStepModel);
        fChemListName = name;
    }
    else 
    {
        G4ExceptionDescription msg;
        msg <<"ChemPhysicsList::RegisterChemListConstructor: <" << name << ">"
            <<" fails - name is not defined";
        G4Exception("ChemPhysicsList::RegisterChemListConstructor", "Phys_WrongName", 
        JustWarning, msg); 
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String ChemPhysicsList::ExtractPhysDNAName()
{
    G4String out ="G4EmDNAPhysics_option2";
    G4fs::path thisP = G4fs::current_path();
    for (const auto &entry : G4fs::directory_iterator(thisP)){
        if (entry.path().filename() == "imp.info") {
            std::ifstream file(entry.path().c_str());
            if(!file.good() ){
                G4ExceptionDescription msg;
                msg<<"**** Fatal Error *****\n";
                msg<<"ChemPhysicsList::ExtractPhysDNAName(): File corupted: "
                <<entry.path()<<"\n";
                msg<<"*************** *****";
                G4Exception("ChemPhysicsList::ExtractPhysDNAName", "Phys_WrongParse", 
                FatalException, msg);
            }

            G4String line;
            while(std::getline(file, line) ){
                std::istringstream iss(line);
                G4String flag;
                iss >> flag;
                if ( flag == "_physList") {
                    iss >> out;
                }
            }
            file.close();
        }
    }
    return out;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......