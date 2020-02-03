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
// Phys. Med. Biol. 63(10) (2018) 105014-12pp
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class

#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4EmDNAChemistry.hh"
#include "G4EmDNAChemistry_option1.hh"
#include "G4EmDNAPhysics.hh"
#include "G4EmDNAPhysics_option1.hh"
#include "G4EmDNAPhysics_option2.hh"
#include "G4EmDNAPhysics_option3.hh"
#include "G4EmDNAPhysics_option4.hh"
#include "G4EmDNAPhysics_option5.hh"
#include "G4EmDNAPhysics_option6.hh"
#include "G4EmDNAPhysics_option7.hh"
#include "G4EmDNAPhysics_option8.hh"
#include "G4PhysicsConstructorRegistry.hh"
#include "CommandLineParser.hh"
#include "G4EmParameters.hh"

using namespace G4DNAPARSER;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
  : G4VModularPhysicsList(),
    fEmDNAPhysicsList(nullptr),fEmDNAChemistryList(nullptr),
    fEmDNAChemistryList1(nullptr),fPhysDNAName("")
{
  G4double currentDefaultCut = 1.*nanometer;
  // fixe lower limit for cut
  G4ProductionCutsTable::GetProductionCutsTable()->
    SetEnergyRange(100*eV, 1*GeV);
  SetDefaultCutValue(currentDefaultCut);
  SetVerboseLevel(1);
  
  RegisterConstructor("G4EmDNAPhysics");
  if(CommandLineParser::GetParser()->GetCommandIfActive("-chemOFF")==0)
    {
      RegisterConstructor("G4EmDNAChemistry");
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete fEmDNAPhysicsList;
  delete fEmDNAChemistryList;
  delete fEmDNAChemistryList1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  if(fEmDNAPhysicsList)    { fEmDNAPhysicsList->ConstructParticle(); }
  if(fEmDNAChemistryList)  { fEmDNAChemistryList->ConstructParticle(); }
  if(fEmDNAChemistryList1) { fEmDNAChemistryList1->ConstructParticle(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  if(fEmDNAPhysicsList)    { fEmDNAPhysicsList->ConstructProcess(); }
  if(fEmDNAChemistryList)  { fEmDNAChemistryList->ConstructProcess(); }
  if(fEmDNAChemistryList1) { fEmDNAChemistryList1->ConstructProcess(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::RegisterConstructor(const G4String& name)
{
  if(name == fPhysDNAName) { return; }
  if(verboseLevel > 0) {
    G4cout << "===== Register constructor ==== " << name << G4endl; 
  }
  if(name == "G4EmDNAPhysics") {
    delete fEmDNAPhysicsList;
    fEmDNAPhysicsList = new G4EmDNAPhysics(verboseLevel);
    fPhysDNAName = name;
  } else if(name == "G4EmDNAPhysics_option1") {
    delete fEmDNAPhysicsList;
    fEmDNAPhysicsList = new G4EmDNAPhysics_option1(verboseLevel);
    fPhysDNAName = name;
  } else if(name == "G4EmDNAPhysics_option2") {
    delete fEmDNAPhysicsList;
    fEmDNAPhysicsList = new G4EmDNAPhysics_option2(verboseLevel);
    fPhysDNAName = name;
  } else if(name == "G4EmDNAPhysics_option3") {
    delete fEmDNAPhysicsList;
    fEmDNAPhysicsList = new G4EmDNAPhysics_option3(verboseLevel);
    fPhysDNAName = name;
  } else if(name == "G4EmDNAPhysics_option4") {
    delete fEmDNAPhysicsList;
    fEmDNAPhysicsList = new G4EmDNAPhysics_option4(verboseLevel);
    fPhysDNAName = name;
  } else if(name == "G4EmDNAPhysics_option5") {
    delete fEmDNAPhysicsList;
    fEmDNAPhysicsList = new G4EmDNAPhysics_option5(verboseLevel);
    fPhysDNAName = name;
  } else if(name == "G4EmDNAPhysics_option6") {
    delete fEmDNAPhysicsList;
    fEmDNAPhysicsList = new G4EmDNAPhysics_option6(verboseLevel);
    fPhysDNAName = name;
  } else if(name == "G4EmDNAPhysics_option7") {
    delete fEmDNAPhysicsList;
    fEmDNAPhysicsList = new G4EmDNAPhysics_option7(verboseLevel);
    fPhysDNAName = name;
  } else if(name == "G4EmDNAPhysics_option8") {
    delete fEmDNAPhysicsList;
    fEmDNAPhysicsList = new G4EmDNAPhysics_option8(verboseLevel);
    fPhysDNAName = name;
  } else if(name == "G4EmDNAChemistry") {
    if(fEmDNAChemistryList || fEmDNAChemistryList1) { return; }
    fEmDNAChemistryList = new G4EmDNAChemistry();
    fEmDNAChemistryList->SetVerboseLevel(verboseLevel);
  } else if(name == "G4EmDNAChemistry_option1") {
    if(fEmDNAChemistryList || fEmDNAChemistryList1) { return; }
    fEmDNAChemistryList1 = new G4EmDNAChemistry_option1();
    fEmDNAChemistryList1->SetVerboseLevel(verboseLevel);
  } else {
    G4cout << "PhysicsList::RegisterConstructor: <" << name << ">"
           << " fails - name is not defined"
           << G4endl;    
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
