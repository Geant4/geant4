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
// -------------------------------------------------------------------
// $Id$
// -------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "PhysicsList.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4EmDNAPhysicsChemistry.hh"
#include "G4EmDNAPhysics.hh"
#include "G4ITStepManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularStepByStepModel.hh"
#include "G4DNASmoluchowskiReactionModel.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhysicsList::PhysicsList():  G4VModularPhysicsList()
{
    defaultCutValue = 1*nanometer;
    cutForGamma     = defaultCutValue;
    cutForElectron  = defaultCutValue;
    cutForPositron  = defaultCutValue;
    cutForProton    = defaultCutValue;

    SetVerboseLevel(1);
    chemistryFlag = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhysicsList::~PhysicsList()
{
    delete emDNAPhysicsList;
    G4ITStepManager::DeleteInstance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructParticle()
{
        if(chemistryFlag)
                emDNAPhysicsList = new G4EmDNAPhysicsChemistry();
        else
                emDNAPhysicsList = new G4EmDNAPhysics();

	emDNAPhysicsList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructProcess()
{
    AddTransportation();
    emDNAPhysicsList->ConstructProcess();

    //__________________________________________________________________
    // Diffusion controlled reaction model

    /**
      * The reaction model defines how to compute the reaction range between molecules
      */

    G4VDNAReactionModel* reactionRadiusComputer = new G4DNASmoluchowskiReactionModel();
    G4DNAMolecularReactionTable::GetReactionTable() -> PrintTable(reactionRadiusComputer);

    /**
      * The StepByStep model tells the step manager how to behave before and after each step,
      * how to compute the time steps.
      */
    G4DNAMolecularStepByStepModel* sbs = new G4DNAMolecularStepByStepModel();
    G4ITStepManager::Instance()->GetModelHandler()->RegisterModel(sbs, 0);
    sbs->SetReactionTable(G4DNAMolecularReactionTable::GetReactionTable());
    sbs->SetReactionModel(reactionRadiusComputer);

    //__________________________________________________________________
    map<double,double>* steps = new map<double, double> ;

    /**
      * Give to G4ITStepManager the user defined time steps
      * eg : from 1 picosecond to 10 picosecond, the minimum time
      * step that the TimeStepper can returned is 0.1 picosecond.
      * Those time steps are used for the chemistry of G4DNA
      */

    (*steps)[1*picosecond] = 0.1*picosecond;
    (*steps)[10*picosecond] = 1*picosecond;
    (*steps)[100*picosecond] = 3*picosecond;
    (*steps)[1000*picosecond] = 10*picosecond;
    (*steps)[10000*picosecond] = 100*picosecond;

    G4ITStepManager::Instance()-> SetTimeSteps(steps);

    /**
      * Tells to the chemistry manager whether the chemistry
      * needs to be activated.
      * WARNING : if you don't use the chemistry do not activate it
      * otherwise it might generate memory leaks with tracks created but
      * not destroyed.
      */
    G4DNAChemistryManager::Instance()->SetChemistryActivation(true);

    G4ITStepManager::Instance()->Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::SetCuts()
{
    if (verboseLevel >0)
    {
        G4cout << "PhysicsList::SetCuts:";
        G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
    }

    // set cut values for gamma at first and for e- second and next for e+,
    // because some processes for e+/e- need cut values for gamma
    SetCutValue(cutForGamma, "gamma");
    SetCutValue(cutForElectron, "e-");
    SetCutValue(cutForPositron, "e+");

    // set cut values for proton and anti_proton before all other hadrons
    // because some processes for hadrons need cut values for proton/anti_proton
    SetCutValue(cutForProton, "proton");
    SetCutValue(cutForProton, "anti_proton");

    if (verboseLevel>0) DumpCutValuesTable();
}


