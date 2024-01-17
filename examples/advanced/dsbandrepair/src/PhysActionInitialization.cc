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
/// \file PhysActionInitialization.hh
/// \brief Definition of the PhysActionInitialization class

#include "PhysActionInitialization.hh"

#include "PhysPrimaryGeneratorAction.hh"
#include "PhysEventAction.hh"
#include "PhysRunAction.hh"
#include "PhysSteppingAction.hh"
#include "PhysChemIO.hh"
#include "G4DNAChemistryManager.hh"
#include "G4Threading.hh"

#include <memory>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysActionInitialization::BuildForMaster() const
{
    SetUserAction(new PhysRunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysActionInitialization::Build() const
{
    PhysPrimaryGeneratorAction* primGenAction = new PhysPrimaryGeneratorAction();
    SetUserAction(primGenAction);

    PhysEventAction* eventAction = new PhysEventAction;
    SetUserAction(eventAction);

    SetUserAction(new PhysRunAction);

    PhysSteppingAction* steppingAction = new PhysSteppingAction(eventAction);
    SetUserAction(steppingAction);
    //pass- PhysChemIO to G4DNAChemistryManager
    std::unique_ptr<G4VPhysChemIO> fPhysChemIO = std::make_unique<PhysChemIO>(steppingAction);
    G4DNAChemistryManager::Instance()->SetPhysChemIO(std::move(fPhysChemIO));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......