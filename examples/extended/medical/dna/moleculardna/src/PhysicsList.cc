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
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class

#include "PhysicsList.hh"
#include "G4SystemOfUnits.hh"
#include "ParallelWorldPhysics.hh"
#include "G4EmDNAPhysics.hh"
#include "G4EmDNAPhysics_option2.hh"
#include "G4EmDNAPhysics_option4.hh"
#include "G4EmDNAPhysics_option6.hh"

#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

#include "ChemistryList.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList(const G4int& phylist)
  : G4VModularPhysicsList()
{
  SetDefaultCutValue(1.0 * micrometer);
  SetVerboseLevel(1);

  UseCoupledTransportation();

  if(phylist == 0)
  {
    RegisterPhysics(new G4EmDNAPhysics());
  }
  else if(phylist == 2)
  {
    RegisterPhysics(new G4EmDNAPhysics_option2());
  }
  else if(phylist == 4)
  {
    RegisterPhysics(new G4EmDNAPhysics_option4());
  }
  else if(phylist == 6)
  {
    RegisterPhysics(new G4EmDNAPhysics_option6());
  }
  else
  {
    G4ExceptionDescription errmsg;
    errmsg << "Recommend to option 2, 4, 6 or default" << G4endl;
    G4Exception("PhysicsList::PhysicsList", "", FatalException, errmsg);
  }
  RegisterPhysics(new G4DecayPhysics());
  RegisterPhysics(new G4RadioactiveDecayPhysics());

  G4bool useParallelPhysicsWorld = false;
  if (useParallelPhysicsWorld) {
    RegisterPhysics(new ParallelWorldPhysics("DNAWorld", true));
  }

  RegisterPhysics(new ChemistryList());

  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100 * eV,
                                                                  1 * GeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
