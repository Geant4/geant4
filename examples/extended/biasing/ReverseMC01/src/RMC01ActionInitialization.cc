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
/// \file biasing/ReverseMC01/src/RMC01ActionInitilaization.cc
/// \brief Definition of the RMC01ActionInitialization class
//
//
//////////////////////////////////////////////////////////////
//  Class Name: RMC01ActionInitilaization
//        Author: L. Desorgher
//        Organisation: Radiation Physics Institute, Lausanne University Hospital
//        Customer: ESA/ESTEC
//////////////////////////////////////////////////////////////
//  ChangeHistory:
//        14-10-2025 creation by L. Desorgher
//
//-------------------------------------------------------------
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RMC01ActionInitialization.hh"
#include "RMC01ActionInitialization.hh"
#include "RMC01PrimaryGeneratorAction.hh"
#include "RMC01EventAction.hh"
#include "RMC01RunAction.hh"

#include "RMC01AnalysisManager.hh"

//Adjoint Simulation
#include "G4AdjointSimManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01ActionInitialization::RMC01ActionInitialization()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01ActionInitialization::~RMC01ActionInitialization()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01ActionInitialization::BuildForMaster() const
{
//  RMC01AnalysisManager::GetInstance();

  RMC01RunAction* theRunAction = new RMC01RunAction();
  SetUserAction(theRunAction);
  //Adjoint Simulation
  G4AdjointSimManager* theAdjointSimManager = G4AdjointSimManager::GetInstance();
  theAdjointSimManager->SetAdjointRunAction(theRunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01ActionInitialization::Build() const
{
  SetUserAction(new RMC01PrimaryGeneratorAction());

  RMC01RunAction* theRunAction = new RMC01RunAction();
  SetUserAction(theRunAction);
  //SetUserAction(new RMC01RunAction);

  //G4UserEventAction *myUserEventAction = new RMC01EventAction();
  RMC01EventAction*  theEventAction = new RMC01EventAction();
  SetUserAction(theEventAction);

  //Adjoint Simulation
  G4AdjointSimManager *theAdjointSimManager = G4AdjointSimManager::GetInstance();
  theAdjointSimManager->SetAdjointRunAction(theRunAction);
  theAdjointSimManager->SetAdjointEventAction(theEventAction);
}
