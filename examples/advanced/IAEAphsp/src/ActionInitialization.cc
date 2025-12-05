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

#include "ActionInitialization.hh"
#include "ActionInitializationMessenger.hh"

#include "globals.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "G4RunManager.hh"
#include "G4Threading.hh"

#include "G4IAEAphspReader.hh"
#include "G4IAEAphspWriterStack.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::ActionInitialization()
  : G4VUserActionInitialization()
{
  // Messenger
  fMessenger = new ActionInitializationMessenger(this);

  // IAEAphsp source file name (including path) for primary generator
  fIAEAphspReaderName = "";
  fNumberOfThreads = G4RunManager::GetRunManager()->GetNumberOfThreads();

  // Name prefix, including path, of IAEAphsp output files (default, nothing).
  fIAEAphspWriterNamePrefix = "";

  // Vector to register phsp planes (Z=const)
  fZphspVec = new std::vector<G4double>;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::~ActionInitialization()
{
  if (fZphspVec) {
    fZphspVec->clear();
    delete fZphspVec;
  }
  if (fMessenger)           delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildForMaster() const
{
  // G4cout << "ActionInitialization::BuildForMaster() started" << G4endl;
  SetUserAction(new RunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::Build() const
{
  // G4cout << "ActionInitialization::Build() started" << G4endl;
  G4cout << "IAEAphsp file to read is \"" << fIAEAphspReaderName << "\"."
	 << G4endl;

  PrimaryGeneratorAction* prim = new PrimaryGeneratorAction(fNumberOfThreads);
  if ( !fIAEAphspReaderName.empty() )   // never true in sequential mode
    prim->SetIAEAphspReader(fIAEAphspReaderName);
  SetUserAction(prim);

  RunAction* runAct = new RunAction();
  if (fIAEAphspWriterNamePrefix != "") {  // never true in sequential mode
    // Set G4IAEAphspWriterStack object for the local thread
    // and register zphsp values to it
    runAct->SetIAEAphspWriterStack(fIAEAphspWriterNamePrefix);

    if (fZphspVec->size() > 0) {
      for (const auto& zphsp : (*fZphspVec))
	runAct->AddZphsp(zphsp);
    }
    else {
      G4ExceptionDescription msg;
      msg << "IAEAphsp output file name provided, but no zphsp values "
	  << "have been registered!" << G4endl;
      G4Exception("ActionInitialization::Build()",
		  "ActionInit001", FatalException, msg);
    }
  }
  SetUserAction(runAct);

  SteppingAction *steppingAction = new SteppingAction();
  SetUserAction(steppingAction);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::SetIAEAphspReader(const G4String& name)
{
  fIAEAphspReaderName = name;

  if ( !(G4Threading::IsMultithreadedApplication()) ) {
    // In sequential mode, when this command is issued, Build() has been
    // called already. Thus, we must set G4IAEAphspReader object here

    const G4VUserPrimaryGeneratorAction* basePrim =
      G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
    if (!basePrim) return; // No PGA yet (very unlikely in sequential)

    // 1) cast while preserving constness
    const auto* myConstPrim =
      dynamic_cast<const PrimaryGeneratorAction*>(basePrim);
    if (!myConstPrim) return;

    // 2) Drop constness to set G4IAEAphspReader object
    auto* myPrim = const_cast<PrimaryGeneratorAction*>(myConstPrim);
    myPrim->SetIAEAphspReader(name);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::SetIAEAphspWriterPrefix(const G4String& prefix)
{
  fIAEAphspWriterNamePrefix = prefix;

  if ( !(G4Threading::IsMultithreadedApplication()) ) {
    // In sequential mode, when this command is issued, Build() has been
    // called already. Thus, we must set G4IAEAphspWriterStack object here
    const G4UserRunAction* baseRA =
      G4RunManager::GetRunManager()->GetUserRunAction();
    if (!baseRA) return; // No run action defined

    // 1) cast while preserving constness
    const auto* myConstRA = dynamic_cast<const RunAction*>(baseRA);
    if (!myConstRA) return;

    // 2) Drop constness to modify RunAction object status
    auto* myRA = const_cast<RunAction*>(myConstRA);
    myRA->SetIAEAphspWriterStack(prefix);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::AddZphsp(const G4double zphsp)
{
  fZphspVec->push_back(zphsp);

  if ( !(G4Threading::IsMultithreadedApplication()) ) {
    // In sequential mode, when this command is issued, Build() has been
    // called already. Thus, we must set G4IAEAphspWriterStack object here
    const G4UserRunAction* baseRA =
      G4RunManager::GetRunManager()->GetUserRunAction();
    if (!baseRA) return; // No run action defined

    // 1) cast while preserving constness
    const auto* myConstRA = dynamic_cast<const RunAction*>(baseRA);
    if (!myConstRA) return;

    // 2) Drop constness to modify RunAction object status
    auto* myRA = const_cast<RunAction*>(myConstRA);
    myRA->AddZphsp(zphsp);
  }
}
