// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: NTSTPrimaryGeneratorAction.cc,v 1.1 2003-11-07 21:30:29 japost Exp $ 
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "NTSTPrimaryGeneratorAction.hh"
#include "NTSTPrimaryGeneratorMessenger.hh"
#include "G4Event.hh"

//
// ctor
//
NTSTPrimaryGeneratorAction::NTSTPrimaryGeneratorAction()
{
  // ctor holds a pointer to the messenger which manages which type of 
  // generator to use
  messenger = new NTSTPrimaryGeneratorMessenger();
}

NTSTPrimaryGeneratorAction::~NTSTPrimaryGeneratorAction()
{
  delete messenger;
}

void
NTSTPrimaryGeneratorAction::GeneratePrimaries(G4Event* evt)
{
  //
  // the type of generator being used is determined by the messenger
  //
  messenger->GetGenerator()->GeneratePrimaryVertex( evt );
  //
  Print(evt);
}

void
NTSTPrimaryGeneratorAction::Print(const G4Event* evt) const
{
  //
  // print an entire event
  //

  if (messenger->PrintState()) {
    G4cout << "Generator \"" << messenger->GetName() << "\" used for " 
	   << " Event " << evt->GetEventID() << " : " << G4endl;
    for (int iv=0; iv < evt->GetNumberOfPrimaryVertex(); iv++) {
      evt->GetPrimaryVertex(iv)->Print();
    }
  }
}





