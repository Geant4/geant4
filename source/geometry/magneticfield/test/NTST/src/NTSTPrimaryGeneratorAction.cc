//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: NTSTPrimaryGeneratorAction.cc,v 1.2 2003-12-09 15:35:40 gunter Exp $ 
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





