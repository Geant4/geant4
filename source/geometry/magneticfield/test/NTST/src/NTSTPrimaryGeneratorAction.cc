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
// $Id: NTSTPrimaryGeneratorAction.cc,v 1.3 2006-06-29 18:26:27 gunter Exp $ 
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





