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
//
// $Id: Tst26EventAction.cc,v 1.3 2003-02-06 11:53:27 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
/////////////////////////////////////////////////////////////////////////
//
// test26: Cut per region physics
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst26EventAction.hh"

#include "Tst26RunAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26EventAction::Tst26EventAction(Tst26RunAction* run)
:Tst26Run(run),
 printModulo(100),
 Eth(1.0*keV)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26EventAction::~Tst26EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26EventAction::BeginOfEventAction(const G4Event* evt)
{
  G4int evtNb = evt->GetEventID();     
 
  E1    = 0.0;
  E9    = 0.0;
  E25   = 0.0;
  Eabs1 = 0.0;
  Eabs2 = 0.0;
  Eabs3 = 0.0;
  Eabs4 = 0.0;
  Evert.clear();
  Nvert.clear();

  //printing survey
  if (evtNb%printModulo == 0) 
    G4cout << "\n---> Begin of Event: " << evtNb << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26EventAction::EndOfEventAction(const G4Event* evt)
{  
  G4int n = Nvert.size();
  G4int nPad = 0;
  if (n > 0) {
    for(G4int i=0; i<n; i++) {
      if (Evert[i] > Eth) nPad++;
    }
  }
  Tst26Run->AddEvent(E1, E9, E25, Eabs1, Eabs2, Eabs3, Eabs4, nPad);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26EventAction::AddEnergy(G4double edep, G4int volIndex, G4int copyNo) 
{  
  if(0 == volIndex) {
    E25 += edep;
    if( (6<=copyNo&&8>=copyNo) || (11<=copyNo&&13>=copyNo) || 
        (16<=copyNo&&18>=copyNo)) {
      E9 += edep;
      if(12 == copyNo) E1 += edep;
    }
  } else if (1 == volIndex) {
    Eabs1 += edep;  
  } else if (2 == volIndex) {
    Eabs2 += edep;  
  } else if (3 == volIndex) {
    Eabs3 += edep;  
  } else if (4 == volIndex) {
    Eabs4 += edep;  
  } else if (5 == volIndex) {
    G4int n = Nvert.size();
    G4bool newPad = true;
    if (n > 0) {
      for(G4int i=0; i<n; i++) {
        if (copyNo == Nvert[i]) {
          newPad = false;
          Evert[i] += edep;
          break;
	}
      }
    }
    if(newPad) {
      Nvert.push_back(copyNo);
      Evert.push_back(edep);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


