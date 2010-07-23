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
// $Id: Test2Run.cc,v 1.1 2010-07-23 06:15:41 akimura Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Test2Run.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"

Test2Run::Test2Run()
{
  G4String detName = "MassWorld";
  G4String primNameSum[7] = {"eDep","trackLengthGamma","nStepGamma",
                   "trackLengthElec","nStepElec","trackLengthPosi","nStepPosi"};
  G4SDManager* SDMan = G4SDManager::GetSDMpointer();
  G4String fullName;
  for(G4int i=0;i<7;i++)
  {
    fullName = detName+"/"+primNameSum[i];
    colIDSum[i] = SDMan->GetCollectionID(fullName);
  }
}

Test2Run::~Test2Run()
{;}

void Test2Run::RecordEvent(const G4Event* evt)
{
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
  numberOfEvent++;
  for(G4int i=0;i<7;i++)
  {
    G4THitsMap<G4double>* evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(colIDSum[i]));
    mapSum[i] += *evtMap;
  }
}

G4double Test2Run::GetTotal(const G4THitsMap<G4double> &map) const
{
  G4double tot = 0.;
  std::map<G4int,G4double*>::iterator itr = map.GetMap()->begin();
  for(; itr != map.GetMap()->end(); itr++) 
  { tot += *(itr->second); }
  return tot;
}


  
