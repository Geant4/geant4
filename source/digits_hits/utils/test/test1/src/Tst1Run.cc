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
// $Id: Tst1Run.cc,v 1.1 2007-07-13 05:55:34 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst1Run.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"

Tst1Run::Tst1Run()
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

Tst1Run::~Tst1Run()
{;}

void Tst1Run::RecordEvent(const G4Event* evt)
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

G4double Tst1Run::GetTotal(const G4THitsMap<G4double> &map) const
{
  G4double tot = 0.;
  std::map<G4int,G4double*>::iterator itr = map.GetMap()->begin();
  for(; itr != map.GetMap()->end(); itr++) 
  { tot += *(itr->second); }
  return tot;
}


  
