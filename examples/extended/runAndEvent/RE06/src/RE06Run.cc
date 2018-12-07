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
/// \file RE06/src/RE06Run.cc
/// \brief Implementation of the RE06Run class
//
//

#include "RE06Run.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RE06Run::RE06Run()
 : G4Run()
{
  G4String detName[6] 
    = {"Calor-A_abs","Calor-A_gap",
       "Calor-B_abs","Calor-B_gap",
       "Calor-C_abs","Calor-C_gap"};
       
  G4String primNameSum[6] 
    = {"eDep","nGamma","nElectron","nPositron","trackLength","nStep"};
    
  G4String primNameMin[3] 
    = {"minEkinGamma","minEkinElectron","minEkinPositron"};

  G4String paraName[3] 
    = {"Calor-AP_para","Calor-BP_para","Calor-CP_para"};

  G4SDManager* SDMan = G4SDManager::GetSDMpointer();
  G4String fullName;
  size_t i,j;
  for(i=0;i<6;i++)
  {
    for(j=0;j<6;j++)
    {
      fullName = detName[i]+"/"+primNameSum[j];
      fColIDSum[i][j] = SDMan->GetCollectionID(fullName);
    }
    for(j=0;j<3;j++)
    {
      fullName = detName[i]+"/"+primNameMin[j];
      fColIDMin[i][j] = SDMan->GetCollectionID(fullName);
    }
  }
  for(i=0;i<3;i++)
  {
    for(j=0;j<6;j++)
    {
      fullName = paraName[i]+"/"+primNameSum[j];
      fColIDPara[i][j] = SDMan->GetCollectionID(fullName);
    }
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RE06Run::~RE06Run()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE06Run::RecordEvent(const G4Event* evt)
{
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
  numberOfEvent++;
  size_t i,j;
  for(i=0;i<6;i++)
  {
    for(j=0;j<6;j++)
    {
      G4THitsMap<G4double>* evtMap 
        = (G4THitsMap<G4double>*)(HCE->GetHC(fColIDSum[i][j]));
      fMapSum[i][j] += *evtMap;
    }

    for(j=0;j<3;j++)
    {
      G4THitsMap<G4double>* evtMap 
        = (G4THitsMap<G4double>*)(HCE->GetHC(fColIDMin[i][j]));
      std::map<G4int,G4double*>::iterator itr = evtMap->GetMap()->begin();
      for(; itr != evtMap->GetMap()->end(); itr++)
      {
        G4int key = (itr->first);
        G4double val = *(itr->second);
        G4double* mapP = fMapMin[i][j][key];
        if( mapP && (val>*mapP) ) continue;
        fMapMin[i][j].set(key,val);
      }
    }

  }
  for(i=0;i<3;i++)
  {
    for(j=0;j<6;j++)
    {
      G4THitsMap<G4double>* evtMap 
        = (G4THitsMap<G4double>*)(HCE->GetHC(fColIDPara[i][j]));
      fMapPara[i][j] += *evtMap;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE06Run::Merge(const G4Run * aRun) {
  const RE06Run * localRun = static_cast<const RE06Run *>(aRun);
  
  for(G4int i = 0; i < 6; i++) {
    for(G4int j = 0; j < 6; j++) {
      fMapSum[i][j] += localRun->fMapSum[i][j];
    }

    for(G4int j = 0; j < 3; j++) {
      std::map<G4int, G4double*>::iterator itr = localRun->fMapMin[i][j].GetMap()->begin();
      for(; itr != localRun->fMapMin[i][j].GetMap()->end(); itr++) {
        G4int key = itr->first;
        G4double val = *(itr->second);
        G4double * mapP = fMapMin[i][j][key];
        if(!mapP || val < *mapP) fMapMin[i][j].set(key, val);
      }
    }
  }

  for(G4int i = 0; i < 3; i++) {
      for(G4int j = 0; j < 6; j++) {
        fMapPara[i][j] += localRun->fMapPara[i][j];
      }
  }
  
  G4Run::Merge(aRun);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double RE06Run::GetTotal(const G4THitsMap<G4double> &map) const
{
  G4double tot = 0.;
  if(map.GetSize()==0) return tot;
  std::map<G4int,G4double*>::iterator itr = map.GetMap()->begin();
  for(; itr != map.GetMap()->end(); itr++) 
  { tot += *(itr->second); }
  return tot;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double RE06Run::FindMinimum(const G4THitsMap<G4double> &map) const
{
  G4double val = DBL_MAX;

  if(map.GetSize()==0) return val;
  std::map<G4int,G4double*>::iterator itr = map.GetMap()->begin();
  for(; itr != map.GetMap()->end(); itr++) 
  { if(val>*(itr->second)) val = *(itr->second); }
  return val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
