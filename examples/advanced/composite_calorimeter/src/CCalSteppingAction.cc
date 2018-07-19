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
///////////////////////////////////////////////////////////////////////////////
// File: CCalSeppingAction.cc
// Description: Study profiling during the steps
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cmath>

#include "CCalSteppingAction.hh"
#include "CCalRunAction.hh"
#include "CCalAnalysis.hh"

#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4StepPoint.hh"
#include "G4ThreeVector.hh"
#include "G4RunManager.hh"

#include "CCalAnalysis.hh"

CCalSteppingAction::CCalSteppingAction() : 
  runAct(nullptr)
{
  for (int i=0; i<200; i++) {timeDeposit[i] = 0.;}
  for (int i=0; i<70;  i++) {LateralProfile[i] = 0.;}

}


CCalSteppingAction::~CCalSteppingAction(){
  G4cout << "CCalSteppingAction deleted" << G4endl;
}
  

void CCalSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  //thread-local run action
  if (!runAct) 
    runAct = 
      dynamic_cast<const CCalRunAction*>
      (G4RunManager::GetRunManager()->GetUserRunAction());

  
  timeHistoMaxBin=runAct->maxbin();
  G4StepPoint*  PostStepPoint= aStep->GetPostStepPoint(); 
  G4StepPoint*  PreStepPoint= aStep->GetPreStepPoint(); 
  int TSliceID;

  if ( PostStepPoint->GetGlobalTime() / nanosecond > 1.0E9 ) TSliceID = 999999999;
  else TSliceID = static_cast<int>( PostStepPoint->GetGlobalTime() / nanosecond );
  TSliceID = TSliceID<timeHistoMaxBin ? TSliceID : timeHistoMaxBin-1;
  timeDeposit[TSliceID] += aStep->GetTotalEnergyDeposit() / GeV;

  G4ThreeVector HitPoint = 0.5*(PostStepPoint->GetPosition()+
				PreStepPoint->GetPosition());	
  // Because the beam axis has been defined as the x-axis, 
  // the lateral displacement is given in terms of the y and z positions. 
  double perp = std::sqrt(HitPoint.y()*HitPoint.y()+HitPoint.z()*HitPoint.z());
  int radialPosition = std::min(69,int(perp/cm));
  LateralProfile[radialPosition] += aStep->GetTotalEnergyDeposit() / GeV;
  
}


void CCalSteppingAction::endOfEvent(){

  G4AnalysisManager* man = G4AnalysisManager::Instance();
   G4double totalFilledProfileHcal = 0.0;

   static G4int IDlateralProfile = -1;
    if (IDlateralProfile < 0)
      IDlateralProfile = man->GetH1Id("h500");

  for (G4int i=0; i<70; i++) {    
    man->FillH1(IDlateralProfile+i,LateralProfile[i]);
#ifdef debug
    G4cout << "Fill Profile Hcal histo " << i << " with " << LateralProfile[i] << G4endl;
#endif
    totalFilledProfileHcal += LateralProfile[i];          
  }
 
#ifdef debug
    G4cout << "CCalAnalysis::InsertLateralProfile: Total filled Profile Hcal"
	   << " histo " << totalFilledProfileHcal << G4endl;
#endif
    
    static G4int IDTimeHist = -1;
    if (IDTimeHist < 0)
      IDTimeHist = man->GetH1Id("h300");
    G4double totalFilledTimeProfile = 0.0;
    for (G4int j=0; j<timeHistoMaxBin; j++) 
    {     
      man->FillH1(IDTimeHist+j,timeDeposit[j]);
#ifdef debug
      G4cout << "Fill Time slice histo " << j << " with " << timeDeposit[j] << G4endl;
#endif
      totalFilledTimeProfile += timeDeposit[j];
      
      static G4int IDTimeProfile = -1;
      if (IDTimeProfile < 0)
	IDTimeProfile = man->GetH1Id("h901");
      G4double t = j + 0.5;
      man->FillH1(IDTimeProfile+1,t,timeDeposit[j]);
#ifdef debug
      G4cout << "Fill Time profile histo 1 with " << t << " " << x << G4endl;
#endif    
    }
  #ifdef debug
    G4cout << "CCalAnalysis::InsertTime: Total filled Time profile histo " 
	   << totalFilledTimeProfile << G4endl;
#endif

  int i=0;
  for (i=0; i<70; i++){LateralProfile[i] = 0.;}
  for (i=0; i<200; i++){timeDeposit[i] = 0.;}

  G4cout << " --- End of event --- " << G4endl;
  
}  
