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

#include "CCalSteppingAction.hh"
#include "CCalRunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4StepPoint.hh"
#include "G4ThreeVector.hh"
#include "G4RunManager.hh"

CCalSteppingAction::CCalSteppingAction()
{
  timeHistoMaxBin = 200;
  for (G4int i=0; i<200; i++) {timeDeposit[i] = 0.f;}
  for (G4int i=0; i<70;  i++) {LateralProfile[i] = 0.f;}
}

CCalSteppingAction::~CCalSteppingAction(){
}
  
void CCalSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4double de = aStep->GetTotalEnergyDeposit();
  if(de < CLHEP::eV) { return; }

  const G4StepPoint*  PostStepPoint= aStep->GetPostStepPoint(); 
  const G4StepPoint*  PreStepPoint= aStep->GetPreStepPoint(); 
  G4double time = 
    (PostStepPoint) ? PostStepPoint->GetGlobalTime()/nanosecond : 0.;

  G4int it = (G4int)time;
  it = std::min(it, 10000);
  //G4cout << "## time= " << time << " it= " << it << G4endl;
  G4int TSliceID = std::max(0,std::min(it,timeHistoMaxBin-1));

  G4float fde = (G4float)(de/CLHEP::GeV);   
  //G4cout << "  TSliceID= " << TSliceID << " de= " << fde << G4endl;

  timeDeposit[TSliceID] += fde;

  G4ThreeVector HitPoint = 0.5*(PostStepPoint->GetPosition()+
                                PreStepPoint->GetPosition());        
  // Because the beam axis has been defined as the x-axis, 
  // the lateral displacement is given in terms of the y and z positions. 
  G4double perp = 
    std::sqrt(HitPoint.y()*HitPoint.y()+HitPoint.z()*HitPoint.z())/cm;
  G4int ir = (G4int)perp;
  //G4cout << "  perp= " << perp << " ir= " << ir << G4endl;
  G4int radialPosition = std::max(0,std::min(69,ir));
  LateralProfile[radialPosition] += fde;
  //G4cout << "  done " << G4endl;
}

void CCalSteppingAction::endOfEvent(){

  G4AnalysisManager* man = G4AnalysisManager::Instance();
#ifdef debug
  G4double totalFilledProfileHcal = 0.0;
#endif
  static G4int IDlateralProfile = -1;
  if (IDlateralProfile < 0) {
    IDlateralProfile = man->GetH1Id("h500");
  }
  for (G4int i=0; i<70; ++i) {
    man->FillH1(IDlateralProfile+i,LateralProfile[i]);
#ifdef debug
    G4cout << "Fill Profile Hcal histo " << i << " with " << LateralProfile[i] << G4endl;
    totalFilledProfileHcal += LateralProfile[i];          
#endif
  }
 
#ifdef debug
  G4cout << "CCalAnalysis::InsertLateralProfile: Total filled Profile Hcal"
         << " histo " << totalFilledProfileHcal << G4endl;
#endif
    
  static G4int IDTimeHist = -1;
  if (IDTimeHist < 0) {
    IDTimeHist = man->GetH1Id("h300");
  }
  static G4int IDTimeProfile = -1;
  if (IDTimeProfile < 0) {
    IDTimeProfile = man->GetH1Id("h901");
  }
#ifdef debug
  G4double totalFilledTimeProfile = 0.0;
#endif
  for (G4int j=0; j<timeHistoMaxBin; ++j) {     
    man->FillH1(IDTimeHist+j,timeDeposit[j]);
#ifdef debug
    G4cout << "Fill Time slice histo " << j << " with " << timeDeposit[j] << G4endl;
    totalFilledTimeProfile += timeDeposit[j];
#endif
      
    G4double t = j + 0.5;
    man->FillH1(IDTimeProfile+1,t,timeDeposit[j]);
#ifdef debug
      G4cout << "Fill Time profile histo 1 with " << t << " " << timeDeposit[j] << G4endl;
#endif    
  }
#ifdef debug
  G4cout << "CCalAnalysis::InsertTime: Total filled Time profile histo " 
         << totalFilledTimeProfile << G4endl;
#endif

  for (G4int i=0; i<70; i++){LateralProfile[i] = 0.f;}
  for (G4int i=0; i<200; i++){timeDeposit[i] = 0.f;}

  G4cout << " --- End of event --- " << G4endl;
  
}  
