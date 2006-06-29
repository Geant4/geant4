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
#include "CCalSteppingAction.hh"

#include "G4SDManager.hh"
#include "G4StepPoint.hh"
#include <iostream>
#include "G4ThreeVector.hh"
#include <cmath>

#ifdef G4ANALYSIS_USE  
#include "CCalAnalysis.hh"
#endif


CCalSteppingAction::CCalSteppingAction(){

#ifdef G4ANALYSIS_USE  
  CCalAnalysis* analysis = CCalAnalysis::getInstance();
  timeHistoMaxBin=analysis->maxbin();
#else
  timeHistoMaxBin=1;
#endif

  int i; 
  for (i=0; i<200; i++) {timeDeposit[i] = 0.;}
  for (i=0; i<70;  i++) {LateralProfile[i] = 0.;}

}


CCalSteppingAction::~CCalSteppingAction(){
  G4cout << "CCalSteppingAction deleted" << G4endl;
}
  

void CCalSteppingAction::UserSteppingAction(const G4Step* aStep){

  G4StepPoint*  PostStepPoint= aStep->GetPostStepPoint(); 
  G4StepPoint*  PreStepPoint= aStep->GetPreStepPoint(); 
  int TSliceID;

  TSliceID = static_cast<int>( (PostStepPoint->GetGlobalTime() ) / nanosecond);
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

#ifdef G4ANALYSIS_USE  
  CCalAnalysis* analysis = CCalAnalysis::getInstance();
  analysis->InsertLateralProfile(LateralProfile);  
  analysis->InsertTime(timeDeposit); 
#endif
  
  int i=0;
  for (i=0; i<70; i++){LateralProfile[i] = 0.;}
  for (i=0; i<200; i++){timeDeposit[i] = 0.;}
  
}  
