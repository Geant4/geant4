///////////////////////////////////////////////////////////////////////////////
// File: CCalSeppingAction.cc
// Description: Study profiling during the steps
///////////////////////////////////////////////////////////////////////////////
#include "CCalSteppingAction.hh"

#include "G4SDManager.hh"
#include "G4StepPoint.hh"
#include "g4std/iostream"
#include "G4ThreeVector.hh"
#include <math.h>

#ifdef G4ANALYSIS_USE  
#include "CCalAnalysis.hh"
#endif


CCalSteppingAction::CCalSteppingAction(){

#ifdef G4ANALYSIS_USE  
  CCalAnalysis* analysis = CCalAnalysis::getInstance();
  timeHistoMaxBin=analysis->maxbin();
#endif

  int i; 
  for (i=0; i<40; i++){timeDeposit[i] = 0.;}
  for (i=0; i<70; i++){LateralProfile[i] = 0.;}

}


CCalSteppingAction::~CCalSteppingAction(){
  G4cout <<"Deleting CCalSteppingAction"<<G4endl;
}
  

void CCalSteppingAction::UserSteppingAction(const G4Step* aStep){

  G4StepPoint*  PostStepPoint= aStep->GetPostStepPoint(); 
  G4StepPoint*  PreStepPoint= aStep->GetPreStepPoint(); 
  int TSliceID;

  TSliceID = (int) (PostStepPoint->GetGlobalTime() ) / nanosecond;
  TSliceID = TSliceID<timeHistoMaxBin ? TSliceID : timeHistoMaxBin-1;
  timeDeposit[TSliceID] += aStep->GetTotalEnergyDeposit() / GeV;

  G4ThreeVector HitPoint = 0.5*(PostStepPoint->GetPosition()+
				PreStepPoint->GetPosition());	
  // Because the beam axis has been defined as the x-axis, 
  // the lateral displacement is given in terms of the y and z positions. 
  double perp = sqrt(HitPoint.y()*HitPoint.y()+HitPoint.z()*HitPoint.z());
  int radialPosition = G4std::min(69,int(perp/cm));
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
  for (i=0; i<40; i++){timeDeposit[i] = 0.;}
  
}  
