#include "HcalTB96SteppingAction.hh"
#include "HcalTB96Analysis.hh"
#include "G4CaloSD.hh"

#include "G4SDManager.hh"
#include "G4StepPoint.hh"
#include <iostream>
#include "G4ThreeVector.hh"

HcalTB96SteppingAction::HcalTB96SteppingAction(){

  HcalTB96Analysis* analysis = HcalTB96Analysis::getInstance();
  timeHistoMaxBin=analysis->maxbin();
  int i; 
  for (i=0; i<50; i++){timeDeposit[i] = 0.;}
  for (i=0; i<28; i++){LateralProfile[i] = 0.;}

}

HcalTB96SteppingAction::~HcalTB96SteppingAction(){

  cout <<"Deleting HcalTB96SteppingAction"<<endl;
}
  

void HcalTB96SteppingAction::UserSteppingAction(const G4Step* aStep){

  G4StepPoint*  PostStepPoint= aStep->GetPostStepPoint(); 
  G4StepPoint*  PreStepPoint= aStep->GetPreStepPoint(); 
  int TSliceID;

  TSliceID = (int) (PostStepPoint->GetGlobalTime() )*nanosecond;
  TSliceID = TSliceID<timeHistoMaxBin ? TSliceID : timeHistoMaxBin-1;
  timeDeposit[TSliceID] += aStep->GetTotalEnergyDeposit();


  G4ThreeVector HitPoint = 0.5*(PostStepPoint->GetPosition()+
				PreStepPoint->GetPosition());	
  int radialPosition = min(27,int(HitPoint.perp()/cm));
  LateralProfile[radialPosition] += aStep->GetTotalEnergyDeposit();
  
}


void HcalTB96SteppingAction::endOfEvent(){

  HcalTB96Analysis* analysis = HcalTB96Analysis::getInstance();
  analysis->InsertLateralProfile(LateralProfile);  
  analysis->InsertTime(timeDeposit); 
  
  int i=0;
  for (i=0; i<28; i++){LateralProfile[i] = 0.;}
  for (i=0; i<50; i++){timeDeposit[i] = 0.;}
  
}  
