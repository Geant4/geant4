//
//
//
// **********************************************************************
// *                                                                    *
// *                    ThyroidSteppingActionGEANT                      *
// *                                                                    *
// **********************************************************************
// 


#include "G4ios.hh"
#include "G4Track.hh"
#include "G4SteppingManager.hh"

#include "globals.hh"

#include <assert.h>
#include "g4std/fstream"
#include "g4std/iomanip"
#include "g4std/iostream"
#include "g4std/vector"

#include "ThyroidSteppingAction.hh"
//
//  Set Flag to write all the relevant parameters of the step (2)
//   or to write only Volume name, Energy loss, and Process Mechanism
// 
// G4int StoringType=1;                        //   dump tracking information
 G4int StoringType=2;                      //    Output selected parameters
 G4int VolumeCode;
 G4int ProcessCode;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ThyroidSteppingAction::ThyroidSteppingAction()

{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ThyroidSteppingAction::~ThyroidSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ThyroidSteppingAction::UserSteppingAction(const G4Step*)
 

{
   const G4SteppingManager* pSM = fpSteppingManager;
   G4Track* fTrack = pSM->GetTrack();
   G4Step* fStep = pSM->GetStep();
   G4int TrackID = fTrack->GetTrackID();
   G4int StepNo = fTrack->GetCurrentStepNumber();
   if(StoringType==1){
   if (StepNo==1){
 
    
     G4cout << G4std::setw( 5) << fTrack->GetCurrentStepNumber()-1 << "\t"
            << G4std::setw( 6) << fTrack-> GetDefinition()->GetParticleName()<< "\t"
	    << G4std::setw( 6) << fStep->GetPreStepPoint()->GetPosition().x()*mm << "\t "
	    << G4std::setw( 6) << fStep->GetPreStepPoint()->GetPosition().y()*mm << "\t "
	    << G4std::setw( 6) << fStep->GetPreStepPoint()->GetPosition().z()*mm << "\t "
	    << G4std::setw( 6) << fStep->GetPreStepPoint()->GetKineticEnergy() << "\t "
	    << G4std::setw( 6) << 0 << "\t "
	    << G4std::setw( 6) << 0 << "\t "
	    << G4std::setw( 6) << 0  << "\t "
            << G4std::setw( 6) <<  fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() << "\t"
            << G4std::setw( 6) << "initStep" << "\t";
     G4cout << G4endl;
   }
 

       
      G4cout << G4std::setw( 5) << fTrack->GetCurrentStepNumber() << "\t"
                << G4std::setw( 6) << fTrack-> GetDefinition()->GetParticleName()<< "\t"
             << G4std::setw( 6) << fTrack->GetPosition().x() << "\t "
             << G4std::setw( 6) << fTrack->GetPosition().y() << "\t "
             << G4std::setw( 6) << fTrack->GetPosition().z() << "\t "
             << G4std::setw( 6) << fTrack->GetKineticEnergy() << "  \t "
             << G4std::setw( 6) << fStep->GetTotalEnergyDeposit() << "\t "
             << G4std::setw( 6) << fStep->GetStepLength() << "\t "
             << G4std::setw( 6) << fTrack->GetTrackLength()  << "\t ";
     G4cout << G4std::setw(10) << fTrack->GetVolume()->GetName() << "\t";
     if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
     G4cout << " \t " 
            << G4std::setw(10) << fStep->GetPostStepPoint()->GetProcessDefinedStep()
  	                                ->GetProcessName();  }
     else { G4cout << "   UserLimit"; }
         G4cout << G4endl;} else 
       {if(fStep->GetTotalEnergyDeposit()>0.)
 {
   
     VolumeCode=10;

     
     if(fTrack->GetVolume()->GetName()=="UnThyroidPhys"){ 
         VolumeCode=1; 
         if(fTrack->GetPosition().y()>2.2){ VolumeCode=5;} 
      }
       if(fTrack->GetVolume()->GetName()=="TracheaPhys"){ VolumeCode=2;}
       if(fTrack->GetVolume()->GetName()=="IntTracheaPhys"){ VolumeCode=3;}
       if(fTrack->GetVolume()->GetName()=="WaterNeckPhys"){ VolumeCode=4;} 
       if(fTrack->GetVolume()->GetName()=="IodineNodulePhys"){ VolumeCode=6;}
    //
             ProcessCode=5;
   //
       if(fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="LowEnCompton"){ProcessCode=6;}
       if(fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="LowEnergyIoni"){ProcessCode=7;}
       if(fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="LowEnPhotoElec"){ProcessCode=8;}
       if(fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="LowEnBrem"){ProcessCode=9;}
   //        
      G4cout << G4std::setw(10) << VolumeCode << "\t "              
             << G4std::setw( 6) << fStep->GetTotalEnergyDeposit() << "\t "
             << G4std::setw(10) << ProcessCode << "\t";
      G4cout << G4endl;
    }     
     }
}




































