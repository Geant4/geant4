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
// $Id: Tst50SteppingAction.cc,v 1.1 2002-11-26 17:57:49 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst50SteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "Tst50EventAction.hh"
#include <math.h> // standard c math library
#include "Tst50AnalysisManager.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst50SteppingAction::Tst50SteppingAction(Tst50EventAction* EA):
IDold(-1),eventaction (EA) 
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst50SteppingAction::UserSteppingAction(const G4Step* Step)
{ Tst50AnalysisManager* analysis = Tst50AnalysisManager::getInstance();

G4int evno = eventaction->GetEventno() ;//mi dice a che evento siamo 
 
  //IDnow identifica univocamente la particella//
IDnow = evno+10000*(Step->GetTrack()->GetTrackID())+
          100000000*(Step->GetTrack()->GetParentID()); 
  //G4String CurV;
  //G4String NexV;

    G4String PTyp;

    G4double KinE;

//    G4double EDep; 
// commented beacause in class track method
// GetTotalEnergyDeposit in declared but not implemented 

    G4double XPos;
    G4double YPos;
    G4double ZPos;
    
    G4double XMoD;
    G4double YMoD;
    G4double ZMoD;

    G4double ThMD;
    G4double PhMD;
    G4double XYRd;

    //CurV = Step->GetTrack()->GetVolume()->GetName();
    //NexV = Step->GetTrack()->GetNextVolume()->GetName();

    PTyp = Step->GetTrack()->GetDefinition()->GetParticleType();    
   	
    KinE = Step->GetTrack()->GetKineticEnergy();

//    EDep = Step->GetTrack()->GetTotalEnergyDeposit();

    XPos = Step->GetTrack()->GetPosition().x();
    YPos = Step->GetTrack()->GetPosition().y();
    ZPos = Step->GetTrack()->GetPosition().z();

    XMoD = Step->GetTrack()->GetMomentumDirection().x();
    YMoD = Step->GetTrack()->GetMomentumDirection().y();
    ZMoD = Step->GetTrack()->GetMomentumDirection().z();

    if( XMoD != 0){
        ThMD = atan(YMoD/XMoD);
    } else {
        ThMD = YMoD/fabs(YMoD)*3.1428; 
    }

    XYRd = sqrt(pow(XMoD,2.0)+pow(YMoD,2.0));

    if (XYRd != 0) {
        PhMD = atan(ZMoD/XYRd);
    } else {
        PhMD = ZMoD/fabs(ZMoD)*3.1428;
    }
  
    //susanna, istogrammo i processi delle particelle primarie// 
 
        
      if(0 == Step->GetTrack()->GetParentID() ) //primary particle 
	{ 
         G4String process=Step->GetPostStepPoint()->GetProcessDefinedStep()
	   ->GetProcessName();

          G4cout<< "processo del fotone-------------------------"
                      << process<<G4endl;
          if (process=="LowEnPhotoElec") 
	    { G4cout<< "processo del fotone-------------------------"
                      <<" 1"<<G4endl; 
	  analysis->primary_processes(3);
	}
      if (process=="Transportation") 
           { 
	     analysis->primary_processes(1);
	}
 if (process=="LowEnRayleigh") 
           { 
	     analysis->primary_processes(2);
 G4cout<< "processo del fotone-------------------------"
                      << 2<<G4endl;
	}
 if (process=="LowEnCompton") 
           { 
	     analysis->primary_processes(4);
	   }

if (process=="LowEnConversion") 
           { 
	     analysis->primary_processes(5);
	   }


	}

	   

    /*
    G4cout << "UserSteppingAction:"
//	   << " CurV " << CurV
//	   << " NexV " << NexV
      //   << " PTyp " << PTyp  
	   << " KinE " << KinE
//         << " EDep " << EDep 
      //<< " XPox " << XPos 
      //	   << " YPos " << YPos 
      //	   << " ZPos " << ZPos 	   
	   << " XMoD " << XMoD 
	   << " YMoD " << YMoD 
	   << " ZMoD " << ZMoD 
           << " ThMD " << ThMD
           << " PhMD " << PhMD
	   << G4endl;
    */
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......













