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
// $Id: Tst50SteppingAction.cc,v 1.31 2003-05-17 11:32:54 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"
#include <math.h> // standard c math library
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "Tst50SteppingAction.hh"
#include "Tst50DetectorConstruction.hh"
#include "Tst50EventAction.hh"
#ifdef G4ANALYSIS_USE
#include "Tst50AnalysisManager.hh"
#endif
#include "Tst50PrimaryGeneratorAction.hh"
#include "Tst50RunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst50SteppingAction::Tst50SteppingAction(Tst50EventAction* EA, Tst50PrimaryGeneratorAction* PG, Tst50RunAction* RA, Tst50DetectorConstruction* DET):
  IDold(-1),eventaction (EA), p_Primary(PG), runaction(RA), detector(DET) 
{ 

  
}

Tst50SteppingAction::~Tst50SteppingAction()

{;}

void Tst50SteppingAction::UserSteppingAction(const G4Step* Step)
{
#ifdef G4ANALYSIS_USE
Tst50AnalysisManager* analysis = Tst50AnalysisManager::getInstance();
#endif
  
  G4int evno = eventaction->GetEventno() ;
  G4int run_ID= runaction-> GetRun_ID();
  G4bool flag=runaction-> Get_flag();

  G4double  initial_energy= p_Primary->GetInitialEnergy();
 
  G4int IDnow = run_ID+1000*evno+10000*(Step->GetTrack()->GetTrackID())+
          100000000*(Step->GetTrack()->GetParentID()); 

  G4double  KinE = Step->GetTrack()->GetKineticEnergy();
  G4String PrimaryName= p_Primary->GetParticle(); 
  G4double  ZMoD = Step->GetTrack()->GetMomentumDirection().z();
    
// ---------- energy of primary transmitted particles-----------// 

  if(PrimaryName=="e-"|| PrimaryName=="e+"|| PrimaryName=="proton" ||
     PrimaryName=="gamma" )
    { 
      if(Step->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Target")
	{
	  if(Step->GetTrack()->GetNextVolume()->GetName() == "World" ) 
	    { 
	      if(0 == Step->GetTrack()->GetParentID() ) 
		{ if(IDnow != IDold) 
		  {
		    IDold=IDnow ;  
		    if (ZMoD >= 0. )
		      {
			runaction->Trans_number();
		      }
		    else
		      {
			runaction->Back_number(); 
		      }
		    
		    if( (initial_energy == KinE) &&( ZMoD==1.))
		      {
			runaction->gamma_transmitted();
		      }
		  }
		}
	    }
	}
    }

// Stoppint Power test

  if (flag==false)
   {
     if(PrimaryName=="e-" || PrimaryName=="proton" || PrimaryName=="e+")
     {
      if(0 ==Step->GetTrack()->GetParentID() ) 
       {
        if(IDnow != IDold)
         { 
          IDold=IDnow ; 
          G4double  XPos_pre=Step->GetPreStepPoint()->GetPosition().x();
          G4double YPos_pre=Step->GetPreStepPoint()->GetPosition().y();
          G4double ZPos_pre=Step->GetPreStepPoint()->GetPosition().z();
          if (0== XPos_pre && 0== YPos_pre && 0== ZPos_pre)
           { 
            G4double energyLost=abs(Step->GetDeltaEnergy());
            G4double stepLength= Step->GetTrack()-> GetStepLength();
            if(stepLength!=0) 
             {
              G4double TotalStoppingPower=(energyLost/stepLength);
              G4double SP=TotalStoppingPower/(detector->GetDensity());
              //#ifdef ANALYSIS_USE
	      analysis->StoppingPower(run_ID,initial_energy/MeV,SP/(MeV*(cm2/g)));
	      //#endif
                                   
              }
	     }
	   }
       
      //CSDA range
         range=0.; 
         if(0.0 == KinE) 
          {  
	   G4double  xend= Step->GetTrack()->GetPosition().x()/mm ;
           G4double yend= Step->GetTrack()->GetPosition().y()/mm ;
           G4double  zend= Step->GetTrack()->GetPosition().z()/mm ;
           range=(sqrt(xend*xend+yend*yend+zend*zend)); 
           G4double range2= range*(detector->GetDensity());
          
           analysis->CSDARange(run_ID,initial_energy/MeV,range2/(g/cm2));
	  }
       }
     }
   }
}


