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
// $Id: Tst50SteppingAction.cc,v 1.6 2003-01-07 15:29:40 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst50SteppingAction.hh"
#include "Tst50DetectorConstruction.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "Tst50EventAction.hh"
#include <math.h> // standard c math library
#ifdef G4ANALYSIS_USE
#include "Tst50AnalysisManager.hh"
#endif
#include "Tst50PrimaryGeneratorAction.hh"
#include "Tst50RunAction.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst50SteppingAction::Tst50SteppingAction(Tst50EventAction* EA, Tst50PrimaryGeneratorAction* PG, Tst50RunAction* RA, Tst50DetectorConstruction* DET ):
  IDold(-1),eventaction (EA), p_Primary(PG), runaction(RA), detector(DET) 
{ 

  //initial_energy= p_Primary->GetInitialEnergy();
  //G4cout<<initial_energy<<"form stepping"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Tst50SteppingAction::~Tst50SteppingAction()

{;}

void Tst50SteppingAction::UserSteppingAction(const G4Step* Step)
{ 
#ifdef G4ANALYSIS_USE
Tst50AnalysisManager* analysis = Tst50AnalysisManager::getInstance();
#endif
G4int evno = eventaction->GetEventno() ;//mi dice a che evento siamo 

//prendo l'energia iniziale delle particelle primarie
  initial_energy= p_Primary->GetInitialEnergy();
 
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
  
    //prendo l'energia cinetica della particella  	
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

   G4String particle_name= p_Primary->GetParticle(); 
   
   

if(particle_name=="gamma")
     {    
      if(0 == Step->GetTrack()->GetParentID() ) //primary particle 
	{ 
         G4String process=Step->GetPostStepPoint()->GetProcessDefinedStep()
	   ->GetProcessName();

	 //G4cout<< "processo del fotone-------------------------"
	 //         << process<<G4endl;
          if (process=="LowEnPhotoElec") 
	    { //G4cout<< "processo del fotone-------------------------"
	      // <<" 1"<<G4endl; 
	  runaction->primary_processes(3);
	}
      if (process=="Transportation") 
           { 

runaction->primary_processes(1);	
}
	
 if (process=="LowEnRayleigh") 
           {

	     runaction->primary_processes(2);
	}
 if (process=="LowEnCompton") 
           { 

	     runaction->primary_processes(4);	 
  }

if (process=="LowEnConversion") 
           { 
	     runaction->primary_processes(5);	
	   }


	
      //per il coeff di attenuazione// 
      // G4cout<<Step->GetPreStepPoint()->GetPhysicalVolume()->GetName()<<G4endl; 
 if(Step->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Target"){
   
    if(Step->GetTrack()->GetNextVolume()->GetName() == "World" ) {

   //volume in cui si esce e' il target 
    if(XMoD==0 && YMoD==0 && ZMoD==1){
      
      if(KinE ==initial_energy)
	{


 runaction->Trans_number();
	}}}}}}
 else {if(particle_name=="e-"){
  
   //new particle    
if(IDnow != IDold){ 

   //range
range=0.;
IDold=IDnow ;  
if(0 ==Step->GetTrack()->GetParentID() ) 
  {
  
 
 if(0.0 == KinE) 
                   {  
		   
                     G4double  xend= Step->GetTrack()->GetPosition().x()/mm ;
                     G4double yend= Step->GetTrack()->GetPosition().y()/mm ;
                     G4double  zend= Step->GetTrack()->GetPosition().z()/mm ;
       
		      range=(sqrt(xend*xend+yend*yend+zend*zend)); 
		    
		      G4cout<<"range di e- in cm: "<<range/cm<<G4endl;
		     
 
 
		  
		   }

  }
}
 }
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
   


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......













