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
// $Id: Tst50SteppingAction.cc,v 1.7 2003-01-08 15:37:14 guatelli Exp $
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

    G4double XIMD;
    G4double YIMD;
    G4double ZIMD;


    G4double MMoD;
    G4double MIMD;

    G4double Theta;

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

// albe 8/1/2003 rifaccio distribuzione di theta

    // XIMD = p_primary->GetInitialMomentumDirection().x();
    //YIMD = p_primary->GetInitialMomentumDirection().y();
    //ZIMD = p_primary->GetInitialMomentumDirection().z();

    /*	
    MMoD=sqrt(pow(XMoD,2.0)+pow(YMoD,2.0)+pow(ZMoD,2.0));
    MIMD=sqrt(pow(XIMD,2.0)+pow(YIMD,2.0)+pow(ZIMD,2.0));

    if (MMoD == 0) {
	Theta = 0;
    } else if (MIMD == 0) {
	Theta = 0;
    } else {   
	Theta = acos((XMoD*XIMD+YMoD*YIMD+ZMoD*ZIMD)/(MMoD*MIMD))
    }
    */
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
 else {

if(particle_name=="e-"){

    //new particle    
if(IDnow != IDold){ 

 IDold=IDnow ; 



if(0 ==Step->GetTrack()->GetParentID() )  
{ 
 G4double  XPos_pre=Step->GetPreStepPoint()->GetPosition().x();
 G4double YPos_pre=Step->GetPreStepPoint()->GetPosition().y();
 G4double ZPos_pre=Step->GetPreStepPoint()->GetPosition().z();


 if (0== XPos_pre && 0== YPos_pre && 0== ZPos_pre)

    { 
      G4double energyLost=abs(Step->GetDeltaEnergy());
      G4double stepLength= Step->GetTrack()-> GetStepLength();
  
      G4cout<<"lunghezza dello step in cm:"<<stepLength/cm<<G4endl;
   
if(stepLength!=0) 
{
G4double TotalStoppingPower=(energyLost/stepLength);
 G4cout<<"TotalStoppingPower in MeV/cm:"<<TotalStoppingPower/(MeV/cm)<<G4endl;}}
}}

range=0.;

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
 //Radiation yield
 G4String process=Step->GetPostStepPoint()->GetProcessDefinedStep()
	   ->GetProcessName();


if (process=="LowEnBrem") 
  { G4cout <<Step->GetDeltaEnergy()/MeV<<G4endl;
 G4double energyLostforBremm=abs(Step->GetDeltaEnergy());
 eventaction->RadiationYield( energyLostforBremm);

  }

}}}}






















