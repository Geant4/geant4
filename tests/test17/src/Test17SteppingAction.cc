// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Test17DetectorConstruction.hh"
#include "G4EnergyLossTables.hh"
#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"
#include "Test17SteppingAction.hh"
#include "Test17PrimaryGeneratorAction.hh"
#include "Test17EventAction.hh"
#include "Test17RunAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "Test17SteppingMessenger.hh"
#include "G4ios.hh"
#include "g4std/iomanip"
#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17SteppingAction::Test17SteppingAction(Test17DetectorConstruction* DET,
                                     Test17EventAction* EA,
                                     Test17RunAction* RA)
:detector (DET),eventaction (EA),runaction (RA),steppingMessenger(NULL),
 IDold(-1) ,evnoold(-1)
{
  steppingMessenger = new Test17SteppingMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17SteppingAction::~Test17SteppingAction()
{
  delete steppingMessenger ;
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 

  G4double Edep,Theta,Thetaback,Ttrans,Tback,Tsec,Egamma,xend,yend,zend,rend ;
  G4double Tkin ;
  G4int evno = eventaction->GetEventno() ; 

  IDnow = evno+10000*(aStep->GetTrack()->GetTrackID())+
          100000000*(aStep->GetTrack()->GetParentID()); 

  Tkin = aStep->GetTrack()->GetKineticEnergy() ; 
    //   Edep = aStep->GetTotalEnergyDeposit() ;         
  Edep = aStep->GetDeltaEnergy() ;         
  Tsec = Tkin - Edep ;
  Theta = acos(aStep->GetTrack()->GetMomentumDirection().x()) ;

  eventaction->AddE(abs(Edep)) ;

  // new particle
  if(IDnow != IDold) {
    IDold=IDnow ;

    // primary
    if(0 == aStep->GetTrack()->GetParentID() ) {
      runaction->SaveToTuple("TKIN",Tsec/MeV);      
      runaction->SaveToTuple("MASS",(aStep->GetTrack()->
                 GetDynamicParticle()->GetDefinition()->GetPDGMass())/MeV);      
      runaction->SaveToTuple("CHAR",(aStep->GetTrack()->
                 GetDynamicParticle()->GetDefinition()->GetPDGCharge()));      

    // secondary
    } else {
      if(((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->
            GetParticleName()) == "e-") 
                 ||
         ((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->
            GetParticleName()) == "e+") ) {

            runaction->Fillvertexz(aStep->GetTrack()->GetVertexPosition().x());
            eventaction->AddCharged() ;
            runaction->FillTsec(Tsec) ;

      } else {
            eventaction->AddNeutral() ;
      }
    }
  }

    // Primary 
  if(0 == aStep->GetTrack()->GetParentID() ) {
    eventaction->CountStepsCharged(aStep->GetStepLength()) ;

    // Stopping or end of track
    if((0.0 == Tkin) || 
         (aStep->GetTrack()->GetNextVolume()->GetName()=="World") ) {

            xend= aStep->GetTrack()->GetPosition().x()/mm ;
            yend= aStep->GetTrack()->GetPosition().y()/mm ;
            zend= aStep->GetTrack()->GetPosition().z()/mm ;
            runaction->SaveToTuple("XEND",xend,1000.0);      
            runaction->SaveToTuple("YEND",yend,1000.0);      
            runaction->SaveToTuple("ZEND",zend,1000.0);      
            runaction->SaveToTuple("TEND",Tkin/MeV,1000.0);
            runaction->SaveToTuple("TET",Theta/deg,1000.0);      
            runaction->AddnStepsCharged(xend) ;
    }
    rend = sqrt(xend*xend + yend*yend) ;

  //secondary
  } else {

    if ( aStep->GetTrack()->GetNextVolume()->GetName() == "World" ) {

      // charged secondaries forward
      if(aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->
                  GetPDGCharge() != 0.0 ) { 

        eventaction->SetTr();

        if(0.5*pi > Theta) {
            runaction->FillTh(Theta) ;
            runaction->FillTt(Tkin) ;
            yend= aStep->GetTrack()->GetPosition().y() ;
            zend= aStep->GetTrack()->GetPosition().z() ;
            rend = sqrt(yend*yend+zend*zend) ;
            runaction->FillR(rend);

  	 // charged secondaries backword
        } else {
            eventaction->SetRef();
            Thetaback = pi - Theta ;
            runaction->FillThBack(Thetaback) ;
            runaction->FillTb(Tsec) ;
        }
 
        // gammas
      } else {
          if (0.5*pi > Theta) runaction->FillGammaSpectrum(Tkin) ;
      }
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....





