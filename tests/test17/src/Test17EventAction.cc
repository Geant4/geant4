// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Test17EventAction.hh"

#include "Test17RunAction.hh"

#include "Test17CalorHit.hh"
#include "Test17EventActionMessenger.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17EventAction::Test17EventAction(Test17RunAction* Test17RA)
:calorimeterCollID(-1),eventMessenger(NULL),
 verboselevel(0),runaction(Test17RA),drawFlag("all")
{
  eventMessenger = new Test17EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17EventAction::~Test17EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17EventAction::BeginOfEventAction(const G4Event* evt)
{  


  if(calorimeterCollID==-1)
  {
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    calorimeterCollID = SDman->GetCollectionID("CalCollection");
  } 

  EnergyDeposition = 0.0;

  if(verboselevel>1)
    G4cout << "<<< Event  " << evt->GetEventID() << " started." << G4endl;
  nstep = 0. ;
  nstepCharged = 0. ;
  nstepNeutral = 0. ;
  Nch = 0. ;
  Nne = 0. ;
  NE=0.;
  NP=0.;
  Transmitted=0.;
  Reflected  =0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17EventAction::EndOfEventAction(const G4Event* evt)
{

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  Test17CalorHitsCollection* CHC = NULL;
  if (HCE)
      CHC = (Test17CalorHitsCollection*)(HCE->GetHC(calorimeterCollID));

  if (CHC)
   {
    int n_hit = CHC->entries();
   // if(verboselevel==2)
   // G4cout << "     " << n_hit
   //      << " hits are stored in Test17CalorHitsCollection." << G4endl;

    G4double totEAbs=0, totLAbs=0;
    for (int i=0;i<n_hit;i++)
      { totEAbs += (*CHC)[i]->GetEdepAbs(); 
        totLAbs += (*CHC)[i]->GetTrakAbs();
      }
  if(verboselevel==2)
    G4cout
       << "   Absorber: total energy: " << G4std::setw(7) << 
                             G4BestUnit(totEAbs,"Energy")
       << "       total track length: " << G4std::setw(7) <<
                             G4BestUnit(totLAbs,"Length")
       << G4endl;           

   // count event, add deposits to the sum ...
    runaction->CountEvent() ;
    runaction->AddTrackLength(totLAbs) ;
    runaction->AddnStepsNeutral(nstepNeutral) ;
    if(verboselevel==2)
      G4cout << " Ncharged=" << Nch << "  ,   Nneutral=" << Nne << G4endl;
    runaction->CountParticles(Nch,Nne);
    runaction->AddEP(NE,NP);
    runaction->AddTrRef(Transmitted,Reflected) ;
    runaction->AddEdeps(totEAbs) ;
    runaction->FillEn(totEAbs) ;

    nstep=nstepCharged+nstepNeutral ;
    runaction->FillNbOfSteps(nstep);
    runaction->SaveToTuple("EDEP",EnergyDeposition);      

    runaction->SaveEvent();
  }


  if(verboselevel>0)
    G4cout << "<<< Event  " << evt->GetEventID() << " ended." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int Test17EventAction::GetEventno()
{
  G4int evno = fpEventManager->GetConstCurrentEvent()->GetEventID() ;
  return evno ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17EventAction::setEventVerbose(G4int level)
{
  verboselevel = level ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17EventAction::CountStepsCharged()
{
  nstepCharged += 1. ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17EventAction::CountStepsNeutral()
{
  nstepNeutral += 1. ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17EventAction::AddCharged() 
{
  Nch += 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17EventAction::AddNeutral() 
{
  Nne += 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17EventAction::AddE(G4double En) 
{
  EnergyDeposition += En;
  NE += 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17EventAction::AddP() 
{
  NP += 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17EventAction::SetTr()
{
  Transmitted = 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17EventAction::SetRef()
{
  Reflected   = 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  









