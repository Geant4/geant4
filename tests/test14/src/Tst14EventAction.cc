// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst14EventAction.cc,v 1.1 1999-05-29 14:12:10 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Tst14EventAction.hh"

#include "Tst14RunAction.hh"

#include "Tst14CalorHit.hh"
#include "Tst14EventActionMessenger.hh"

#include <rw/tvordvec.h>

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst14EventAction::Tst14EventAction(Tst14RunAction* Tst14RA)
:calorimeterCollID(-1),eventMessenger(NULL),
 verboselevel(0),runaction (Tst14RA)
{
  eventMessenger = new Tst14EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst14EventAction::~Tst14EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst14EventAction::BeginOfEventAction(const G4Event* evt)
{  if(calorimeterCollID==-1)
  {
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    calorimeterCollID = SDman->GetCollectionID("CalCollection");
  } 

  if(verboselevel>1)
    G4cout << "<<< Event  " << evt->GetEventID() << " started." << endl;
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

void Tst14EventAction::EndOfEventAction(const G4Event* evt)
{

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  Tst14CalorHitsCollection* CHC = NULL;
  if (HCE)
      CHC = (Tst14CalorHitsCollection*)(HCE->GetHC(calorimeterCollID));

  if (CHC)
   {
    int n_hit = CHC->entries();
   // if(verboselevel==2)
   // G4cout << "     " << n_hit
   //      << " hits are stored in Tst14CalorHitsCollection." << endl;

    G4double totEAbs=0, totLAbs=0;
    for (int i=0;i<n_hit;i++)
      { totEAbs += (*CHC)[i]->GetEdepAbs(); 
        totLAbs += (*CHC)[i]->GetTrakAbs();
      }
  if(verboselevel==2)
    G4cout
       << "   Absorber: total energy: " << setw(7) << 
                             G4BestUnit(totEAbs,"Energy")
       << "       total track length: " << setw(7) <<
                             G4BestUnit(totLAbs,"Length")
       << endl;           

   // count event, add deposits to the sum ...
    runaction->CountEvent() ;
    runaction->AddTrackLength(totLAbs) ;
    runaction->AddnStepsCharged(nstepCharged) ;
    runaction->AddnStepsNeutral(nstepNeutral) ;
    if(verboselevel==2)
      G4cout << " Ncharged=" << Nch << "  ,   Nneutral=" << Nne << endl;
    runaction->CountParticles(Nch,Nne);
    runaction->AddEP(NE,NP);
    runaction->AddTrRef(Transmitted,Reflected) ;
    runaction->AddEdeps(totEAbs) ;
    runaction->FillEn(totEAbs) ;

    nstep=nstepCharged+nstepNeutral ;
    runaction->FillNbOfSteps(nstep);
  }

  if(verboselevel>0)
    G4cout << "<<< Event  " << evt->GetEventID() << " ended." << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int Tst14EventAction::GetEventno()
{
  G4int evno = fpEventManager->GetConstCurrentEvent()->GetEventID() ;
  return evno ;
}

void Tst14EventAction::setEventVerbose(G4int level)
{
  verboselevel = level ;
}
void Tst14EventAction::CountStepsCharged()
{
  nstepCharged += 1. ;
}
void Tst14EventAction::CountStepsNeutral()
{
  nstepNeutral += 1. ;
}
void Tst14EventAction::AddCharged() 
{
  Nch += 1.;
}
void Tst14EventAction::AddNeutral() 
{
  Nne += 1.;
}
void Tst14EventAction::AddE() 
{
  NE += 1.;
}
void Tst14EventAction::AddP() 
{
  NP += 1.;
}
void Tst14EventAction::SetTr()
{
  Transmitted = 1.;
}
void Tst14EventAction::SetRef()
{
  Reflected   = 1.;
}
  

