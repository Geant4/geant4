// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T07EventAction.cc,v 1.2 1999-04-17 07:00:41 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "T07EventAction.hh"

#include "T07CalorHit.hh"
#include "T07EventActionMessenger.hh"

#include <rw/tvordvec.h>

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


T07EventAction::T07EventAction()
:calorimeterCollID(-1),drawFlag("all"),eventMessenger(NULL)
#ifndef MAKEHBOOK
,EnergyAbsorber(250,0.,750.),EnergyGap(100,0.,100.),
TrackLengthAbsorber(100,0.,500.),TrackLengthGap(100,0.,500.)
#endif
{
  eventMessenger = new T07EventActionMessenger(this);
#ifdef MAKEHBOOK
    EnergyAbsorber=new HbookHistogram("Energy Absorber",100,0.,500.);
    EnergyGap=new HbookHistogram("Energy Gap",100,0.,50.);
    TrackLengthAbsorber=new HbookHistogram("Track Length Absorber",
					   100,0.,300.);
    TrackLengthGap=new HbookHistogram("Track Gap Absorber",
				      100,0.,200.);
#endif    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

T07EventAction::~T07EventAction()
{
  delete eventMessenger;
#ifndef MAKEHBOOK
  ofstream o("test07.plt");
  o << "# Histos 4" << endl;
  EnergyAbsorber.output(o);
  o << endl;
  o << endl;
  EnergyGap.output(o);
  o << endl;
  o << endl;
  TrackLengthAbsorber.output(o);
  o << endl;
  o << endl;
  TrackLengthGap.output(o);
  o.close();
#endif
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void T07EventAction::BeginOfEventAction(const G4Event*)
{  if(calorimeterCollID==-1)
  {
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    calorimeterCollID = SDman->GetCollectionID("CalCollection");
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void T07EventAction::EndOfEventAction( const G4Event* evt )
{
 
#ifdef MAKEHBOOK
  if(evt->GetEventID()%100==0)
    G4cout << ">>> Event " << evt->GetEventID() << endl;
#else  
  // G4cout << ">>> Event " << evt->GetEventID() << endl;
#endif

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  T07CalorHitsCollection* CHC = NULL;
  if (HCE)
      CHC = (T07CalorHitsCollection*)(HCE->GetHC(calorimeterCollID));

  if (CHC)
   {
    int n_hit = CHC->entries();
#ifndef MAKEHBOOK
    //G4cout << "     " << n_hit
    //     << " hits are stored in T07CalorHitsCollection." << endl;
#endif
    G4double totEAbs=0, totLAbs=0, totEGap=0, totLGap=0;
    for (int i=0;i<n_hit;i++)
      { totEAbs += (*CHC)[i]->GetEdepAbs(); 
        totLAbs += (*CHC)[i]->GetTrakAbs();
        totEGap += (*CHC)[i]->GetEdepGap(); 
        totLGap += (*CHC)[i]->GetTrakGap();
        
      }
#ifdef MAKEHBOOK
    EnergyAbsorber->accumulate(totEAbs*MeV);
    EnergyGap->accumulate(totEGap*MeV);
    TrackLengthAbsorber->accumulate(totLAbs*mm);
    TrackLengthGap->accumulate(totLGap*mm);
#else    
    EnergyAbsorber.accumulate(totEAbs*MeV);
    EnergyGap.accumulate(totEGap*MeV);
    TrackLengthAbsorber.accumulate(totLAbs*mm);
    TrackLengthGap.accumulate(totLGap*mm);
    /*
    G4cout
       << "   Absorber: total energy: " << setw(7) << G4BestUnit(totEAbs,"Energy")
       << "       total track length: " << setw(7) << G4BestUnit(totLAbs,"Length")
       << endl
       << "        Gap: total energy: " << setw(7) << G4BestUnit(totEGap,"Energy")
       << "       total track length: " << setw(7) << G4BestUnit(totLGap,"Length")
       << endl;
    */
#endif           
    }

#ifndef MAKEHBOOK
  /*
  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if(trajectoryContainer)
  { n_trajectories = trajectoryContainer->entries(); }
  G4cout << "    " << n_trajectories 
       << " trajectories stored in this event." << endl;


  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
    for(G4int i=0; i<n_trajectories; i++) 
         { G4Trajectory* trj = (*(evt->GetTrajectoryContainer()))[i];
           if (drawFlag == "all") trj->DrawTrajectory(50);
           else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
                                  trj->DrawTrajectory(50); 
         }
  }
  */
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


