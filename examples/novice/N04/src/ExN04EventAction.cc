
#include "ExN04EventAction.hh"

#include "ExN04TrackerHit.hh"
#include "ExN04CalorimeterHit.hh"
#include "ExN04MuonHit.hh"

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

ExN04EventAction::ExN04EventAction()
{
  trackerCollID = -1;
  calorimeterCollID = -1;
  muonCollID = -1;
}

ExN04EventAction::~ExN04EventAction()
{;}

void ExN04EventAction::BeginOfEventAction()
{
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  if(trackerCollID<0||calorimeterCollID<0||muonCollID<0)
  {
    G4String colNam;
    trackerCollID = SDman->GetCollectionID(colNam="trackerCollection");
    calorimeterCollID = SDman->GetCollectionID(colNam="calCollection");
    muonCollID = SDman->GetCollectionID(colNam="muonCollection");
  }
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/draw/current");
  }
}

void ExN04EventAction::EndOfEventAction()
{
  const G4Event* evt = fpEventManager->GetConstCurrentEvent();

  G4cout << ">>> Event " << evt->GetEventID() << endl;

  if(trackerCollID<0||calorimeterCollID<0||muonCollID<0) return;

  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  ExN04TrackerHitsCollection* THC = NULL;
  ExN04CalorimeterHitsCollection* CHC = NULL;
  ExN04MuonHitsCollection* MHC = NULL;
  if(HCE)
  {
    THC = (ExN04TrackerHitsCollection*)(HCE->GetHC(trackerCollID));
    CHC = (ExN04CalorimeterHitsCollection*)(HCE->GetHC(calorimeterCollID));
    MHC = (ExN04MuonHitsCollection*)(HCE->GetHC(muonCollID));
  }

  if(THC)
  {
    int n_hit = THC->entries();
    G4cout << "     " << n_hit
         << " hits are stored in ExN04TrackerHitsCollection." << endl;
  }
  if(CHC)
  {
    int n_hit = CHC->entries();
    G4cout << "     " << n_hit
         << " hits are stored in ExN04CalorimeterHitsCollection." << endl;
    G4double totE = 0;
    for(int i=0;i<n_hit;i++)
    { totE += (*CHC)[i]->GetEdep(); }
    G4cout << "     Total energy deposition in calorimeter : "
         << totE / GeV << " (GeV)" << endl;
  }
  if(MHC)
  {
    int n_hit = MHC->entries();
    G4cout << "     " << n_hit
         << " hits are stored in ExN04MuonHitsCollection." << endl;
  }

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    if(THC) THC->DrawAllHits();
    if(CHC) CHC->DrawAllHits();
    if(MHC) MHC->DrawAllHits();
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/show/view");
  }
}



