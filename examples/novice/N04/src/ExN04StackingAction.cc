
#include "ExN04StackingAction.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "ExN04StackingActionMessenger.hh"
#include "G4ios.hh"

ExN04StackingAction::ExN04StackingAction()
:stage(0),trkHits(NULL),muonHits(NULL)
{ 
  angRoI = 30.0*deg; 
  reqMuon = 2;
  reqIso = 10;
  theMessenger = new ExN04StackingActionMessenger(this);
}

ExN04StackingAction::~ExN04StackingAction()
{ delete theMessenger; }

G4ClassificationOfNewTrack 
ExN04StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  G4ClassificationOfNewTrack classification = fWaiting;
  switch(stage)
  {
  case 0: // Stage 0 : Primary muons only
    if(aTrack->GetParentID()==0)
    {
      G4ParticleDefinition * particleType = aTrack->GetDefinition();
      if((particleType==G4MuonPlus::MuonPlusDefinition())
       ||(particleType==G4MuonMinus::MuonMinusDefinition()))
      { classification = fUrgent; }
    }
    break;

  case 1: // Stage 1 : Charged primaries only
          //           Suspended tracks will be sent to the waiting stack
    if(aTrack->GetParentID()!=0) { break; }
    if(aTrack->GetTrackStatus()==fSuspend) { break; }
    if(aTrack->GetDefinition()->GetPDGCharge()==0.) { break; }
    classification = fUrgent;
    break;

  default: // Stage 2 : Accept all primaries
           //           Accept all secondaries in RoI
           //           Kill secondaries outside RoI
    if(aTrack->GetParentID()==0)
    { 
      classification = fUrgent;
      break;
    }
    if(InsideRoI(aTrack,angRoI))
    { 
      classification = fUrgent;
      break;
    }
    classification = fKill;
  }
  return classification;
}

G4bool ExN04StackingAction::InsideRoI(G4Track *const aTrack,G4double ang)
{
  if(!muonHits)
  { muonHits = (ExN04MuonHitsCollection*)GetCollection("muonCollection"); }
  if(!muonHits)
  { G4cerr << "muonCollection NOT FOUND" << endl;
    return true; }

  G4int nhits = muonHits->entries();

  const G4ThreeVector trPos = aTrack->GetPosition();
  for(G4int i=0;i<nhits;i++)
  {
    G4ThreeVector muHitPos = (*muonHits)[i]->GetPos();
    G4double angl = muHitPos.angle(trPos);
    if(angl<ang) { return true; }
  }

  return false;
}

G4VHitsCollection* ExN04StackingAction::GetCollection(G4String colName)
{
  G4SDManager* SDMan = G4SDManager::GetSDMpointer();
  G4RunManager* runMan = G4RunManager::GetRunManager();
  int colID = SDMan->GetCollectionID(colName);
  if(colID>=0)
  {
    const G4Event* currentEvent = runMan->GetCurrentEvent();
    G4HCofThisEvent* HCE = currentEvent->GetHCofThisEvent();
    return HCE->GetHC(colID);
  }
  return NULL;
}

void ExN04StackingAction::NewStage()
{
  stage++;
  G4int nhits;
  if(stage==1)
  {
  // Stage 0->1 : check if at least "reqMuon" hits on muon chamber
  //              otherwise abort current event
    if(!muonHits)
    { muonHits = (ExN04MuonHitsCollection*)GetCollection("muonCollection"); }
    if(!muonHits)
    { G4cerr << "muonCollection NOT FOUND" << endl;
      return; }
    nhits = muonHits->entries();
    G4cout << "Stage 0->1 : " << nhits << " hits found in the muon chamber."
         << endl;
    if(nhits<reqMuon)
    { 
      stackManager->clear();
      G4cout << "++++++++ event aborted" << endl;
      return;
    }
    stackManager->ReClassify();
    return;
  }

  else if(stage==2)
  {
  // Stage 1->2 : check the isolation of muon tracks
  //              at least "reqIsoMuon" isolated muons
  //              otherwise abort current event.
  //              Isolation requires "reqIso" or less hits
  //              (including own hits) in the RoI region
  //              in the tracker layers.
    nhits = muonHits->entries();
    if(!trkHits)
    { trkHits = (ExN04TrackerHitsCollection*)GetCollection("trackerCollection"); }
    if(!trkHits)
    { G4cerr << "trackerCollection NOT FOUND" << endl;
      return; }
    G4int nTrkhits = trkHits->entries();
    G4int isoMuon = 0;
    for(G4int j=0;j<nhits;j++)
    {
      G4ThreeVector hitPos = (*muonHits)[j]->GetPos();
      G4int nhitIn = 0;
      for(G4int jj=0;(jj<nTrkhits)&&(nhitIn<=reqIso);jj++)
      {
        G4ThreeVector trkhitPos = (*trkHits)[jj]->GetPos();
        if(trkhitPos.angle(hitPos)<angRoI) nhitIn++;
      }
      if(nhitIn<=reqIso) isoMuon++;
    }
    G4cout << "Stage 1->2 : " << isoMuon << " isolated muon found." << endl;
    if(isoMuon<reqIsoMuon)
    {
      stackManager->clear();
      G4cout << "++++++++ event aborted" << endl;
      return;
    }
    stackManager->ReClassify();
    return;
  }

  else
  {
  // Other stage change : just re-classify
    stackManager->ReClassify();
  }
}
    
void ExN04StackingAction::PrepareNewEvent()
{ 
  stage = 0; 
  trkHits = NULL;
  muonHits = NULL;
}


