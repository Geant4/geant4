#include "G4PIScorer.hh"
#include "G4Step.hh"
#include "G4PStep.hh"
#include "G4Sigma.hh"
#include "G4StateManager.hh"
#include "G4VIStore.hh"
#include "G4PStepStream.hh"
#include "G4StepPoint.hh"

G4PIScorer::G4PIScorer(const G4VIStore &IStore) : fIStore(IStore),
  fNerr(0),fCorrectWeight(true){}
G4PIScorer::~G4PIScorer(){}

void G4PIScorer::Score(const G4Step &aStep, const G4PStep &aPstep){
  G4Track *track = aStep.GetTrack();


  if (track->GetTrackStatus()==fStopAndKill) {
    G4cout << "      track status is StopAndKill -> do nothing" << G4endl;
  }
  // do the scoring
  else {    
    G4double weight = track->GetWeight();
    G4double import = fIStore.GetImportance(aPstep.fPostTouchableKey);
    G4double i_w =import  * weight;
    if (i_w != 1.) {
      fNerr++;
      G4StepPoint* prepoint = aStep.GetPreStepPoint();
      G4StepPoint* postpoint = aStep.GetPostStepPoint();
      
      G4cout << "G4PIScorer::Score:" << G4endl;
      G4cout << "  Stepnumber = " << aStep.GetTrack()->GetCurrentStepNumber() 
	     << G4endl;
      G4cout << "  TrackID = " << aStep.GetTrack()->GetTrackID() << G4endl;
      G4cout << "  ParentID = " << aStep.GetTrack()->GetParentID() << G4endl;
      G4cout << "  Track wieght = " << weight << G4endl;
      G4cout << "  Importance = " << import << G4endl;
      G4cout << "  i_w = " << i_w << G4endl;
      G4cout << "  aPstep: " << aPstep << G4endl;
      G4cout << "  steplength = " << aStep.GetStepLength() << G4endl;
      G4cout << "  pre_pos = " << prepoint->GetPosition()
	     << ", pre_dir = " << prepoint->GetMomentumDirection() << G4endl;
      G4cout << "  post_pos = " << postpoint->GetPosition()
	     << ", post_dir = " << postpoint->GetMomentumDirection() << G4endl;
      G4cout << "  vxt mom: " << track->GetVertexMomentumDirection() << G4endl;
      fCorrectWeight = false;
    }
    fPScorer.Score(aStep, aPstep);
  }
}

ostream& operator<<(ostream &out, const G4PIScorer &ps) {
  out << ps.GetMapPtkTallys();
  return out;
}










