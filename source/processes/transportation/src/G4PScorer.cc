#include "G4PScorer.hh"
#include "G4Step.hh"
#include "G4PStep.hh"
#include "G4Sigma.hh"
#include "G4StateManager.hh"

G4PScorer::G4PScorer(){}
G4PScorer::~G4PScorer(){}

void G4PScorer::Score(const G4Step &aStep, const G4PStep &aPstep){
  G4Track *track = aStep.GetTrack();

  if (track->GetTrackStatus()==fStopAndKill) {
    G4cout << "      track status is StopAndKill -> do nothing" << G4endl;
  }
  // do the scoring
  else {    
    G4double weight = track->GetWeight();
    // the map fPtkTallys, the map nametallys and the "G4Sigma" will
    // be setup the first time they are accesed. 
    G4PTouchableKey post_ptk(aPstep.fPostTouchableKey); 
    if (aPstep.fCrossBoundary) { 
      // Pstep crosses boundary
      fPtkTallys[post_ptk]["HistorysEntering"].Xin(1);
      fPtkTallys[post_ptk]["WeighteOfHistorysEntering"].Xin(weight);
      fPtkTallys[post_ptk]["EnergyEnteringHistory"].
	Xin(track->GetKineticEnergy());
      fPtkTallys[post_ptk]["WeightedEnergyEnteringHistory"].
	Xin(track->GetKineticEnergy(), weight);
    } 
    else { 
      // Pstep with both points in the same I volume
      if (aStep.GetPostStepPoint()->GetStepStatus() != fGeomBoundary) {
	fPtkTallys[post_ptk]["Collisions"].Xin(1);
	fPtkTallys[post_ptk]["WeighteOfCollisions"].Xin(weight);
      }
    }
  }
}

ostream& operator<<(ostream &out, const G4PScorer &ps) {
  out << ps.GetMapPtkTallys();
  return out;
}










