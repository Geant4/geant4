#include "B01Scorer.hh"
#include "G4Step.hh"
#include "G4PStep.hh"
#include "G4Sigma.hh"
#include "G4StateManager.hh"


B01Scorer::B01Scorer(){}
B01Scorer::~B01Scorer(){}

void B01Scorer::Score(const G4Step &aStep, const G4PStep &aPstep){
  G4Track *track = aStep.GetTrack();

  if (track->GetTrackStatus()==fStopAndKill) {
    G4cout << "      track status is StopAndKill -> do nothing" << G4endl;
  }
  // do the scoring
  else {    
    // the map fPtkTallys, the map nametallys and the "G4Sigma" will
    // be setup the first time they are accesed. 
    // the user may choos any other place to store the results e.g. 
    // histogramms
    G4PTouchableKey post_ptk(aPstep.fPostTouchableKey); 
    if (aPstep.fCrossBoundary) { 
      // Pstep crosses boundary
      fPtkTallys[post_ptk]["HistorysEntering"].Xin(1);
      fPtkTallys[post_ptk]["EnergyEnteringHistory"].
	Xin(track->GetKineticEnergy());
    } 
    else { 
      fPtkTallys[post_ptk]["Collisions"].Xin(1);
    }
  }
}

G4std::ostream& operator<<(G4std::ostream &out, const B01Scorer &ps) {
  out << ps.GetMapPtkTallys();
  return out;
}










