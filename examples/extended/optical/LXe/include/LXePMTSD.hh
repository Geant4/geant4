#ifndef LXePMTSD_h
#define LXePMTSD_h 1

#include "G4DataVector.hh"
#include "G4VSensitiveDetector.hh"
#include "LXePMTHit.hh"
class G4Step;
class G4HCofThisEvent;

class LXePMTSD : public G4VSensitiveDetector
{

public:
  LXePMTSD(G4String name);
  ~LXePMTSD();
  
  void Initialize(G4HCofThisEvent* HCE);
  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  
  //A version of processHits that keeps aStep constant
  G4bool ProcessHits_constStep(const G4Step* aStep,
			       G4TouchableHistory* ROhist);
  void EndOfEvent(G4HCofThisEvent* HCE);
  void clear();
  void DrawAll();
  void PrintAll();
  
  //Initialize the arrays to store pmt possitions
  inline void InitPMTs(G4int nPMTs){
    if(pmtPositionsX)delete pmtPositionsX;
    if(pmtPositionsY)delete pmtPositionsY;
    if(pmtPositionsZ)delete pmtPositionsZ;
    pmtPositionsX=new G4DataVector(nPMTs);
    pmtPositionsY=new G4DataVector(nPMTs);
    pmtPositionsZ=new G4DataVector(nPMTs);
  }

  //Store a pmt position
  inline void SetPMTPos(G4int n,G4double x,G4double y,G4double z){
    if(pmtPositionsX)pmtPositionsX->insertAt(n,x);
    if(pmtPositionsY)pmtPositionsY->insertAt(n,y);
    if(pmtPositionsZ)pmtPositionsZ->insertAt(n,z);
  }
  
private:
  LXePMTHitsCollection* pmtHitCollection;
  
  G4DataVector* pmtPositionsX;
  G4DataVector* pmtPositionsY;
  G4DataVector* pmtPositionsZ;
};

#endif

