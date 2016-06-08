
#ifndef ExN04StackingAction_H
#define ExN04StackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "G4ThreeVector.hh"
#include <rw/tpordvec.h>

class G4Track;

#include "ExN04TrackerHit.hh"
#include "ExN04MuonHit.hh"
class ExN04StackingActionMessenger;

class ExN04StackingAction : public G4UserStackingAction
{
  public:
    ExN04StackingAction();
    virtual ~ExN04StackingAction();

  public:
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    virtual void NewStage();
    virtual void PrepareNewEvent();

  private:
    G4bool InsideRoI(const G4Track * aTrack,G4double ang);
    G4VHitsCollection* GetCollection(G4String colName);
    
    ExN04TrackerHitsCollection* trkHits;
    ExN04MuonHitsCollection* muonHits;
    ExN04StackingActionMessenger* theMessenger;

    G4int stage;
    G4int reqMuon;
    G4int reqIsoMuon;
    G4int reqIso;
    G4double angRoI;
  
  public:
    inline void SetNRequestMuon(G4int val) { reqMuon = val; }
    inline G4int GetNRequestMuon() const { return reqMuon; }
    inline void SetNRequestIsoMuon(G4int val) { reqIsoMuon = val; }
    inline G4int GetNRequestIsoMuon() const { return reqIsoMuon; }
    inline void SetNIsolation(G4int val) { reqIso = val; }
    inline G4int GetNIsolation() const { return reqIso; }
    inline void SetRoIAngle(G4double val) { angRoI = val; }
    inline G4double GetRoIAngle() const { return angRoI; }
};

#endif

