//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

#ifndef ExN04StackingAction_H
#define ExN04StackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "G4ThreeVector.hh"

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

