//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: RE05StackingAction.hh 66526 2012-12-19 13:41:33Z ihrivnac $
//
/// \file RE05/include/RE05StackingAction.hh
/// \brief Definition of the RE05StackingAction class
//

#ifndef RE05StackingAction_H
#define RE05StackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "G4ThreeVector.hh"

class G4Track;

#include "RE05TrackerHit.hh"
#include "RE05MuonHit.hh"
class RE05StackingActionMessenger;

class RE05StackingAction : public G4UserStackingAction
{
  public:
    RE05StackingAction();
    virtual ~RE05StackingAction();

  public:
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    virtual void NewStage();
    virtual void PrepareNewEvent();

  private:
    G4bool InsideRoI(const G4Track * aTrack,G4double ang);
    G4VHitsCollection* GetCollection(G4String colName);
    
    RE05TrackerHitsCollection* trkHits;
    RE05MuonHitsCollection* muonHits;
    RE05StackingActionMessenger* theMessenger;

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

