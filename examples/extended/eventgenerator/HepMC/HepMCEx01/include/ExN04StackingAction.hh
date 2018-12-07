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
/// \file eventgenerator/HepMC/HepMCEx01/include/ExN04StackingAction.hh
/// \brief Definition of the ExN04StackingAction class
//
//

#ifndef ExN04StackingAction_H
#define ExN04StackingAction_H 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4UserStackingAction.hh"
#include "ExN04TrackerHit.hh"
#include "ExN04MuonHit.hh"

class G4Track;
class ExN04StackingActionMessenger;

class ExN04StackingAction : public G4UserStackingAction {
public:
  ExN04StackingAction();
  virtual ~ExN04StackingAction();

  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
  virtual void NewStage();
  virtual void PrepareNewEvent();

  inline void SetNRequestMuon(G4int val) { fReqMuon = val; }
  inline G4int GetNRequestMuon() const { return fReqMuon; }
  inline void SetNRequestIsoMuon(G4int val) { fReqIsoMuon = val; }
  inline G4int GetNRequestIsoMuon() const { return fReqIsoMuon; }
  inline void SetNIsolation(G4int val) { fReqIso = val; }
  inline G4int GetNIsolation() const { return fReqIso; }
  inline void SetRoIAngle(G4double val) { fAngRoI = val; }
  inline G4double GetRoIAngle() const { return fAngRoI; }

private:
  G4bool InsideRoI(const G4Track * aTrack,G4double ang);
  G4VHitsCollection* GetCollection(G4String colName);

  ExN04TrackerHitsCollection* fTrkHits;
  ExN04MuonHitsCollection* fMuonHits;
  ExN04StackingActionMessenger* fMessenger;

  G4int fStage;
  G4int fReqMuon;
  G4int fReqIsoMuon;
  G4int fReqIso;
  G4double fAngRoI;
};

#endif
