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
///////////////////////////////////////////////////////////////////////////////
// File: CCalStackingAction.hh
// Description: Action needed to keep track of the ID of particle generating 
//              showers in any calorimeter part.
//              For the moment it is requested only by CCalSensAssign
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalStackingAction_h
#define CCalStackingAction_h 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class CCaloSD;

class CCalStackingAction : public G4UserStackingAction
{
  friend class CCalSensAssign;
  
public:

  CCalStackingAction();
  ~CCalStackingAction();

  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
  virtual void NewStage();
  virtual void PrepareNewEvent();
      
  enum stageLevel {firstStage, end} ; 
  enum {maxNumberOfSD = 30};

private:
  void initialize();    
  G4bool trackStartsInCalo(const G4Track* atrack);
  void setPrimaryID(G4int id);
      
  G4double fTimeLimit;
  stageLevel stage;
  G4int numberOfSD;
  G4String SDName[maxNumberOfSD];
  G4int nurgent;
  G4int acceptSecondaries;
  CCaloSD* theCaloSD[maxNumberOfSD];
  G4bool isInitialized;
};

#endif
