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
//              showers in any calorimeter part.
//              For the moment it is requested only by CCalSensAssign
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalStackingAction_h
#define CCalStackingAction_h 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class CCaloSD;

class CCalStackingAction : public G4UserStackingAction {

  friend class CCalSensAssign;
  
private:
  CCalStackingAction();
public:
  ~CCalStackingAction();

public:
  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
  virtual void NewStage();
  virtual void PrepareNewEvent();
      
public:
  enum stageLevel {firstStage, end} ; 
  enum {maxNumberOfSD = 30};

private:
  stageLevel stage;
  int numberOfSD;
  G4String SDName[maxNumberOfSD];
  int nurgent;
  int acceptSecondaries;
  CCaloSD* theCaloSD[maxNumberOfSD];
  G4bool isInitialized;
private:
  void initialize();    
  G4bool trackStartsInCalo(const G4Track* atrack);
  void setPrimaryID(G4int id);
      
};

#endif
