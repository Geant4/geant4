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

class G4CaloSD;

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
  G4CaloSD* theCaloSD[maxNumberOfSD];
  G4bool isInitialized;
private:
  void initialize();    
  G4bool trackStartsInCalo(const G4Track* atrack);
  void setPrimaryID(G4int id);
      
};

#endif
