///////////////////////////////////////////////////////////////////////////////
// File: CMSStackingAction.h
// Author: Veronique Lefebure
// Last modified: 
// 07/10/99 V.L.
//              -->Action needed to keep track of the ID of
//                 particle generating showers
//		   in any calorimeter part.
//                 For the moment it is requested only by G4CaloSD
//		   ---> private constructor.
// Modifications:
///////////////////////////////////////////////////////////////////////////////
// Defines actions to be taken in each step
//

#ifndef CMSStackingAction_h
#define CMSStackingAction_h 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class G4CaloSD;

class CMSStackingAction : public G4UserStackingAction {

  friend class CMSSensAssign;
  
private:
  CMSStackingAction();
public:
  ~CMSStackingAction();

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
