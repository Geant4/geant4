#ifndef DetectorSliceRunAction_h
#define DetectorSliceRunAction_h

#include "G4UserRunAction.hh"

class G4Run;


class DetectorSliceRunAction: public G4UserRunAction {

public:

  DetectorSliceRunAction(){};
  virtual ~DetectorSliceRunAction();
      
  virtual void BeginOfRunAction( const G4Run* aRun );    
  virtual void EndOfRunAction( const G4Run* aRun );    

};

#endif
