// Rich advanced example for Geant4
// RichTbEventAction.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbEventAction_h
#define RichTbEventAction_h 1
#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "RichTbRunConfig.hh"
#include "RichTbAnalysisManager.hh"

class G4Event;
class RichTbVisManager;
class RichTbAnalysisManager;
class RichTbIOData;

class RichTbEventAction : public G4UserEventAction {

 public:
  RichTbEventAction();

  RichTbEventAction(RichTbRunConfig* ,
		    RichTbVisManager* , RichTbIOData* );
  virtual ~RichTbEventAction();
public:
    void BeginOfEventAction(const G4Event* );
    void EndOfEventAction(const G4Event* );
    G4int GetRichCollID() {return RichTbCollID;};
  RichTbAnalysisManager* getAnalysisM()
  {return ranalysisManager; }
  RichTbVisManager* getVisM()
  {return pVisManager; }
 private:
  RichTbRunConfig* runConfiguration;
  G4int RichTbCollID;
  RichTbAnalysisManager* ranalysisManager;
  RichTbVisManager* pVisManager;
  RichTbIOData* rTbIOData;

};
#endif
