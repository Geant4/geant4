// Rich advanced example for Geant4
// RichTbStackingAction.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbStackingAction_h
#define RichTbStackingAction_h 1

#include "G4UserStackingAction.hh"
class RichTbStackingAction : public G4UserStackingAction {

 public:
    RichTbStackingAction();
  virtual ~RichTbStackingAction();
  void NewStage(){ ;}
  void PrepareNewEvent(){;}


private:
};
#endif
