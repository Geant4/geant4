//
//
// **********************************************************************
// *                                                                    *         
// *                    ThyriodSteppingAction.hh                        *
// *                                                                    *
// **********************************************************************
//


#ifndef ThyroidSteppingAction_h
#define ThyroidSteppingAction_h 1

// class ThyroidAnalysisManager;

#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
#include "g4std/vector"

//

class ThyroidSteppingAction : public G4UserSteppingAction
{
public:
  ThyroidSteppingAction();
  //			G4std::vector<G4double*> *enEnergy, 
  //		G4std::vector<G4ThreeVector*> *enDirect,
  //		G4bool* dEvent,
  //		XrayTelAnalysisManager* = 0);
  virtual ~ThyroidSteppingAction();

  virtual void UserSteppingAction(const G4Step*);
  
private:
  G4bool* drawEvent;
  //  G4std::vector<G4double*>* enteringEnergy;
  // G4std::vector<G4ThreeVector*>* enteringDirection;

  // ThyroidAnalysisManager* fAnalysisManager;
};

#endif
