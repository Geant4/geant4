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
#include <vector>

//

class ThyroidSteppingAction : public G4UserSteppingAction
{
public:
  ThyroidSteppingAction();
  //			std::vector<G4double*> *enEnergy, 
  //		std::vector<G4ThreeVector*> *enDirect,
  //		G4bool* dEvent,
  //		XrayTelAnalysisManager* = 0);
  virtual ~ThyroidSteppingAction();

  virtual void UserSteppingAction(const G4Step*);
  
private:
  G4bool* drawEvent;
  //  std::vector<G4double*>* enteringEnergy;
  // std::vector<G4ThreeVector*>* enteringDirection;

  // ThyroidAnalysisManager* fAnalysisManager;
};

#endif
