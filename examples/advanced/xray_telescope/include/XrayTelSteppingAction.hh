// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelSteppingAction.hh                        *
// * -------                                                            *
// *                                                                    *
// * Version:           0.4                                             *
// * Date:              06/11/00                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
//
// CHANGE HISTORY
// --------------
//
// 06.11.2000 R.Nartallo
// - First implementation of X-ray Telescope advanced example.
// - Based on Chandra and XMM models
//
//
// **********************************************************************

#ifndef XrayTelSteppingAction_h
#define XrayTelSteppingAction_h 1

class XrayTelHistogram;

#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
#include "g4std/vector"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayTelSteppingAction : public G4UserSteppingAction
{
public:
  XrayTelSteppingAction(  
			G4std::vector<G4double*> *enEnergy, 
			G4std::vector<G4ThreeVector*> *enDirect,
			G4bool* dEvent);
  virtual ~XrayTelSteppingAction();

  virtual void UserSteppingAction(const G4Step*);
  
private:
  //    XrayTelHistogram* histoManager;
  G4bool* drawEvent;
  G4std::vector<G4double*>* EnteringEnergy;
  G4std::vector<G4ThreeVector*>* EnteringDirection;
};

#endif
