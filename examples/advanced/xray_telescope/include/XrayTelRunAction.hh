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
// * MODULE:            XrayTelRunAction.hh                             *
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

#ifndef XrayTelRunAction_h
#define XrayTelRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "g4std/vector"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;
class XrayTelAnalysisManager;

class XrayTelRunAction : public G4UserRunAction
{
public:
  XrayTelRunAction(XrayTelAnalysisManager* = 0);
  //		   G4std::vector<G4double*> *enEnergy,
  //	   G4std::vector<G4ThreeVector*> *enDirect,
  //	   G4bool* dEvent);
  ~XrayTelRunAction();

public:
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

private:
  XrayTelAnalysisManager* fAnalysisManager;
  //G4bool* drawEvent;
  //G4std::vector<G4double*>* EnteringEnergy;
  //G4std::vector<G4ThreeVector*>* EnteringDirection;

};

#endif

