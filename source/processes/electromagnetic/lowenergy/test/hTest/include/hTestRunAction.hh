#ifndef hTestRunAction_h
#define hTestRunAction_h 1

// -------------------------------------------------------------
//
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------
//      GEANT4 hTest
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- hTestCalorimeterSD -------------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestRunAction : public G4UserRunAction
{
public: // Without description

    hTestRunAction();
   ~hTestRunAction();

public: // With description
 
    void BeginOfRunAction(const G4Run*);
  // In this method histogramms are booked

    void EndOfRunAction(const G4Run*);
  // In this method bookHisto method is called in which histogramms are filled

};

#endif

