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
#include "hTestDetectorConstruction.hh"
#include "G4Run.hh"
#include "globals.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4IonC12.hh"
#include "G4IonAr40.hh"
#include "G4IonFe56.hh"
#include "G4Electron.hh"
#include "g4std/iostream"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestRunMessenger;
class HepTupleManager;
class HepTuple;
class HepHistogram;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestRunAction : public G4UserRunAction
{
public: // Without description

    hTestRunAction(hTestDetectorConstruction*);
   ~hTestRunAction();

public: // With description
 
    void BeginOfRunAction(const G4Run*);
  // In this method histogramms are booked

    void EndOfRunAction(const G4Run*);
  // In this method bookHisto method is called in which histogramms are filled

public: // Without description

    void SethistName(G4String name) {histName = name;};
    void bookHisto();
    HepTuple* GetNtuple() const {return ntup;};
    void SaveToTuple(G4String, G4double);
    void SaveToTuple(G4String, G4double, G4double);
    void SaveEvent();
    void AddEnergy(G4double, G4double);
    void SetVerbose(G4int val) {verbose = val;};
    G4int GetVerbose() const {return verbose;};
    void SetHistoNumber(G4int val) {nHisto = val;};

  private:

    hTestDetectorConstruction* theDet;
    hTestRunMessenger* runMessenger;

    G4String histName ;
    G4std::vector<HepHistogram*> histo;
    HepTupleManager* hbookManager;
    HepTuple* ntup;
    G4int nHisto;
    G4int verbose; 
};

#endif

