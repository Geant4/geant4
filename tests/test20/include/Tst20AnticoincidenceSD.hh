// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20AnticoincidenceSD.hh,v 1.1 2001-05-24 19:49:21 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ Tst20AnticoincidenceSD  ------
// ************************************************************

#ifndef Tst20AnticoincidenceSD_h
#define Tst20AnticoincidenceSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class Tst20DetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "Tst20AnticoincidenceHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst20AnticoincidenceSD : public G4VSensitiveDetector
{
public:
  
  Tst20AnticoincidenceSD(G4String, Tst20DetectorConstruction* );
  ~Tst20AnticoincidenceSD();
  
  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step* astep,G4TouchableHistory* ROHist);
  void EndOfEvent(G4HCofThisEvent*);
  void clear();
  void DrawAll();
  void PrintAll();
  
private:
  
  Tst20AnticoincidenceHitsCollection*  AnticoincidenceCollection;      
  Tst20DetectorConstruction* Detector;
  G4int *HitLateralID;
  G4int *HitTopID;
  G4int NbOfACDLateralTiles;
  G4int NbOfACDTopTiles; 
};

#endif






