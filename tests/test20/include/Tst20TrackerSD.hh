// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20TrackerSD.hh,v 1.1 2001-05-24 19:49:21 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ ComptonTelTrackerSD  ------
// ************************************************************

#ifndef Tst20TrackerSD_h
#define Tst20TrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class Tst20DetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "Tst20TrackerHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst20TrackerSD : public G4VSensitiveDetector
{
public:
  
  Tst20TrackerSD(G4String, Tst20DetectorConstruction* );
  ~Tst20TrackerSD();
  
  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step* astep,G4TouchableHistory* ROHist);
  void EndOfEvent(G4HCofThisEvent*);
  void clear();
  void DrawAll();
  void PrintAll();
  
private:
  
  Tst20TrackerHitsCollection*  TrackerCollection;      
  Tst20DetectorConstruction* Detector;
  G4int (*HitID)[30];
  G4int NbOfTKRLayers;
  G4int NbOfTKRPixels;

};

#endif






