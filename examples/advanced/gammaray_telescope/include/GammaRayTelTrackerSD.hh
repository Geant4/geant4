// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelTrackerSD.hh,v 1.1 2001-03-05 13:58:21 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelTrackerSD  ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************

#ifndef GammaRayTelTrackerSD_h
#define GammaRayTelTrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class GammaRayTelDetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "GammaRayTelTrackerHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelTrackerSD : public G4VSensitiveDetector
{
public:
  
  GammaRayTelTrackerSD(G4String, GammaRayTelDetectorConstruction* );
  ~GammaRayTelTrackerSD();
  
  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step* astep,G4TouchableHistory* ROHist);
  void EndOfEvent(G4HCofThisEvent*);
  void clear();
  void DrawAll();
  void PrintAll();
  
private:
  
  GammaRayTelTrackerHitsCollection*  TrackerCollection;      
  GammaRayTelDetectorConstruction* Detector;
  G4int (*ThitXID)[30];
  G4int (*ThitYID)[30];
  G4int NbOfTKRLayers;
  G4int NbOfTKRStrips;

};

#endif






