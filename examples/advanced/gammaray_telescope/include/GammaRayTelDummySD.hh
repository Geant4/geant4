
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelDummySD.hh,v 1.1 2000-11-15 20:27:39 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelDummySD  ------
//           by F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************
//
// Dummy sensitive used only to flag sensitivity
// in cells of RO geometry.
//

#ifndef GammaRayTelDummySD_h
#define GammaRayTelDummySD_h 1

#include "G4VSensitiveDetector.hh"
class G4Step;

class GammaRayTelDummySD : public G4VSensitiveDetector
{
public:
  GammaRayTelDummySD();
  ~GammaRayTelDummySD() {};
  
  void Initialize(G4HCofThisEvent*HCE) {};
  G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist) {return false;}
  void EndOfEvent(G4HCofThisEvent*HCE) {};
  void clear() {};
  void DrawAll() {};
  void PrintAll() {};
};
GammaRayTelDummySD::GammaRayTelDummySD()
  : G4VSensitiveDetector("dummySD")
{}
#endif
