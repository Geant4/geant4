// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20DummySD.hh,v 1.1 2001-05-24 19:49:21 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ Tst20DummySD  ------
// ************************************************************
//
// Dummy sensitive used only to flag sensitivity
// in cells of RO geometry.
//

#ifndef Tst20DummySD_h
#define Tst20DummySD_h 1

#include "G4VSensitiveDetector.hh"
class G4Step;

class Tst20DummySD : public G4VSensitiveDetector
{
public:
  Tst20DummySD();
  ~Tst20DummySD() {};
  
  void Initialize(G4HCofThisEvent*HCE) {};
  G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist) {return false;}
  void EndOfEvent(G4HCofThisEvent*HCE) {};
  void clear() {};
  void DrawAll() {};
  void PrintAll() {};
};
Tst20DummySD::Tst20DummySD()
  : G4VSensitiveDetector("dummySD")
{}
#endif
