// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyTrackerSD.hh,v 1.3 2000-05-26 13:11:39 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyTrackerSD_h
#define MyTrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"

#include "MyTrackerHit.hh"

class MyTrackerSD : public G4VSensitiveDetector
{

  public:
      MyTrackerSD(G4String);
      ~MyTrackerSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      MyTrackerHitsCollection* TrackerCollection;
};




#endif

