// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyTrackerSD.hh,v 1.1 1999-01-08 16:35:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyTrackerSD_h
#define MyTrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "MyTrackerHitsCollection.hh"
#include "MyTrackerHit.hh"
#include "G4Step.hh"

class MyTrackerSD : public G4VSensitiveDetector
{

  public:
      MyTrackerSD(G4String name);
      ~MyTrackerSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      MyTrackerHitsCollection *TrackerCollection;
};




#endif

