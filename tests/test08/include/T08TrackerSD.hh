// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08TrackerSD.hh,v 1.1 1999-01-08 16:35:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef T08TrackerSD_h
#define T08TrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4ThreeVector.hh"
#include "T08TrackerHit.hh"
class G4Step;
class G4HCofThisEvent;

class T08TrackerSD : public G4VSensitiveDetector
{

  public:
      T08TrackerSD(G4String name);
      ~T08TrackerSD();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      T08TrackerHitsCollection *trackerCollection;

};




#endif

