// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02TrackerSD.hh,v 1.1 1999-01-07 16:05:48 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef ExN02TrackerSD_h
#define ExN02TrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4ThreeVector.hh"
#include "ExN02TrackerHit.hh"
class G4Step;
class G4HCofThisEvent;

class ExN02TrackerSD : public G4VSensitiveDetector
{

  public:
      ExN02TrackerSD(G4String name);
      ~ExN02TrackerSD();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      ExN02TrackerHitsCollection *trackerCollection;

};




#endif

