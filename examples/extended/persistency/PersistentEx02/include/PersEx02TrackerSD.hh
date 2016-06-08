// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: PersEx02TrackerSD.hh,v 1.3 1999/11/29 18:23:33 morita Exp $
// GEANT4 tag $Name: geant4-01-01 $
//

#ifndef PersEx02TrackerSD_h
#define PersEx02TrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "PersEx02TrackerHitsCollection.hh"
#include "G4Step.hh"
#include "HepODBMS/odbms/HepODBMS.h"

class G4PersistentHitMan;

class PersEx02TrackerSD : public G4VSensitiveDetector
{

  public:
      PersEx02TrackerSD(G4String name);
      ~PersEx02TrackerSD();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      G4PersistentHitMan* f_hitMan;
      HepRef(PersEx02TrackerHitsCollection) EvenCollection;
      HepRef(PersEx02TrackerHitsCollection)  OddCollection;
      HepContainerRef hitContainer;
};

#endif

