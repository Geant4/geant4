// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyCalorimeterSD.hh,v 1.1 1999-01-08 16:35:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyCalorimeterSD_h
#define MyCalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "MyCalorimeterHitsCollection.hh"
#include "MyCalorimeterHit.hh"
#include "G4Step.hh"

class MyCalorimeterSD : public G4VSensitiveDetector
{

  public:
      MyCalorimeterSD(G4String name);
      ~MyCalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      MyCalorimeterHitsCollection *CalCollection;

      int CellID[3];
};




#endif

