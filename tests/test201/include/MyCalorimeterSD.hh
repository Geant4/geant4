// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyCalorimeterSD.hh,v 1.3 2000-05-26 13:11:39 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyCalorimeterSD_h
#define MyCalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"

#include "MyCalorimeterHit.hh"

class MyCalorimeterSD : public G4VSensitiveDetector
{

  public:
      MyCalorimeterSD(G4String);
      ~MyCalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      MyCalorimeterHitsCollection *CalCollection;

      int CellID[3];
};




#endif

