// Rich advanced example for Geant4
// RichTbSD.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbSD_h
#define RichTbSD_h 1
#include "globals.hh"
#include "G4VSensitiveDetector.hh"
#include "G4ThreeVector.hh"
#include "RichTbHit.hh"
#include "RichTbGeometryParameters.hh"
class G4Step;
class G4HCofThisEvent;


class RichTbSD : public G4VSensitiveDetector
{

  public:
      RichTbSD(G4String , G4int, G4int, G4int,G4String );
      virtual ~RichTbSD();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      RichTbHitsCollection * RichTbHitCollection;
      G4int NumberOfSensitiveHpds;
      G4int NumberOfSensitiveSectorsInHpd;
      G4int NumberOfSensitivePixelsInSector;
      vector<G4int> HpdSDID;

     
      G4int HCID;
};



#endif
