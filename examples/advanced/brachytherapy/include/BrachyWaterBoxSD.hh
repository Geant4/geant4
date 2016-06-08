//    ********************************
//    *                              *  
//    *     BrachyWaterBoxSD.hh      *
//    *                              *
//    ********************************

#ifndef BrachyWaterBoxSD_h
#define BrachyWaterBoxSD_h 1

#include "G4VSensitiveDetector.hh"
#include "BrachyWaterBoxHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class BrachyWaterBoxSD : public G4VSensitiveDetector
{
 public:
      BrachyWaterBoxSD(G4String name, G4int NumVoxelX, G4int NumVoxelZ);
      ~BrachyWaterBoxSD();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();
 private:
      BrachyWaterBoxHitsCollection *m_pWaterBoxHitsCollection;
	
      const G4int m_NumVoxelX;
      const G4int m_NumVoxelZ;
      G4int *m_pVoxelID;
};
#endif

