#ifndef BrachyRunAction_h
#define BrachyRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Timer;
class G4Run;
class BrachyEventAction;
class BrachyDetectorConstruction;
class BrachyAnalysisManager;
class BrachyRunAction : public G4UserRunAction
{
  public:
    BrachyRunAction(G4String& );
   ~BrachyRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run* );
  private:
  G4Timer* timer;
   G4String      SDname;
  G4float* pVoxel;
  G4int  * NumberHits;
   G4int m_NumVoxelX;
   G4int m_NumVoxelZ;
  G4int j;
   G4double VoxelWidth_X;
   G4double VoxelWidth_Z;
  G4double x,z;
  const BrachyEventAction* pEvent;
   BrachyDetectorConstruction *pDetector;
 
};

#endif


