//
//    **********************************
//    *                                *
//    *      ThyroidEventAction.hh     *
//    *                                *
//    **********************************

#ifndef ThyroidEventAction_h
#define ThyroidEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class ThyroidEventAction : public G4UserEventAction
{
  public:
    ThyroidEventAction();
    ~ThyroidEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

 
  //G4float *m_pVoxel;

  private:
    G4double  EnergyDep;
  G4int NumVoxelX;
  G4int NumVoxelY;
  G4int NumVoxelZ;
  // ThyroidDetectorConstruction *pDetector;
  G4double VoxelWidth_X;
  G4double VoxelWidth_Y;
  G4double VoxelWidth_Z;

 G4int m_HitsCollectionID;
 G4int i;
  G4int w;
  G4int k;
  G4double x;
  G4double y;
  G4double z;
  G4String      SDname;
};

#endif

    




