//    **********************************
//    *                                *
//    *      BrachyEventAction.hh      *
//    *                                *
//    **********************************

#ifndef BrachyEventAction_h
#define BrachyEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class BrachyEventAction : public G4UserEventAction
{
  public:
    BrachyEventAction(float *pVoxel,G4int NumVoxelX,G4int NumVoxelZ);
    ~BrachyEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

  public:
    const G4int m_NumVoxelX;
    const G4int m_NumVoxelZ;
    float *m_pVoxel;

  private:
    G4int m_HitsCollectionID;
};

#endif

    

