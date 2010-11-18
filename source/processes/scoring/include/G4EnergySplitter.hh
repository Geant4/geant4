#ifndef G4EnergySplitter_h
#define G4EnergySplitter_h 1

////////////////////////////////////////////////////////////////////////////////
// (Description)
// 
// Class to calculate the split of energy in voxels when G4RegularNavigation is used
// It takes into account energy loss and multiple scattering corrections as the particles loses energies from voxel to voxel
//
// Created: 2010-11-09 Pedro Arce
// 
///////////////////////////////////////////////////////////////////////////////
#include "globals.hh"
class G4PhantomParameterisation;
class G4EnergyLossForExtrapolator;
class G4VPhysicalVolume;
class G4Material;
class G4Step;
#include <vector>

class G4EnergySplitter 
{
public: // with description
  G4EnergySplitter();
  virtual ~G4EnergySplitter();
  
  G4int SplitEnergyInVolumes(const G4Step* aStep );    // Calculates the energy spliting, and dumps it into theEnergies. Returns number of steps

  inline void GetLastVoxelID( G4int& voxelID);
  inline void GetFirstVoxelID( G4int& voxelID);
  inline void GetVoxelID( G4int stepNo, G4int& voxelID );
  inline void GetVoxelIDAndLength( G4int stepNo, G4int& voxelID, G4double& stepLength ); 
  inline void GetLengthAndEnergyDeposited( G4int stepNo, G4int& voxelID, G4double& stepLength, G4double &energyLoss);
  inline void GetLengthAndInitialEnergy( G4double &preStepEnergy, G4int stepNo, G4int& voxelID, G4double& stepLength, G4double &initialEnergy);

  inline void SetNIterations( G4int niter );

  inline G4Material* GetVoxelMaterial( G4int stepNo );

private:
  inline void GetStepLength( G4int stepNo, G4double& stepLength );

  void GetPhantomParam(G4bool mustExist);
  G4bool IsPhantomVolume( G4VPhysicalVolume* pv );

  G4EnergyLossForExtrapolator* theElossExt;
  
  G4int theNIterations;

  std::vector<G4double> theEnergies;

  G4PhantomParameterisation* thePhantomParam;

};

#include "G4EnergySplitter.icc"

#endif
