#ifndef inputDataH
#define inputDataH

#include "globals.hh"
#include <vector>
#include "G4RotationMatrix.hh"

enum idParticleSource
{
	id_randomTarget=1,
	id_phaseSpace=2
};
enum idTypeOfSensitiveDetector
{
	idSD_ComponentROG=1,
	idSD_PhaseSpace=2,
	idSD_KillerPlane=3
};
struct SStartInputData
{
	G4String fileInputData;
	G4int seed;
};
struct SGeneralData
{
	G4String WorldName, fileExperimentalData, StartFileInputData;
	G4int seed, nBeam, nMaxParticlesInRamPlanePhaseSpace;
	G4bool bSaveROG, bCompareExp;
	G4String PhaseSpaceOutFile, ROGOutFile;

	G4bool bSavePhaseSpace;
	G4bool bStopAtPhaseSpace;
	G4ThreeVector centrePhaseSpace, halfSizePhaseSpace;

	G4int minNumberOfEvents;
	int saving_in_Selected_Voxels_every_events;
	int saving_in_ROG_Voxels_every_events;
	int max_N_particles_in_PhSp_File; 
};

struct Sparticle
{
	G4ThreeVector pos, dir;
	G4double kinEnergy;
	G4int nPrimaryPart, primaryParticlePDGE, partPDGE, volumeId;
	G4String volumeName;
};
struct SPrimaryParticle
{
	G4int partPDGE, nPrimaryParticle;
	G4int nParticlesInPhSp;
};
struct SInputData
{
	SGeneralData generalData;
	SPrimaryParticle primaryParticleData;
};
struct Svoxel
{
	G4ThreeVector pos, halfSize;
	G4double depEnergy, depEnergy2, expDose, depEnergyNorm, depEnergyNormError;
	G4int nEvents;
	G4String volumeName;
};
#endif
