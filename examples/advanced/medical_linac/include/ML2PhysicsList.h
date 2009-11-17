#ifndef CML2PhysicsListH
#define CML2PhysicsListH

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LossTableManager.hh" //
#include "G4EmProcessOptions.hh" //
#include "G4ComptonScattering.hh" //
#include "G4GammaConversion.hh" //
#include "G4PhotoElectricEffect.hh" //
#include "G4LowEnergyRayleigh.hh" 
#include "G4eMultipleScattering.hh" //
#include "G4hMultipleScattering.hh" //
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4hBremsstrahlung.hh" // *
#include "G4hPairProduction.hh" // *
#include "G4hIonisation.hh" // *
#include "G4ionIonisation.hh" // *
#include "G4Gamma.hh" // *
#include "G4Electron.hh" // *
#include "G4Positron.hh" // *


#include "G4VUserPhysicsList.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"
#include "G4VUserPhysicsList.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessVector.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4UserLimits.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"


class G4LowEnergyIonisation;
class G4LowEnergyPhotoElectric;
class G4LowEnergyBremsstrahlung;

class CML2PhysicsList : public G4VUserPhysicsList
{
public:
	CML2PhysicsList(void);
	~CML2PhysicsList(void);
protected:
    void ConstructParticle();
    void ConstructProcess();
    void ConstructEM();
	void SetCuts();
	void ShowCutsValues();
};

#endif
