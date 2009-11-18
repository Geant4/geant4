//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// The code was written by :
//	^Claudio Andenna claudio.andenna@iss.infn.it, claudio.andenna@ispesl.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//
// ^ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#ifndef CML2PhysicsListH
#define CML2PhysicsListH

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LossTableManager.hh" //
#include "G4EmProcessOptions.hh"
#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4RayleighScattering.hh" 
#include "G4LivermoreRayleighModel.hh"
#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4eMultipleScattering.hh"
#include "G4UniversalFluctuation.hh"
#include "G4hMultipleScattering.hh" //
#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4eplusAnnihilation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4hBremsstrahlung.hh" // *
#include "G4hPairProduction.hh" // *
#include "G4hIonisation.hh" // *
#include "G4ionIonisation.hh" // *
#include "G4Gamma.hh" // *
#include "G4Electron.hh" // *
#include "G4Positron.hh" // *
#include "G4StepLimiter.hh"

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
