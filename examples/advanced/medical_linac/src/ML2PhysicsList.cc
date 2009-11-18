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


#include "ML2PhysicsList.hh"


CML2PhysicsList::CML2PhysicsList(void):  G4VUserPhysicsList() 
{
	this->SetDefaultCutValue(1.0*mm);
	this->SetCutsWithDefault();
}

CML2PhysicsList::~CML2PhysicsList(void)
{}
void CML2PhysicsList::ConstructParticle()
{
	G4Electron::ElectronDefinition();
	G4Gamma::GammaDefinition();
	G4Positron::PositronDefinition();
}
void CML2PhysicsList::ConstructProcess()
{
	AddTransportation();
	ConstructEM();
}
void CML2PhysicsList::ConstructEM()
{

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {

	G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
	G4LivermorePhotoElectricModel* theLivermorePhotoElectricModel = new G4LivermorePhotoElectricModel();
	thePhotoElectricEffect->SetModel(theLivermorePhotoElectricModel);
	pmanager->AddDiscreteProcess(thePhotoElectricEffect);

	G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
	G4LivermoreComptonModel* theLivermoreComptonModel = new G4LivermoreComptonModel();
	theComptonScattering->SetModel(theLivermoreComptonModel);
	pmanager->AddDiscreteProcess(theComptonScattering);
	
	G4GammaConversion* theGammaConversion = new G4GammaConversion();
	G4LivermoreGammaConversionModel* theLivermoreGammaConversionModel = new G4LivermoreGammaConversionModel();
	theGammaConversion->SetModel(theLivermoreGammaConversionModel);
	pmanager->AddDiscreteProcess(theGammaConversion);
	
	G4RayleighScattering* theRayleigh = new G4RayleighScattering();
	G4LivermoreRayleighModel* theRayleighModel = new G4LivermoreRayleighModel();
	theRayleigh->SetModel(theRayleighModel);
	pmanager->AddDiscreteProcess(theRayleigh);


    } else if (particleName == "e-") {

      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      pmanager->AddProcess(msc,                   -1, 1, 1);
      
      // Ionisation
      G4eIonisation* eIoni = new G4eIonisation();
      G4LivermoreIonisationModel* theLivermoreIonisationModel = new G4LivermoreIonisationModel();
      eIoni->SetEmModel(theLivermoreIonisationModel);
      eIoni->SetStepFunction(0.2, 100*um); //     
      pmanager->AddProcess(eIoni,                 -1, 2, 2);
      
      // Bremsstrahlung
      G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
      eBrem->SetEmModel(new G4LivermoreBremsstrahlungModel());
      pmanager->AddProcess(eBrem,         -1,-3, 3);

      pmanager->AddProcess(new G4StepLimiter(), -1, -1, 4);

    } else if (particleName == "e+") {

      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      pmanager->AddProcess(msc,                   -1, 1, 1);
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.2, 100*um);      
      pmanager->AddProcess(eIoni,                 -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung, -1,-3, 3);
      pmanager->AddProcess(new G4eplusAnnihilation,0,-1, 4);

	}
	}
}
void CML2PhysicsList::ShowCutsValues()
{
	G4Region *reg;
	G4String regName;
	G4double cutValue;
	G4RegionStore *regStore=G4RegionStore::GetInstance();
	G4int nRegions=regStore->size();
	std::vector <G4Region*>::iterator reg_Iter;
	reg_Iter=regStore->begin();
	for (int i=0; i< nRegions;i++)
	{
		reg=reg_Iter[i];
		regName = reg->GetName();
		cutValue=reg->GetProductionCuts()->GetProductionCut(0);
		std::cout << regName<<" cut Value: " << cutValue/mm<<" [mm]"<< G4endl;
	}
}
void CML2PhysicsList::SetCuts()
{
	this->ShowCutsValues();
}

