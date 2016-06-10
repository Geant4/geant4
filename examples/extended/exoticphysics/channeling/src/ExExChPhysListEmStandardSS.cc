//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#include "ExExChPhysListEmStandardSS.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4LossTableManager.hh"
#include "G4EmProcessOptions.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4PEEffectFluoModel.hh"
#include "G4KleinNishinaModel.hh"
#include "G4LowEPComptonModel.hh"
#include "G4PenelopeGammaConversionModel.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"
#include "G4eSingleCoulombScatteringModel.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4Generator2BS.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4UniversalFluctuation.hh"

#include "G4eplusAnnihilation.hh"
#include "G4UAtomicDeexcitation.hh"


#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"

#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4ProcessManager.hh"

// Wrapper
#include "XWrapperDiscreteProcess.hh"
#include "XWrapperContinuousDiscreteProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExExChPhysListEmStandardSS::ExExChPhysListEmStandardSS(
                                                const G4String& name)
:G4VPhysicsConstructor(name){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExExChPhysListEmStandardSS::~ExExChPhysListEmStandardSS()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChPhysListEmStandardSS::ConstructProcess()
{
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    
    // Add standard EM Processes
    //
    aParticleIterator->reset();
    while( (*aParticleIterator)() ){
        G4ParticleDefinition* particle = aParticleIterator->value();
        G4String particleName = particle->GetParticleName();
        
        if (particleName == "gamma") {
            //G4cout << particleName << G4endl;
            // Compton scattering
            G4ComptonScattering* cs = new G4ComptonScattering;
            cs->SetEmModel(new G4KleinNishinaModel(),1);
            G4VEmModel* theLowEPComptonModel = new G4LowEPComptonModel();
            theLowEPComptonModel->SetHighEnergyLimit(20*MeV);
            cs->AddEmModel(0, theLowEPComptonModel);
            ph->RegisterProcess(cs, particle);

            // Photoelectric
            G4PhotoElectricEffect* pe = new G4PhotoElectricEffect();
            G4VEmModel* theLivermorePEModel =
                new G4LivermorePhotoElectricModel();
            theLivermorePEModel->SetHighEnergyLimit(10*GeV);
            pe->SetEmModel(theLivermorePEModel,1);
            ph->RegisterProcess(pe, particle);

            // Gamma conversion
            G4GammaConversion* gc = new G4GammaConversion();
            G4VEmModel* thePenelopeGCModel =
                new G4PenelopeGammaConversionModel();
            thePenelopeGCModel->SetHighEnergyLimit(1*GeV);
            gc->SetEmModel(thePenelopeGCModel,1);
            ph->RegisterProcess(gc, particle);
            
            // Rayleigh scattering
            ph->RegisterProcess(new G4RayleighScattering(), particle);
        } else if (particleName == "e-") {
            
            //G4cout << particleName << G4endl;
            // ionisation
            G4eIonisation* eIoni = new G4eIonisation();
            eIoni->SetStepFunction(0.1, 100*um);
            G4VEmModel* theIoniPenelope = new G4PenelopeIonisationModel();
            theIoniPenelope->SetHighEnergyLimit(0.1*MeV);
            eIoni->AddEmModel(0, theIoniPenelope, new G4UniversalFluctuation());
            
            XWrapperContinuousDiscreteProcess *eIoni_wrapper =
                new XWrapperContinuousDiscreteProcess();
            eIoni_wrapper->RegisterProcess(eIoni,-1);
            ph->RegisterProcess(eIoni_wrapper, particle);

            // bremsstrahlung
            G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
            XWrapperContinuousDiscreteProcess *eBrem_wrapper =
                new XWrapperContinuousDiscreteProcess();
            eBrem_wrapper->RegisterProcess(eBrem,-1);
            ph->RegisterProcess(eBrem_wrapper, particle);

            // coulomb scattering
            G4CoulombScattering* ecs = new G4CoulombScattering();
            ecs->SetBuildTableFlag(false);
            G4eSingleCoulombScatteringModel* ecsmodel =
                new G4eSingleCoulombScatteringModel();
            ecsmodel->SetPolarAngleLimit(0.0);
            ecs->AddEmModel(0, ecsmodel);
            XWrapperDiscreteProcess *ecs_wrapper =
                new XWrapperDiscreteProcess();
            ecs_wrapper->RegisterProcess(ecs,1,1);
            ph->RegisterProcess(ecs_wrapper, particle);


            // multiple scattering
            G4eMultipleScattering* ems = new G4eMultipleScattering();
            XWrapperContinuousDiscreteProcess *ems_wrapper =
            new XWrapperContinuousDiscreteProcess();
            ems_wrapper->RegisterProcess(ems,0,2);
            ph->RegisterProcess(ems_wrapper, particle);
            
        } else if (particleName == "e+") {
            //G4cout << particleName << G4endl;
            // ionisation
            G4eIonisation* eIoni = new G4eIonisation();
            eIoni->SetStepFunction(0.2, 100*um);
            G4VEmModel* theIoniPenelope = new G4PenelopeIonisationModel();
            theIoniPenelope->SetHighEnergyLimit(0.1*MeV);
            eIoni->AddEmModel(0, theIoniPenelope, new G4UniversalFluctuation());

            XWrapperContinuousDiscreteProcess *eIoni_wrapper =
                new XWrapperContinuousDiscreteProcess();
            eIoni_wrapper->RegisterProcess(eIoni,-1);
            ph->RegisterProcess(eIoni_wrapper, particle);

            // bremsstrahlung
            G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
            XWrapperContinuousDiscreteProcess *eBrem_wrapper =
                new XWrapperContinuousDiscreteProcess();
            eBrem_wrapper->RegisterProcess(eBrem,-1);
            ph->RegisterProcess(eBrem_wrapper, particle);

            // annihilation at rest and in flight
            G4eplusAnnihilation* eplusAnn = new G4eplusAnnihilation();
            XWrapperDiscreteProcess *eplusAnn_wrapper =
                new XWrapperDiscreteProcess();
            eplusAnn_wrapper->RegisterProcess(eplusAnn,-1);
            ph->RegisterProcess(eplusAnn_wrapper, particle);

            // coulomb scattering
            G4CoulombScattering* ecs = new G4CoulombScattering();
            ecs->SetBuildTableFlag(false);
            G4eSingleCoulombScatteringModel* ecsmodel =
                new G4eSingleCoulombScatteringModel();
            ecsmodel->SetPolarAngleLimit(0.0);
            ecs->AddEmModel(0, ecsmodel);
            XWrapperDiscreteProcess *ecs_wrapper =
                new XWrapperDiscreteProcess();
            ecs_wrapper->RegisterProcess(ecs,1,1);
            ph->RegisterProcess(ecs_wrapper, particle);

            // multiple scattering
            G4eMultipleScattering* ems = new G4eMultipleScattering();
            XWrapperContinuousDiscreteProcess *ems_wrapper =
            new XWrapperContinuousDiscreteProcess();
            ems_wrapper->RegisterProcess(ems,0,2);
            ph->RegisterProcess(ems_wrapper, particle);
            
        } else if ((particleName == "mu+" ||
                   particleName == "mu-")) {
            //G4cout << particleName << G4endl;
           // ionisation
            G4MuIonisation* muIoni = new G4MuIonisation();
            muIoni->SetStepFunction(0.2, 50*um);
            XWrapperContinuousDiscreteProcess *muIoni_wrapper =
                new XWrapperContinuousDiscreteProcess();
            muIoni_wrapper->RegisterProcess(muIoni,-1);
            ph->RegisterProcess(muIoni_wrapper, particle);

            // bremsstrahlung
            G4MuBremsstrahlung* muBrem = new G4MuBremsstrahlung();
            XWrapperContinuousDiscreteProcess *muBrem_wrapper =
                new XWrapperContinuousDiscreteProcess();
            muBrem_wrapper->RegisterProcess(muBrem,-1);
            ph->RegisterProcess(muBrem_wrapper, particle);

            // pair production
            G4MuPairProduction* muPair = new G4MuPairProduction();
            XWrapperContinuousDiscreteProcess* muPair_wrapper =
                new XWrapperContinuousDiscreteProcess();
            muPair_wrapper->RegisterProcess(muPair,-1);
            ph->RegisterProcess(muPair_wrapper, particle);

            // coulomb scattering
            G4CoulombScattering* ecs = new G4CoulombScattering();
            ecs->SetBuildTableFlag(false);
            G4eCoulombScatteringModel* ecsmodel =
                new G4eCoulombScatteringModel();
            ecsmodel->SetPolarAngleLimit(0.0);
            ecs->AddEmModel(0, ecsmodel);
            XWrapperDiscreteProcess *ecs_wrapper =
                new XWrapperDiscreteProcess();
            ecs_wrapper->RegisterProcess(ecs,1,1);
            ph->RegisterProcess(ecs_wrapper, particle);

            // multiple scattering
            XWrapperContinuousDiscreteProcess *mums_wrapper =
            new XWrapperContinuousDiscreteProcess();
            mums_wrapper->RegisterProcess(new G4MuMultipleScattering(),0,2);
            ph->RegisterProcess(mums_wrapper, particle);

        } else if ((particleName == "alpha" || particleName == "He3") ) {
            //G4cout << particleName << G4endl;
            // ionisation
            G4ionIonisation* ionIoni = new G4ionIonisation();
            ionIoni->SetStepFunction(0.1, 10*um);
            XWrapperContinuousDiscreteProcess *ionIoni_wrapper =
                new XWrapperContinuousDiscreteProcess();
            ionIoni_wrapper->RegisterProcess(ionIoni,-1);
            ph->RegisterProcess(ionIoni_wrapper, particle);

            // coulomb scattering
            G4CoulombScattering* ecs = new G4CoulombScattering();
            ecs->SetBuildTableFlag(false);
            XWrapperDiscreteProcess *ecs_wrapper =
                new XWrapperDiscreteProcess();
            ecs_wrapper->RegisterProcess(ecs,1,1);
            ph->RegisterProcess(ecs_wrapper, particle);


            // multiple scattering
            XWrapperContinuousDiscreteProcess *hms_wrapper =
            new XWrapperContinuousDiscreteProcess();
            hms_wrapper->RegisterProcess(new G4hMultipleScattering(),0,2);
            ph->RegisterProcess(hms_wrapper, particle);
            
        } else if( particleName == "proton" ||
                  particleName == "pi-" ||
                  particleName == "pi+"    ) {
            
            //G4cout << particleName << G4endl;
            // ionisation
            G4hIonisation* hIoni = new G4hIonisation();
            hIoni->SetStepFunction(0.1, 20*um);
            XWrapperContinuousDiscreteProcess *hIoni_wrapper =
            new XWrapperContinuousDiscreteProcess();
            hIoni_wrapper->RegisterProcess(hIoni,-1);
            ph->RegisterProcess(hIoni_wrapper, particle);

            // bremsstrahlung
            G4hBremsstrahlung* hBrem = new G4hBremsstrahlung();
            XWrapperContinuousDiscreteProcess *hBrem_wrapper =
            new XWrapperContinuousDiscreteProcess();
            hBrem_wrapper->RegisterProcess(hBrem,-1);
            ph->RegisterProcess(hBrem_wrapper, particle);

            // pair production
            G4hPairProduction* hPair = new G4hPairProduction();
            XWrapperContinuousDiscreteProcess* hPair_wrapper =
            new XWrapperContinuousDiscreteProcess();
            hPair_wrapper->RegisterProcess(hPair,-1);
            ph->RegisterProcess(hPair_wrapper, particle);

            // coulomb scattering
            G4CoulombScattering* ecs = new G4CoulombScattering();
            ecs->SetBuildTableFlag(false);
            XWrapperDiscreteProcess *ecs_wrapper =
                new XWrapperDiscreteProcess();
            ecs_wrapper->RegisterProcess(ecs,1,1);
            ph->RegisterProcess(ecs_wrapper, particle);

            // multiple scattering
            XWrapperContinuousDiscreteProcess *hms_wrapper =
            new XWrapperContinuousDiscreteProcess();
            hms_wrapper->RegisterProcess(new G4hMultipleScattering(),0,2);
            ph->RegisterProcess(hms_wrapper, particle);

        } else if (particleName == "GenericIon"  ) {
            //G4cout << particleName << G4endl;
            // ionisation
            G4ionIonisation* ionIoni = new G4ionIonisation();
            ionIoni->SetStepFunction(0.1, 1*um);
            XWrapperContinuousDiscreteProcess *ionIoni_wrapper =
                new XWrapperContinuousDiscreteProcess();
            ionIoni_wrapper->RegisterProcess(ionIoni,-1);
            ph->RegisterProcess(ionIoni_wrapper, particle);

            // coulomb scattering
            G4CoulombScattering* ecs = new G4CoulombScattering();
            ecs->SetBuildTableFlag(false);
            XWrapperDiscreteProcess *ecs_wrapper =
                new XWrapperDiscreteProcess();
            ecs_wrapper->RegisterProcess(ecs,1,1);
            ph->RegisterProcess(ecs_wrapper, particle);

            // multiple scattering
            XWrapperContinuousDiscreteProcess *hms_wrapper =
            new XWrapperContinuousDiscreteProcess();
            hms_wrapper->RegisterProcess(new G4hMultipleScattering(),0,2);
            ph->RegisterProcess(hms_wrapper, particle);

        } else if ((!particle->IsShortLived()) &&
                   (particle->GetPDGCharge() != 0.0) &&
                   (particle->GetParticleName() != "chargedgeantino")  ) {
            //G4cout << particleName << G4endl;
           //all others charged particles except geantino
            
            // ionisation
            G4hIonisation* hIoni = new G4hIonisation();
            XWrapperContinuousDiscreteProcess *hIoni_wrapper =
                new XWrapperContinuousDiscreteProcess();
            hIoni_wrapper->RegisterProcess(hIoni,-1);
            ph->RegisterProcess(hIoni_wrapper, particle);

            // coulomb scattering
            G4CoulombScattering* ecs = new G4CoulombScattering();
            ecs->SetBuildTableFlag(false);
            XWrapperDiscreteProcess *ecs_wrapper =
                new XWrapperDiscreteProcess();
            ecs_wrapper->RegisterProcess(ecs,1,1);
            ph->RegisterProcess(ecs_wrapper, particle);

            // multiple scattering
            XWrapperContinuousDiscreteProcess *hms_wrapper =
            new XWrapperContinuousDiscreteProcess();
            hms_wrapper->RegisterProcess(new G4hMultipleScattering(),0,2);
            ph->RegisterProcess(hms_wrapper, particle);
            
        }
    }
    
    // Em options
    //
    // Main options and setting parameters are shown here.
    // Several of them have default values.
    //
    G4EmProcessOptions emOptions;
    
    //physics tables
    //
    emOptions.SetMinEnergy(10*eV);
    emOptions.SetMaxEnergy(10*TeV);
    emOptions.SetDEDXBinning(12*20);
    emOptions.SetLambdaBinning(12*20);
    
    // scattering
    emOptions.SetPolarAngleLimit(0.0);
    
    // Deexcitation
    G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
    G4LossTableManager::Instance()->SetAtomDeexcitation(de);
    de->SetFluo(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

