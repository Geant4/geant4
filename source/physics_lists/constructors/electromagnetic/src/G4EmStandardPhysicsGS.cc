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
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmStandardPhysicsGS
//
// Author:      V.Ivanchenko 05.06.2015
//
// Modified:
//
// 25.10.25    Changed the settings to the most accurate configuration of the
//             Goudsmit-Saunderson MSC model for e-/e+ multiple Coulomb
//             scattering. The class descrition has been changed to reflect the
//             new configuration (Mihaly Novak).
//
// Class Description:
//
// This EM physics constructor utilises the Goudsmit-Saunderson (GS) MSC model
// for e-/e+ multiple Coulomb scattering (below 1 GeV kinetic energies). The
// GS MSC model has been changed in version 11.4 keeping only its accurate
// stepping and boundary crossing algorithms while removeing the other, less
// accurate alternatives. All the corrections offered by the GS model, including
// the Mott, screening and scattering power corrections, are activated. This,
// together with the Penelope model for e-/e+ ionisations (below 1 GeV),
// offers an accurate an accurate e-/e+ simualtion down to few keV kinetic
// enegies independently form the target material and geometry.
//
// The same settings of the GS MSC model has already been used for e-/e+ below
// 100 MeV kinetic energies in the option4, Penelope and Livermore EM physics
// constructors since version 10.6.
//
//----------------------------------------------------------------------------
//

#include "G4EmStandardPhysicsGS.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4LossTableManager.hh"
#include "G4EmParameters.hh"
#include "G4EmStandUtil.hh"
#include "G4EmBuilder.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4RayleighScattering.hh"

#include "G4hMultipleScattering.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"
#include "G4WentzelVIModel.hh"
#include "G4GoudsmitSaundersonMscModel.hh"

#include "G4eIonisation.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4eBremsstrahlung.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4Generator2BS.hh"

#include "G4eplusAnnihilation.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4GenericIon.hh"

#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4EmModelActivator.hh"
#include "G4GammaGeneralProcess.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmStandardPhysicsGS);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmStandardPhysicsGS::G4EmStandardPhysicsGS(G4int ver, const G4String&)
  : G4VPhysicsConstructor("G4EmStandardGS")
{
  SetVerboseLevel(ver);
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetVerbose(ver);
  // use a denser discrete kinetic energy grid for more accurate interpolation
  param->SetNumberOfBinsPerDecade(16);
  // set the continuous step limit to: 0.2*Range that goes to Range below 10 um
  param->SetStepFunction(0.2, 10*CLHEP::um);
  // set the GS MSC model for e-/e+ to be used below 1.0 GeV with its (Mott,
  // screening, scattering power) corrections activated and with the accurate
  // stepping and boundary crossing algorithms (no other options since 11.4)
  // with a skin of 3 elastic MFP near boundary
  param->SetMscEnergyLimit(1.0*CLHEP::GeV);
  param->SetUseMottCorrection(true);
  param->SetMscStepLimitType(fUseSafetyPlus);
  param->SetMscSkin(3);
  param->SetMscRangeFactor(0.08);
  // activate fluoresence, i.e. emission of characteristic X-ray
  param->SetFluo(true);
  // set the energy loss fluctuation type
  param->SetFluctuationType(fUrbanFluctuation);
  SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmStandardPhysicsGS::~G4EmStandardPhysicsGS()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmStandardPhysicsGS::ConstructParticle()
{
  // minimal set of particles for EM physics
  G4EmBuilder::ConstructMinimalEmSet();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmStandardPhysicsGS::ConstructProcess()
{
  if(verboseLevel > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4EmBuilder::PrepareEMPhysics();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // processes used by several particles
  G4hMultipleScattering* hmsc = new G4hMultipleScattering("ionmsc");
  G4NuclearStopping* pnuc(nullptr);

  // high energy limit for e+- scattering models and bremsstrahlung
  G4double mscEnergyLimit = G4EmParameters::Instance()->MscEnergyLimit();

  // Add gamma EM processes

  // gamma
  G4ParticleDefinition* particle = G4Gamma::Gamma();

  G4PhotoElectricEffect* pe = new G4PhotoElectricEffect();
  pe->SetEmModel(new G4LivermorePhotoElectricModel());

  G4ComptonScattering* cs = new G4ComptonScattering;
  G4GammaConversion* gc = new G4GammaConversion;
  G4RayleighScattering* rs = new G4RayleighScattering;

  if (G4EmParameters::Instance()->GeneralProcessActive()) {
    G4GammaGeneralProcess* sp = new G4GammaGeneralProcess();
    sp->AddEmProcess(pe);
    sp->AddEmProcess(cs);
    sp->AddEmProcess(gc);
    sp->AddEmProcess(rs);
    G4LossTableManager::Instance()->SetGammaGeneralProcess(sp);
    ph->RegisterProcess(sp, particle);
  } else {
    ph->RegisterProcess(pe, particle);
    ph->RegisterProcess(cs, particle);
    ph->RegisterProcess(gc, particle);
    ph->RegisterProcess(rs, particle);
  }

  // e-
  particle = G4Electron::Electron();

  // msc: GS[:100 MeV] + WentzelVI[100 MeV:]
  G4GoudsmitSaundersonMscModel* msc1 = new G4GoudsmitSaundersonMscModel();
  G4WentzelVIModel* msc2 = new G4WentzelVIModel();
  msc1->SetHighEnergyLimit(mscEnergyLimit);
  msc2->SetLowEnergyLimit(mscEnergyLimit);
  G4EmBuilder::ConstructElectronMscProcess(msc1, msc2, particle);
  // (WVI is a mixed model, i.e. needs single scattering)
  G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel();
  G4CoulombScattering* ss = new G4CoulombScattering();
  ss->SetEmModel(ssm);
  ss->SetMinKinEnergy(mscEnergyLimit);
  ssm->SetLowEnergyLimit(mscEnergyLimit);
  ssm->SetActivationLowEnergyLimit(mscEnergyLimit);

  // ionisation: Penelope[:1.0 GeV] + Moller[1.0 GeV:]
  G4eIonisation* eioni = new G4eIonisation();
  eioni->SetFluctModel(G4EmStandUtil::ModelOfFluctuations());
  G4VEmModel* theIoniMod = new G4PenelopeIonisationModel();
  theIoniMod->SetHighEnergyLimit(1.0*CLHEP::GeV);
  eioni->AddEmModel(0, theIoniMod);

  // bremsstrahlung: Seltzer-Berger[:1.0 GeV] + extended Bethe–Heitler[1.0 GeV:]
  G4eBremsstrahlung* brem = new G4eBremsstrahlung();
  G4SeltzerBergerModel* br1 = new G4SeltzerBergerModel();
  G4eBremsstrahlungRelModel* br2 = new G4eBremsstrahlungRelModel();
  br1->SetAngularDistribution(new G4Generator2BS());
  br2->SetAngularDistribution(new G4Generator2BS());
  brem->SetEmModel(br1);
  brem->SetEmModel(br2);
  br1->SetHighEnergyLimit(1.0*CLHEP::GeV);

  ph->RegisterProcess(eioni, particle);
  ph->RegisterProcess(brem, particle);
  ph->RegisterProcess(ss, particle);

  // e+
  particle = G4Positron::Positron();

  // msc: GS[:100 MeV] + WentzelVI[100 MeV:]
  msc1 = new G4GoudsmitSaundersonMscModel();
  msc2 = new G4WentzelVIModel();
  msc1->SetHighEnergyLimit(mscEnergyLimit);
  msc2->SetLowEnergyLimit(mscEnergyLimit);
  G4EmBuilder::ConstructElectronMscProcess(msc1, msc2, particle);
  // (WVI is a mixed model, i.e. needs single scattering)
  ssm = new G4eCoulombScatteringModel();
  ss = new G4CoulombScattering();
  ss->SetEmModel(ssm);
  ss->SetMinKinEnergy(mscEnergyLimit);
  ssm->SetLowEnergyLimit(mscEnergyLimit);
  ssm->SetActivationLowEnergyLimit(mscEnergyLimit);

  // ionisation: Penelope[:1.0 GeV] + Bhabha[1.0 GeV:]
  eioni = new G4eIonisation();
  eioni->SetFluctModel(G4EmStandUtil::ModelOfFluctuations());
  G4VEmModel* pen = new G4PenelopeIonisationModel();
  pen->SetHighEnergyLimit(1.0*CLHEP::GeV);
  eioni->AddEmModel(0, pen);

  // bremsstrahlung: Seltzer-Berger[:1.0 GeV] + extended Bethe–Heitler[1.0 GeV:]
  brem = new G4eBremsstrahlung();
  br1 = new G4SeltzerBergerModel();
  br2 = new G4eBremsstrahlungRelModel();
  br1->SetAngularDistribution(new G4Generator2BS());
  br2->SetAngularDistribution(new G4Generator2BS());
  brem->SetEmModel(br1);
  brem->SetEmModel(br2);
  br1->SetHighEnergyLimit(1.0*CLHEP::GeV);

  ph->RegisterProcess(eioni, particle);
  ph->RegisterProcess(brem, particle);
  ph->RegisterProcess(new G4eplusAnnihilation(), particle);
  ph->RegisterProcess(ss, particle);

  // generic ion
  particle = G4GenericIon::GenericIon();
  G4ionIonisation* ionIoni = new G4ionIonisation();
  ph->RegisterProcess(hmsc, particle);
  ph->RegisterProcess(ionIoni, particle);

  // muons, hadrons, ions
  G4EmBuilder::ConstructCharged(hmsc, pnuc);

  // extra configuration
  G4EmModelActivator mact(GetPhysicsName());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
