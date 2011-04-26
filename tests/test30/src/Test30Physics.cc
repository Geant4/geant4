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
// -------------------------------------------------------------
//      GEANT 4 class for test30
//
//      History: based on object model of
//      ---------- Test30Physics -------
//                by Vladimir Ivanchenko, 12 March 2002 
// 
//    Modified:
//  11.10.2007 Added INCL cascade and RPG parameterized model (V.Ivanchenko)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Test30Physics.hh"

#include "G4UnitsTable.hh"
#include "Test30VSecondaryGenerator.hh"
#include "G4PhysicsTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4DecayPhysics.hh"

#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4RPGPiPlusInelastic.hh"
#include "G4RPGPiMinusInelastic.hh"
#include "G4RPGProtonInelastic.hh"
#include "G4RPGNeutronInelastic.hh"

#include "G4QStringChipsParticleLevelInterface.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4BinaryCascade.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4CascadeInterface.hh"
#include "G4InclAblaCascadeInterface.hh"
#include "G4InclLightIonInterface.hh"
#include "G4WilsonAbrasionModel.hh"
#include "G4QMDReaction.hh"

#include "G4TheoFSGenerator.hh"
#include "G4FTFModel.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4LundStringFragmentation.hh"

#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4QuasiElasticChannel.hh"

#include "G4QuasiElasticChannel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QStringChipsParticleLevelInterface.hh"

#include "HsQGSCInterface.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test30Physics::Test30Physics()
{
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test30Physics::~Test30Physics()
{
  //delete theDeExcitation;
  delete thePreCompound;
  delete theProcess;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test30Physics::Initialise()
{
  G4DecayPhysics dp;
  dp.ConstructParticle();
  theProcess = 0;
  theDeExcitation = 0;
  thePreCompound = 0;
  theQuasiElastic = 0;
  hkmod = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VProcess* Test30Physics::GetProcess(const G4String& gen_name,
		                      const G4ParticleDefinition* part,
				      G4Material* mat)
{
  G4cout <<  "Test30Physics entry" << G4endl;
  if(theProcess) { delete theProcess; }
  theProcess = 0;

  G4String part_name = part->GetParticleName();
  G4ProcessManager* man = new G4ProcessManager(part);
  if(!man) { return 0; }

  theProcess = new Test30HadronProduction();
 
  G4cout <<  "Process is created; gen= " << gen_name << G4endl;

  if(!theDeExcitation) {
    theDeExcitation = new G4ExcitationHandler();
    G4cout << "New deexcitation" << G4endl;
  }
  if(!thePreCompound) {
    thePreCompound = new G4PreCompoundModel(theDeExcitation);
    G4cout << "New pre-compound" << G4endl;
  }

  // Physics list for the given run
  Test30VSecondaryGenerator* sg = 0;

  // Choose generator
  if(gen_name == "lhep") {
    if(part_name == "proton")   
      sg = new Test30VSecondaryGenerator(new G4LEProtonInelastic(),mat);
    else if(part_name == "pi+") 
      sg = new Test30VSecondaryGenerator(new G4LEPionPlusInelastic(),mat);
    else if(part_name == "pi-") 
      sg = new Test30VSecondaryGenerator(new G4LEPionMinusInelastic(),mat);
    else if(part_name == "neutron") 
      sg = new Test30VSecondaryGenerator(new G4LENeutronInelastic(),mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "rpg") {
    if(part_name == "proton")   
      sg = new Test30VSecondaryGenerator(new G4RPGProtonInelastic(),mat);
    else if(part_name == "pi+") 
      sg = new Test30VSecondaryGenerator(new G4RPGPiPlusInelastic(),mat);
    else if(part_name == "pi-") 
      sg = new Test30VSecondaryGenerator(new G4RPGPiMinusInelastic(),mat);
    else if(part_name == "neutron") 
      sg = new Test30VSecondaryGenerator(new G4RPGNeutronInelastic(),mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "chips") {
    sg = new Test30VSecondaryGenerator(new G4QStringChipsParticleLevelInterface(),mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "preco") {
    sg = new Test30VSecondaryGenerator(thePreCompound,mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "binary") {
    G4BinaryCascade* hkm = new G4BinaryCascade();
    hkm->SetDeExcitation(thePreCompound);
    sg = new Test30VSecondaryGenerator(hkm, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "binary_ion") {
    G4BinaryLightIonReaction* hkm = new G4BinaryLightIonReaction();
    hkm->SetDeExcitation(theDeExcitation);
    hkm->SetPrecompound(thePreCompound);
    sg = new Test30VSecondaryGenerator(hkm, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "ftfp") {

    G4TheoFSGenerator* theModel = new G4TheoFSGenerator;
    G4FTFModel* theStringModel = new G4FTFModel();
    G4GeneratorPrecompoundInterface* theCascade = 
      new G4GeneratorPrecompoundInterface;
    theCascade->SetDeExcitation(thePreCompound);
    G4ExcitedStringDecay* theStringDecay = 
      new G4ExcitedStringDecay(new G4LundStringFragmentation());
    theStringModel->SetFragmentationModel(theStringDecay);

    theModel->SetTransport(theCascade);
    theModel->SetHighEnergyGenerator(theStringModel);
    theModel->SetMinEnergy(GeV);

    sg = new Test30VSecondaryGenerator(theModel, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "ftfc") {

    G4TheoFSGenerator* theModel = new G4TheoFSGenerator;
    G4FTFModel* theStringModel = new G4FTFModel();
    G4QStringChipsParticleLevelInterface* theCascade = 
      new G4QStringChipsParticleLevelInterface;
    G4ExcitedStringDecay* theStringDecay = 
      new G4ExcitedStringDecay(new G4LundStringFragmentation());
    theStringModel->SetFragmentationModel(theStringDecay);

    theModel->SetTransport(theCascade);
    theModel->SetHighEnergyGenerator(theStringModel);
    theModel->SetMinEnergy(GeV);

    sg = new Test30VSecondaryGenerator(theModel, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "ftfb") {

    G4TheoFSGenerator* theModel = new G4TheoFSGenerator();
    G4FTFModel* theStringModel= new G4FTFModel();
    G4BinaryCascade* theCascade = new G4BinaryCascade();
    theCascade->SetDeExcitation(thePreCompound);
    G4ExcitedStringDecay * theStringDecay = 
      new G4ExcitedStringDecay(new G4LundStringFragmentation());
    theModel->SetHighEnergyGenerator(theStringModel);
    theStringModel->SetFragmentationModel(theStringDecay);

    theModel->SetTransport(theCascade);
    theModel->SetHighEnergyGenerator(theStringModel);
    theModel->SetMinEnergy(GeV);

    sg = new Test30VSecondaryGenerator(theModel, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "qgsp") {
    G4TheoFSGenerator* theModel = new G4TheoFSGenerator();
    G4GeneratorPrecompoundInterface* theCascade = 
      new G4GeneratorPrecompoundInterface();
    theCascade->SetDeExcitation(thePreCompound);  
    G4QGSModel< G4QGSParticipants > * theStringModel = 
      new G4QGSModel< G4QGSParticipants >;
    G4ExcitedStringDecay* theQGStringDecay = 
      new G4ExcitedStringDecay(new G4QGSMFragmentation());
    theStringModel->SetFragmentationModel(theQGStringDecay);
    theModel->SetTransport(theCascade);
    theModel->SetHighEnergyGenerator(theStringModel);
    G4QuasiElasticChannel* theQuasiElastic = new G4QuasiElasticChannel;
    theModel->SetQuasiElasticChannel(theQuasiElastic);
    theModel->SetMinEnergy(0.5*GeV);
    theModel->SetMaxEnergy(100*TeV);
    sg = new Test30VSecondaryGenerator(theModel, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "qgsc") {
    HsQGSCInterface* qgsc = new HsQGSCInterface();
    qgsc->SetMinEnergy(GeV);
    sg = new Test30VSecondaryGenerator(qgsc, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "qgsb") {
    G4TheoFSGenerator* theModel = new G4TheoFSGenerator();
    G4BinaryCascade* theCascade = new G4BinaryCascade();
    theCascade->SetDeExcitation(thePreCompound);  
    G4QGSModel< G4QGSParticipants > * theStringModel = 
      new G4QGSModel< G4QGSParticipants >;
    G4ExcitedStringDecay* theQGStringDecay = 
      new G4ExcitedStringDecay(new G4QGSMFragmentation());
    theStringModel->SetFragmentationModel(theQGStringDecay);
    theModel->SetTransport(theCascade);
    theModel->SetHighEnergyGenerator(theStringModel);
    G4QuasiElasticChannel* theQuasiElastic = new G4QuasiElasticChannel;
    theModel->SetQuasiElasticChannel(theQuasiElastic);
    theModel->SetMinEnergy(0.5*GeV);
    theModel->SetMaxEnergy(100*TeV);
    sg = new Test30VSecondaryGenerator(theModel, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "bertini") {
    G4CascadeInterface* hkm = new G4CascadeInterface();
    sg = new Test30VSecondaryGenerator(hkm, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "bertini_preco") {
    G4CascadeInterface* hkm = new G4CascadeInterface();
    hkm->usePreCompoundDeexcitation();
    sg = new Test30VSecondaryGenerator(hkm, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "incl") {
    G4InclAblaCascadeInterface* hkm = new G4InclAblaCascadeInterface();
    sg = new Test30VSecondaryGenerator(hkm, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "incl_ion") {
    G4InclLightIonInterface* hkm = new G4InclLightIonInterface();
    sg = new Test30VSecondaryGenerator(hkm, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "abrasion") {
    G4WilsonAbrasionModel* wam = new G4WilsonAbrasionModel();
    size_t ne = mat->GetNumberOfElements();
    const G4ElementVector* ev = mat->GetElementVector();
    for(size_t i=0; i<ne; i++) {
      G4Element* elm = (*ev)[i];
      wam->ActivateFor(elm);
    }   
    sg = new Test30VSecondaryGenerator(wam, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);
    
  } else if(gen_name == "qmd") {
    G4QMDReaction* qmd = new G4QMDReaction();
    sg = new Test30VSecondaryGenerator(qmd, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);
    
  } else {
    G4cout << "WARNING: <" << gen_name << "> generator is unknown" << G4endl;
  }

  G4cout <<  "Secondary generator <"
         << gen_name << "> is initialized"
         << G4endl;
  return theProcess;

}	

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  






