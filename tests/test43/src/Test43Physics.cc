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
//
// -------------------------------------------------------------
//      GEANT 4 class 
//
//      ---------- Test43Physics -------
//           created from test30 files originally by Vladimir Ivanchenko
// 
//    Modified:
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Test43Physics.hh"

#include "G4UnitsTable.hh"
#include "Test43VSecondaryGenerator.hh"
#include "G4PhysicsTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4PionZero.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Neutron.hh"
#include "G4AntiNeutron.hh"
#include "G4NeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4GenericIon.hh"
#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4LEKaonMinusInelastic.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4LEKaonZeroSInelastic.hh"



#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4HEProtonInelastic.hh"
#include "G4HENeutronInelastic.hh"
#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4HEProtonInelastic.hh"
#include "G4HEKaonPlusInelastic.hh"
#include "G4HEKaonMinusInelastic.hh"
#include "G4HEKaonZeroLongInelastic.hh"
#include "G4HEKaonZeroShortInelastic.hh"

//#include "G4PionMinusNuclearReaction.hh" 
#include "G4StringChipsInterface.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4BinaryCascade.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4CascadeInterface.hh"
#include "G4WilsonAbrasionModel.hh"

#include "G4TheoFSGenerator.hh"
#include "G4FTFModel.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4LundStringFragmentation.hh"

#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4StringChipsParticleLevelInterface.hh"
//#include ".hh"

#include "G4ElasticHadrNucleusHE.hh"
#include "G4LElastic.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test43Physics::Test43Physics()
{
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test43Physics::~Test43Physics()
{
  //  if(theProcess) delete theProcess;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test43Physics::Initialise()
{
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();

	  // Strange barions
  G4Lambda::LambdaDefinition();
  G4AntiLambda::AntiLambdaDefinition();
  G4SigmaZero::SigmaZeroDefinition();
  G4AntiSigmaZero::AntiSigmaZeroDefinition();
  G4SigmaPlus::SigmaPlusDefinition();
  G4AntiSigmaPlus::AntiSigmaPlusDefinition();
  G4SigmaMinus::SigmaMinusDefinition();
  G4AntiSigmaMinus::AntiSigmaMinusDefinition();
  G4XiZero::XiZeroDefinition();
  G4AntiXiZero::AntiXiZeroDefinition();
  G4XiMinus::XiMinusDefinition();
  G4AntiXiMinus::AntiXiMinusDefinition();
  G4OmegaMinus::OmegaMinusDefinition();
  G4AntiOmegaMinus::AntiOmegaMinusDefinition();


  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();

  G4GenericIon::GenericIonDefinition();
  G4Deuteron::DeuteronDefinition();
  G4Alpha::AlphaDefinition();
  G4Triton::TritonDefinition();
  theProcess = 0;
  theDeExcitation = 0;
  thePreCompound = 0;
  hkmod = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VProcess* Test43Physics::GetProcess(const G4String& gen_name,
		                      const G4String& part_name,
		                            G4Material* mat)
{
  G4cout <<  "Test43Physics entry" << G4endl;
  if(theProcess) delete theProcess;
  theProcess = 0;

  G4ProcessManager* man = 0;

  if(part_name == "proton")   man = new G4ProcessManager(G4Proton::Proton());
  else if(part_name == "pi+") man = new G4ProcessManager(G4PionPlus::PionPlus());
  else if(part_name == "pi-") man = new G4ProcessManager(G4PionMinus::PionMinus());
  else if(part_name == "neutron") man = new G4ProcessManager(G4Neutron::Neutron());
  else if(part_name == "deuteron") man = new G4ProcessManager(G4Deuteron::Deuteron());
  else if(part_name == "alpha") man = new G4ProcessManager(G4Alpha::Alpha());
  else if(part_name == "GenericIon") man = new G4ProcessManager(G4GenericIon::GenericIon());
  else if(part_name == "kaon+") man = new G4ProcessManager(G4KaonPlus::KaonPlus());
  else if(part_name == "kaon-") man = new G4ProcessManager(G4KaonMinus::KaonMinus());
  else if(part_name == "kaon0L") man = new G4ProcessManager(G4KaonZeroLong::KaonZeroLong());
  else if(part_name == "kaon0S") man = new G4ProcessManager(G4KaonZeroShort::KaonZeroShort());
  else if(part_name == "anti_proton") man = new G4ProcessManager(G4AntiProton::AntiProton());
  else {
       G4cerr << "Test43Physics::GetProcess() - Fatal error: particle "<< part_name << " not implemeted" << G4endl;
       G4Exception("Test43Physics::GetProcess() " ) ;
  }

  if(!man) return 0;

  theProcess = new Test43HadronProduction();
  G4cout <<  "Process is created; gen= " << gen_name << G4endl;

  // Physics list for the given run
  Test43VSecondaryGenerator* sg = 0;

  // Choose generator

  if(gen_name == "lepar") {
     if(part_name == "proton")   sg = new Test43VSecondaryGenerator(new G4LEProtonInelastic(),mat);
     else if(part_name == "pi+") sg = new Test43VSecondaryGenerator(new G4LEPionPlusInelastic(),mat);
     else if(part_name == "pi-") sg = new Test43VSecondaryGenerator(new G4LEPionMinusInelastic(),mat);
     else if(part_name == "neutron") sg = new Test43VSecondaryGenerator(new G4LENeutronInelastic(),mat);
     else if(part_name == "kaon+")  sg = new Test43VSecondaryGenerator(new G4LEKaonPlusInelastic(),mat);
     else if(part_name == "kaon-")  sg = new Test43VSecondaryGenerator(new G4LEKaonMinusInelastic(),mat);
     else if(part_name == "kaon0L") sg = new Test43VSecondaryGenerator(new G4LEKaonZeroLInelastic(),mat);
     else if(part_name == "kaon0S") sg = new Test43VSecondaryGenerator(new G4LEKaonZeroSInelastic(),mat);
     theProcess->SetSecondaryGenerator(sg);
     man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "hepar") {
     if(part_name == "proton")   sg = new Test43VSecondaryGenerator(new G4HEProtonInelastic(),mat);
     else if(part_name == "pi+") sg = new Test43VSecondaryGenerator(new G4HEPionPlusInelastic(),mat);
     else if(part_name == "pi-") sg = new Test43VSecondaryGenerator(new G4HEPionMinusInelastic(),mat);
     else if(part_name == "neutron") sg = new Test43VSecondaryGenerator(new G4HENeutronInelastic(),mat);
     else if(part_name == "kaon+")  sg = new Test43VSecondaryGenerator(new G4HEKaonPlusInelastic(),mat);
     else if(part_name == "kaon-")  sg = new Test43VSecondaryGenerator(new G4HEKaonMinusInelastic(),mat);
     else if(part_name == "kaon0L") sg = new Test43VSecondaryGenerator(new G4HEKaonZeroLongInelastic(),mat);
     else if(part_name == "kaon0S") sg = new Test43VSecondaryGenerator(new G4HEKaonZeroShortInelastic(),mat);
     theProcess->SetSecondaryGenerator(sg);
     man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "QGSP") {
     G4TheoFSGenerator * model = new G4TheoFSGenerator;

     G4QGSModel< G4QGSParticipants > *theStringModel = new G4QGSModel< G4QGSParticipants >;
     G4ExcitedStringDecay *theStringDecay = new G4ExcitedStringDecay(new G4QGSMFragmentation);
     theStringModel->SetFragmentationModel(theStringDecay);

     G4GeneratorPrecompoundInterface *theCascade = new G4GeneratorPrecompoundInterface;
     G4PreCompoundModel *thePreEquilib = new G4PreCompoundModel(new G4ExcitationHandler);
     theCascade->SetDeExcitation(thePreEquilib);  

     model->SetTransport(theCascade);
     model->SetHighEnergyGenerator(theStringModel);
  
    sg = new Test43VSecondaryGenerator(model,mat);

    //G4cout <<  "Generator is ready" << G4endl;
    theProcess->SetSecondaryGenerator(sg);
    //G4cout <<  "Generator is set" << G4endl;
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "QGSC") {
     G4TheoFSGenerator * model = new G4TheoFSGenerator;

     G4QGSModel< G4QGSParticipants > *theStringModel = new G4QGSModel< G4QGSParticipants >;
     G4ExcitedStringDecay *theStringDecay = new G4ExcitedStringDecay(new G4QGSMFragmentation);
     theStringModel->SetFragmentationModel(theStringDecay);

     G4StringChipsParticleLevelInterface *theCascade = new G4StringChipsParticleLevelInterface;

     model->SetTransport(theCascade);
     model->SetHighEnergyGenerator(theStringModel);
  
    sg = new Test43VSecondaryGenerator(model,mat);

    //G4cout <<  "Generator is ready" << G4endl;
    theProcess->SetSecondaryGenerator(sg);
    //G4cout <<  "Generator is set" << G4endl;
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "FTFP") {
     G4TheoFSGenerator * model = new G4TheoFSGenerator;

     G4FTFModel *theStringModel = new G4FTFModel;
     G4ExcitedStringDecay *theStringDecay = new G4ExcitedStringDecay(new G4LundStringFragmentation);
     theStringModel->SetFragmentationModel(theStringDecay);

     G4GeneratorPrecompoundInterface *theCascade = new G4GeneratorPrecompoundInterface;
     G4PreCompoundModel *thePreEquilib = new G4PreCompoundModel(new G4ExcitationHandler);
     theCascade->SetDeExcitation(thePreEquilib);  

     model->SetTransport(theCascade);
     model->SetHighEnergyGenerator(theStringModel);
  
    sg = new Test43VSecondaryGenerator(model,mat);

    //G4cout <<  "Generator is ready" << G4endl;
    theProcess->SetSecondaryGenerator(sg);
    //G4cout <<  "Generator is set" << G4endl;
    man->AddDiscreteProcess(theProcess);

/*
  } else if(gen_name == "CHIPS") {

    sg = new Test43VSecondaryGenerator(new G4PionMinusNuclearReaction(),mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);
*/
  } else if(gen_name == "strCHIPS") {
    sg = new Test43VSecondaryGenerator(new G4StringChipsInterface(),mat);

    //G4cout <<  "Generator is ready" << G4endl;
    theProcess->SetSecondaryGenerator(sg);
    //G4cout <<  "Generator is set" << G4endl;
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "preCompound") {
    theDeExcitation = new G4ExcitationHandler();
    G4PreCompoundModel* pcm = new G4PreCompoundModel(theDeExcitation);
    thePreCompound = pcm;
    sg = new Test43VSecondaryGenerator(pcm,mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "binary_pc") {
    theDeExcitation = new G4ExcitationHandler();
    G4PreCompoundModel* pcm = new G4PreCompoundModel(theDeExcitation);
    thePreCompound = pcm;
    G4BinaryCascade* hkm = new G4BinaryCascade();
    sg = new Test43VSecondaryGenerator(hkm, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);
    hkm->SetDeExcitation(pcm);
    hkmod = hkm;

  } else if(gen_name == "binary") {
    G4BinaryCascade* hkm = new G4BinaryCascade();
    sg = new Test43VSecondaryGenerator(hkm, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);
//    hkm->SetDeExcitation(0);

  } else if(gen_name == "qgsbinary") {


    G4TheoFSGenerator * model = new G4TheoFSGenerator;
    G4QGSModel< G4QGSParticipants >  * stringmodel= new G4QGSModel< G4QGSParticipants >;
    G4BinaryCascade* cascade = new G4BinaryCascade();
    G4ExcitedStringDecay * stringDecay = new G4ExcitedStringDecay(new G4QGSMFragmentation);
    model->SetHighEnergyGenerator(stringmodel);
    stringmodel->SetFragmentationModel(stringDecay);
    model->SetTransport(cascade);
    sg = new Test43VSecondaryGenerator(model, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "ftfbinary") {


    G4TheoFSGenerator * model = new G4TheoFSGenerator;
    G4FTFModel * stringmodel= new G4FTFModel;
    G4BinaryCascade* cascade = new G4BinaryCascade();
    G4ExcitedStringDecay * stringDecay = new G4ExcitedStringDecay(new G4LundStringFragmentation());
    model->SetHighEnergyGenerator(stringmodel);
    stringmodel->SetFragmentationModel(stringDecay);
    model->SetTransport(cascade);
    sg = new Test43VSecondaryGenerator(model, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "binary_ion") {
    G4BinaryLightIonReaction* hkm = new G4BinaryLightIonReaction();
    sg = new Test43VSecondaryGenerator(hkm, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "bertini") {
    G4CascadeInterface* hkm = new G4CascadeInterface();
    sg = new Test43VSecondaryGenerator(hkm, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "LElastic") {
    G4LElastic* els = new G4LElastic();
    sg = new Test43VSecondaryGenerator(els, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else if(gen_name == "HElastic") {
    G4ElasticHadrNucleusHE* els = new G4ElasticHadrNucleusHE();
    sg = new Test43VSecondaryGenerator(els, mat);
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
    sg = new Test43VSecondaryGenerator(wam, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } else {
    G4cout << gen_name
           << " generator is unkown - no hadron production" << G4endl;
  }

  G4cout <<  "Secondary generator <"
         << gen_name << "> is initialized"
         << G4endl;
  return theProcess;

}	

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  






