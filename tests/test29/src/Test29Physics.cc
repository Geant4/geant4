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

#include "Test29Physics.hh"

#include "G4UnitsTable.hh"
#include "Test29VSecondaryGenerator.hh"
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
#include "G4QCaptureAtRest.hh"
//#include "G4PionMinusNuclearReaction.hh"
//#include "G4StringChipsInterface.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4BinaryCascade.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4CascadeInterface.hh"

Test29Physics::Test29Physics()
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

  // Nucleons
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();

  // Nuclear Fragments (Ions)
  G4GenericIon::GenericIonDefinition();
  G4Deuteron::DeuteronDefinition();
  G4Alpha::AlphaDefinition();
  G4Triton::TritonDefinition();

  theWDProcess = 0;
  theWRProcess = 0;
  theHDProcess = 0;
  theDeExcitation = 0;
  thePreCompound = 0;
  hkmod = 0;
}

Test29Physics::~Test29Physics()
{
  //  if(theHDProcess) delete theHDProcess;
}

// Definition of the Progectile and the Process for the interaction (onFlight or atRest)
G4VProcess* Test29Physics::GetProcess(const G4String& gen_name, const G4String& part_name,
		                              G4Material* mat, G4bool atrest, G4bool elweak)
{
  G4cout<<"Test29Physics::GetProcess's called Gen="<<gen_name<<",Proj="<<part_name<<G4endl;
  if(theWDProcess) delete theWDProcess;
  if(theWRProcess) delete theWRProcess;
  if(theHDProcess) delete theHDProcess;
  theWDProcess = 0;
  theWRProcess = 0;
  theHDProcess = 0;

  G4ProcessManager* man = 0;

  // Initial particle definition (add more if it is necessary, but initialize them above)
  if     (part_name == "proton")  man = new G4ProcessManager(G4Proton::Proton());
  else if(part_name == "neutron") man = new G4ProcessManager(G4Neutron::Neutron());
  else if(part_name == "antiP")   man = new G4ProcessManager(G4AntiProton::AntiProton());
  else if(part_name == "antiN")   man = new G4ProcessManager(G4AntiNeutron::AntiNeutron());
  else if(part_name == "pi+")     man = new G4ProcessManager(G4PionPlus::PionPlus());
  else if(part_name == "pi-")     man = new G4ProcessManager(G4PionMinus::PionMinus());
  else if(part_name == "K+")      man = new G4ProcessManager(G4KaonPlus::KaonPlus());
  else if(part_name == "K-")      man = new G4ProcessManager(G4KaonMinus::KaonMinus());
  else if(part_name == "Sigma-")  man = new G4ProcessManager(G4SigmaMinus::SigmaMinus());
  else if(part_name == "Xi-")     man = new G4ProcessManager(G4XiMinus::XiMinus());
  else if(part_name == "Omega-")  man = new G4ProcessManager(G4OmegaMinus::OmegaMinus());
  else if(part_name == "aS+") man = new G4ProcessManager(G4AntiSigmaPlus::AntiSigmaPlus());
  else if(part_name == "deuteron")man = new G4ProcessManager(G4Deuteron::Deuteron());
  else if(part_name == "alpha")   man = new G4ProcessManager(G4Alpha::Alpha());

  if(!man)
  {
    G4cout<<"Test29Physics::GetProcess: Projectile="<<part_name<<" is not found"<<G4endl;
    return 0;
  }

  if(elweak) // The ElectroWeek Process is selected for this RUN
  {
    if(atrest)
    {
      theWRProcess = new G4QCaptureAtRest; // Universal CHIPS ElectroWeak G4RestProcess
      man->AddRestProcess(theWRProcess);
      return theWRProcess;
    }
    else
    {
      G4cout<<"Test29Physics::GetProcess: ElectroWeakDiscreteProcess isn't made"<<G4endl;
      exit(1);
      //G4QInteractionOnFlight theWDProcess;// UniversalCHIPS ElectroWeak G4DiscreteProcess
      //man->AddDiscreteProcess(theWDProcess);
      //return theWDProcess;
    }
  }
  else       // The Hadronic Process is selected for this RUN
  {
    if(atrest)
    {
      G4cout<<"Test29Physics::GetProcess:HadronicAtRestProcess is not implemented"<<G4endl;
      exit(1);
    }
    else
    {
      theHDProcess = new Test29HadronProduction();
      G4cout<<"HadronicProcessPrototype is created gen="<<gen_name<<G4endl;

      // Definition of the Model from the Hadronic Collection
      Test29VSecondaryGenerator* sg = 0; //@@ Make similar for the capture at rest (?)

      // Choose the generator from the HadronicModelCollection (add more if necessary)
      if(gen_name == "LEparametrisation")
      {
        if     (part_name == "proton")
          sg = new Test29VSecondaryGenerator(new G4LEProtonInelastic(),mat);
        else if(part_name == "pi+")
          sg = new Test29VSecondaryGenerator(new G4LEPionPlusInelastic(),mat);
        else if(part_name == "pi-")
          sg = new Test29VSecondaryGenerator(new G4LEPionMinusInelastic(),mat);
        else if(part_name == "neutron")
          sg = new Test29VSecondaryGenerator(new G4LENeutronInelastic(),mat);
        theHDProcess->SetSecondaryGenerator(sg);
        man->AddDiscreteProcess(theHDProcess);
      }
      //else if(gen_name == "CHIPS")
      //{
      //  sg = new Test29QSecondaryGenerator(new G4QCaptureAtRest(),mat);
      //  theHDProcess->SetSecondaryGenerator(sg);
      //  man->AddDiscreteProcess(theHDProcess);
      //}
      else if(gen_name == "stringCHIPS")
      {
        sg = new Test29VSecondaryGenerator(new G4StringChipsInterface(),mat);
        //G4cout <<  "Generator is ready" << G4endl;
        theHDProcess->SetSecondaryGenerator(sg);
        //G4cout <<  "Generator is set" << G4endl;
        man->AddDiscreteProcess(theHDProcess);
      }
      else if(gen_name == "preCompound")
      {
        theDeExcitation = new G4ExcitationHandler();
        G4PreCompoundModel* pcm = new G4PreCompoundModel(theDeExcitation);
        thePreCompound = pcm;
        sg = new Test29VSecondaryGenerator(pcm,mat);
        theHDProcess->SetSecondaryGenerator(sg);
        man->AddDiscreteProcess(theHDProcess);
      }
      else if(gen_name == "binary")
      {
        theDeExcitation = new G4ExcitationHandler();
        G4PreCompoundModel* pcm = new G4PreCompoundModel(theDeExcitation);
        thePreCompound = pcm;
        G4BinaryCascade* hkm = new G4BinaryCascade();
        sg = new Test29VSecondaryGenerator(hkm, mat);
        theHDProcess->SetSecondaryGenerator(sg);
        man->AddDiscreteProcess(theHDProcess);
        hkm->SetDeExcitation(pcm);
        hkmod = hkm;
      }
      else if(gen_name == "binary_no_pc")
      {
        G4BinaryCascade* hkm = new G4BinaryCascade();
        sg = new Test29VSecondaryGenerator(hkm, mat);
        theHDProcess->SetSecondaryGenerator(sg);
        man->AddDiscreteProcess(theHDProcess);
        hkm->SetDeExcitation(0);
      }
      else if(gen_name == "binary_ion")
      {
        G4BinaryLightIonReaction* hkm = new G4BinaryLightIonReaction();
        sg = new Test29VSecondaryGenerator(hkm, mat);
        theHDProcess->SetSecondaryGenerator(sg);
        man->AddDiscreteProcess(theHDProcess);
        //hkm->SetDeExcitation(0);
      }
      else if(gen_name == "bertini")
      {
        G4CascadeInterface* hkm = new G4CascadeInterface();
        sg = new Test29VSecondaryGenerator(hkm, mat);
        theHDProcess->SetSecondaryGenerator(sg);
        man->AddDiscreteProcess(theHDProcess);
      }
      else G4cout<<"Test29Physics::GetProcess: "<<gen_name<<" generator is unkown"<<G4endl;
    }
    G4cout<<"Secondary generator "<<gen_name<<" is initialized"<< G4endl;
    return theHDProcess;
  }
  return 0;
}	
  






