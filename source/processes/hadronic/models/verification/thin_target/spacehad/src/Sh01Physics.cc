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
//      Sh01Physics 
//          
// 
//    Modified:
//    06.03.03 V. Grichine (based on test30)of V. Ivanchenko)
//
//

#include "G4CascadeInterface.hh"


#include "Sh01Physics.hh"

#include "G4UnitsTable.hh"
#include "Sh01SecondaryGenerator.hh"
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
// #include "G4PionMinusNuclearReaction.hh"
// #include "G4StringChipsInterface.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4BinaryCascade.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4CascadeInterface.hh"

////////////////////////////////////////////////////////////////////////

Sh01Physics::Sh01Physics()
{
  Initialise();
}

///////////////////////////////////////////////////////////////////////

Sh01Physics::~Sh01Physics()
{
  //  if(theProcess) delete theProcess;
}

//////////////////////////////////////////////////////////////////////////

void Sh01Physics::Initialise()
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

//////////////////////////////////////////////////////////////////////

G4VProcess* Sh01Physics::GetProcess( const G4String& gen_name, 
		                     const G4String& part_name,
		                            G4Material* mat)
{ 
  G4cout <<  "Sh01Physics entry" << G4endl;
  if(theProcess) delete theProcess;
  theProcess = 0;
	
  G4ProcessManager* man = 0;
	  
  if(part_name == "proton")
  {
    man = new G4ProcessManager(G4Proton::Proton());
  }
  else if(part_name == "pi+") 
  {
    man = new G4ProcessManager(G4PionPlus::PionPlus());
  }
  else if(part_name == "pi-") 
  {
    man = new G4ProcessManager(G4PionMinus::PionMinus());
  }
  else if(part_name == "neutron") 
  {
    man = new G4ProcessManager(G4Neutron::Neutron());
  }	  
  if(!man) return 0; 
	
  theProcess = new Sh01HadronProduction();	
  G4cout <<  "Process is created; gen= " << gen_name << G4endl;


  //////////////////////////////////////////////////////////////////////
  /// 	
  /// Physics list for the given run
  /// Choose generator

  Sh01SecondaryGenerator* sg = 0;
   
  if(gen_name == "LEparametrisation") 
  {
    if(part_name == "proton")   
    {
      sg = new Sh01SecondaryGenerator(new G4LEProtonInelastic(),mat);
    }
    else if(part_name == "pi+") 
    {
      sg = new Sh01SecondaryGenerator(new G4LEPionPlusInelastic(),mat);
    }
    else if(part_name == "pi-") 
    {
      sg = new Sh01SecondaryGenerator(new G4LEPionMinusInelastic(),mat);
    }
    else if(part_name == "neutron") 
    {
      sg = new Sh01SecondaryGenerator(new G4LENeutronInelastic(),mat);
    }
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);
  }
  /* 
  else if(gen_name == "CHIPS") 
  {
    sg = new Sh01SecondaryGenerator(new G4PionMinusNuclearReaction(),mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);

  } 
  else if(gen_name == "stringCHIPS") 
  {
    sg = new Sh01SecondaryGenerator(new G4StringChipsInterface(),mat);
		
    //G4cout <<  "Generator is ready" << G4endl;		
    theProcess->SetSecondaryGenerator(sg);
    //G4cout <<  "Generator is set" << G4endl;
    man->AddDiscreteProcess(theProcess);
  }
  */ 
  else if( gen_name == "preCompound" ) 
  {
    sg = new Sh01SecondaryGenerator(
         new G4PreCompoundModel(new G4ExcitationHandler()),mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);
		
  } 
  else if( gen_name == "kinetic" ) 
  {
    G4BinaryCascade*             hkm = new G4BinaryCascade();
    G4ExcitationHandler*      excite = new G4ExcitationHandler();
    G4VPreCompoundModel* precompound = new G4PreCompoundModel(excite);
    hkm->SetDeExcitation(precompound);

    sg = new Sh01SecondaryGenerator(hkm, mat);
    theProcess->SetSecondaryGenerator(sg);

    man->AddDiscreteProcess(theProcess);
  } 
  else if(gen_name == "binary") 
  {
    theDeExcitation = new G4ExcitationHandler();
    G4PreCompoundModel* pcm = new G4PreCompoundModel(theDeExcitation);
    thePreCompound = pcm;
    G4BinaryCascade* hkm = new G4BinaryCascade();
    sg = new Sh01SecondaryGenerator(hkm, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);
    hkm->SetDeExcitation(pcm);
    hkmod = hkm;
  } 
  else if(gen_name == "binary_no_pc") 
  {
    G4BinaryCascade* hkm = new G4BinaryCascade();
    sg = new Sh01SecondaryGenerator(hkm, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);
    hkm->SetDeExcitation(0);

  } 
  else if(gen_name == "binary_ion") 
  {
    G4BinaryLightIonReaction* hkm = new G4BinaryLightIonReaction();
    sg = new Sh01SecondaryGenerator(hkm, mat);
    theProcess->SetSecondaryGenerator(sg);
    man->AddDiscreteProcess(theProcess);
//    hkm->SetDeExcitation(0);

  }  
  else if( gen_name == "bertini" ) 
  {
    // G4CascadeInterface*             hkm = new G4CascadeInterface();

    G4CascadeInterface* hkm=new G4CascadeInterface();

    //G4ExcitationHandler*      excite = new G4ExcitationHandler();
    // G4VPreCompoundModel* precompound = new G4PreCompoundModel(excite);
    // hkm->SetDeExcitation(precompound);

    sg = new Sh01SecondaryGenerator(hkm, mat);
    theProcess->SetSecondaryGenerator(sg);

    man->AddDiscreteProcess(theProcess);
  } 
  else 
  {
    G4cout << gen_name 
           << " generator is unkown - no hadron production" << G4endl;
  }  
  G4cout <<  "Secondary generator <" 
         << gen_name << "> is initialized" 
         << G4endl;

  return theProcess;
}	

/////////////////////////////////////////////////////////////////////////
  






