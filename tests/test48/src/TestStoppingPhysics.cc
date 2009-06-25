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

#include "TestStoppingPhysics.hh"

#include "G4UnitsTable.hh"
//#include "Test30VSecondaryGenerator.hh"

#include "G4VRestProcess.hh"

#include "G4PhysicsTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4DecayPhysics.hh"

#include "G4MuonMinus.hh"
#include "G4PionMinus.hh"
#include "G4KaonMinus.hh"
#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"

/*
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4RPGPiPlusInelastic.hh"
#include "G4RPGPiMinusInelastic.hh"
#include "G4RPGProtonInelastic.hh"
#include "G4RPGNeutronInelastic.hh"

#include "G4StringChipsParticleLevelInterface.hh"
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

#include "G4ElasticHadrNucleusHE.hh"
#include "G4LElastic.hh"
#include "G4HadronElastic.hh"
#include "G4ChargeExchange.hh"
#include "G4CascadeElasticInterface.hh"
#include "G4QuasiElasticChannel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QStringChipsParticleLevelInterface.hh"

#include "HsQGSPInterface.hh"
#include "HsQGSCInterface.hh"
#include "HsQGSBInterface.hh"
#include "G4DiffuseElastic.hh"
*/

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"

#include "G4MuonMinusCaptureAtRest.hh"
#include "G4AntiProtonAnnihilationAtRest.hh"
#include "G4AntiNeutronAnnihilationAtRest.hh"
#include "G4PionMinusAbsorptionAtRest.hh"
#include "G4KaonMinusAbsorption.hh"

#include "G4QCaptureAtRest.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TestStoppingPhysics::TestStoppingPhysics()
{
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TestStoppingPhysics::~TestStoppingPhysics()
{
  //delete theDeExcitation;
  //delete thePreCompound;
  delete theProcess;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TestStoppingPhysics::Initialise()
{
  //G4DecayPhysics dp;
  //dp.ConstructParticle();
  theProcess = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VProcess* TestStoppingPhysics::GetProcess(const G4String& gen_name,
		                            const G4String& part_name,
		                            G4Material* mat)
{

  if(theProcess) delete theProcess;
  theProcess = 0;

  G4ProcessManager* man = 0;

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();  
    
  if (part_name == "anti_proton")   
  {
     // man = G4AntiProton::AntiProton()->GetProcessManager();
     man = new G4ProcessManager(G4AntiProton::AntiProton());
  }
  else if (part_name == "anti_neutron") 
  {
     // man = G4AntiNeutron::AntiNeutron()->GetProcessManager();
     man = new G4ProcessManager(G4AntiNeutron::AntiNeutron());
  }
  else if (part_name == "pi-") 
  {
     //man = G4PionMinus::PionMinus()->GetProcessManager();
     man = new G4ProcessManager(G4PionMinus::PionMinus());
  }
  else if (part_name == "kaon-")  
  {
     //man = G4KaonMinus::KaonMinus()->GetProcessManager();
     man = new G4ProcessManager(G4KaonMinus::KaonMinus());
  }

  if(!man) return 0;

  // Choose generator
  //
  if(gen_name == "stopping") 
  {
    if(part_name == "anti_proton")   
    {
      theProcess = new G4AntiProtonAnnihilationAtRest();
    }
    else if (part_name == "anti_neutron" )
    {
       theProcess = new G4AntiNeutronAnnihilationAtRest();
    }
    else if(part_name == "pi-") 
    {  
      theProcess = new G4PionMinusAbsorptionAtRest();
    }
    else if ( part_name == "kaon-")
    {
       theProcess = new G4KaonMinusAbsorption();
    }
  } 
  else if(gen_name == "CHIPS") 
  {
    theProcess = new G4QCaptureAtRest();
  } 
  else 
  {
    G4cout << gen_name
           << " generator is unkown - no hadron production" << G4endl;
  }
  
  man->AddRestProcess(theProcess);

  G4cout <<  " Model <"
         << gen_name << "> is initialized"
         << G4endl;
  
  return theProcess;

}	

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  






