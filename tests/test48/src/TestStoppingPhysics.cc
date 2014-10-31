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
//
//      History: based on object model of
//      ---------- TestStoppingPhysics -------
//                   by Julia Yarba 
// 
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "TestStoppingPhysics.hh"

#include "G4VRestProcess.hh"

#include "G4PhysicsTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"

#include "G4MuonMinus.hh"
#include "G4PionMinus.hh"
#include "G4KaonMinus.hh"
#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"
#include "G4SigmaMinus.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4BosonConstructor.hh"

#include "G4SystemOfUnits.hh"

//
// (what remains of the) traditional stopping code
//
//#include "G4MuonMinusCaptureAtRest.hh"
#include "G4AntiNeutronAnnihilationAtRest.hh"

//
// new development, in place since g4.9.6.b01
//
#include "G4MuonMinusCapture.hh"

#include "G4PreCompoundModel.hh"

// Bertini cascade
// local code for pi-, K- & Sigma-
// (very original implementation before the standard interface
//  went in place; still kept here for cross-checks)
#include "TestBertiniStopping.hh"
//
// standard interface for pi-
//
#include "G4PiMinusAbsorptionBertini.hh"
//
// standard interface for K- and Sigma- available since g4.9.6.b01
//
#include "G4KaonMinusAbsorptionBertini.hh"
#include "G4SigmaMinusAbsorptionBertini.hh"
 
#if defined (USE_MUCAPTURE)
// new Muon code
#include "G4MuonMinusCapturePhysics.hh"
#endif

// FTF model for pbar and AntiSigma-
// old code up to 4.9.6.b01
//#include "G4FTFCaptureAtRest.hh"
// new code starting 4.9.5-ref09
#include "G4AntiProtonAbsorptionFritiof.hh"
#include "G4AntiSigmaPlusAbsorptionFritiof.hh"

//#include "G4SpecialMuMinusCapturePrecompound.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TestStoppingPhysics::TestStoppingPhysics(G4int verbose):
#if defined (USE_MUCAPTURE)
  theMuonMinusCaptureConstructor(new G4MuonMinusCapturePhysics(verbose)),
#endif
  theProcess(0), theProcessMan(0), verboseLevel(verbose)
{

  //G4DecayPhysics pDecayPhysics;
  //pDecayPhysics.ConstructParticle(); // it is an example; we do not do it here though
  
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();  
  
  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();
  
  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();
  
  G4BosonConstructor pBosonConstructor;
  pBosonConstructor.ConstructParticle();
  
#if defined (USE_MUCAPTURE)
  // for the new mu capture
  theMuonMinusCaptureConstructor->ConstructParticle();
#endif
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TestStoppingPhysics::~TestStoppingPhysics()
{
  //delete theDeExcitation;
  //delete thePreCompound;
  if ( theProcess )                     delete theProcess;
  if ( theProcessMan)                   delete theProcessMan;
#if defined (USE_MUCAPTURE)
  if ( theMuonMinusCaptureConstructor ) delete theMuonMinusCaptureConstructor;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VProcess* TestStoppingPhysics::GetProcess(const G4String& gen_name,
		                            const G4String& part_name)
{

  G4cout << part_name 
         << gen_name 
         << " initializing 0"
         << G4endl;         

  if(theProcess) delete theProcess;
  theProcess = 0;
  if ( theProcessMan) delete theProcessMan;
  theProcessMan = 0;

//  G4ProcessManager* man = 0;
    
  if (part_name == "anti_proton")   
  {
//     man = new G4ProcessManager(G4AntiProton::AntiProton());
     theProcessMan = new G4ProcessManager(G4AntiProton::AntiProton());
  }
  else if (part_name == "anti_neutron") 
  {
//     man = new G4ProcessManager(G4AntiNeutron::AntiNeutron());
     theProcessMan = new G4ProcessManager(G4AntiNeutron::AntiNeutron());
  }
  else if (part_name == "pi-") 
  {
//     man = new G4ProcessManager(G4PionMinus::PionMinus());
     theProcessMan = new G4ProcessManager(G4PionMinus::PionMinus());
  }
  else if (part_name == "kaon-")  
  {
//     man = new G4ProcessManager(G4KaonMinus::KaonMinus());
     theProcessMan = new G4ProcessManager(G4KaonMinus::KaonMinus());
  }
  else if (part_name == "mu-")  
  {

#if defined (USE_MUCAPTURE)

    // in case of new mucapture muons we need more particles, processes and Process Managers
    // creating process managers if not existing
    G4ParticleTable*   partTable = G4ParticleTable::GetParticleTable();
    G4ParticleTable::G4PTblDicIterator* theParticleIterator = partTable->GetIterator();

    theParticleIterator->reset();

    while( (*theParticleIterator)() ) {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();

      if (verboseLevel>1){
        G4cout <<"TestStoppingPhysics::GetProcess"
               <<" : Checking Process Manager for "
               << particle->GetParticleName() 
               << G4endl;
      }

      if ( pmanager == 0) {
        // no process manager
        if (verboseLevel>0){
          G4cout <<"TestStoppingPhysics::GetProcess"
                 <<" : Creating Process Manager for "
                 << particle->GetParticleName() 
                 << G4endl;
        }
        pmanager = new G4ProcessManager(particle); // FIXME leak?
                                                   // something needs
                                                   // to clean it up
        particle->SetProcessManager(pmanager);
      }
    }

    // does the above cover objects in  G4MuAtomTable ??? does it need to ???
    theMuonMinusCaptureConstructor->ConstructProcess();
//    man = G4MuonMinus::MuonMinus()->GetProcessManager();
    theProcessMan = G4MuonMinus::MuonMinus()->GetProcessManager();

#else

//    man = new G4ProcessManager(G4MuonMinus::MuonMinus());
    theProcessMan = new G4ProcessManager(G4MuonMinus::MuonMinus());
	 
    if (verboseLevel>1){   // use verboseLevel to suppress compilation warning
	    G4cout <<"TestStoppingPhysics::GetProcess() not using new mucapture"<<G4endl;
    }		
#endif

  }
  else if ( part_name == "sigma-" )
  {
//     man = new G4ProcessManager(G4SigmaMinus::SigmaMinus());
     theProcessMan = new G4ProcessManager(G4SigmaMinus::SigmaMinus());
  }

//  if(!man) return 0;
  if(!theProcessMan) return 0;

  G4cout << part_name 
         << gen_name 
         << " initializing 1"
         << G4endl;         


  // Choose generator
  //
  if(gen_name == "stopping" || gen_name == "captureUpdate"
#if defined (USE_MUCAPTURE)
     || gen_name == "newmustopping" 
#endif
     ) 
  {

    if (part_name == "anti_neutron" )
    {
       theProcess = new G4AntiNeutronAnnihilationAtRest();
    }
    else if ( part_name == "mu-")
    {

      if (gen_name == "stopping" )
      {
	//   theProcess = new G4MuonMinusCaptureAtRest();
          theProcess = new G4MuonMinusCapture();
      }
      else if ( gen_name == "captureUpdate" )
      {
          theProcess = new G4MuonMinusCapture();
          // theProcess = new G4MuonMinusCapture(new G4SpecialMuMinusCapturePrecompound());
      }
#if defined (USE_MUCAPTURE) 
      else if (gen_name == "newmustopping" )
        {
          //          theMuonMinusCaptureConstructor->ConstructProcess();
          // constructed earlier in the case of mu- (FIXME should have been guarded by an if ???)
          theProcess = theMuonMinusCaptureConstructor->GetMuonMinusAtomicCaptureProcess();
          // new MuonStopping is done differently and we need to exit without "re adding" the process
          G4cout <<  " Model <"
                 << gen_name << "> is initialized"
                 << G4endl;         
          return theProcess;
        }
#endif
    }
  } 
  else if ( gen_name == "Bertini" || gen_name == "BertiniPreCo" )
  {

     if ( part_name == "kaon-" ) 
     {
        // for K-, still use the local interface
	// although it should be noted that PreCo 
	// doesn't do anything for K- on H (the only one we have data for)
	// since there're no excited nuclei
/*
	TestBertiniStopping* proc = new TestBertiniStopping();
        if ( gen_name == "BertiniPreCo" )
        {
           proc->UsePreCompound();
        }
        proc->InitTarget(mat);
        theProcess = proc;
*/
        // NOTE (JVY):
	// Standard interface will pull in an "EMcascade" part 
	// of the capture business (unlike "Test" interface); 
	// secondaries coming from the "EMcascade" will be
	// filtered out for K- by the analysis part of the code, 
	// but Sigma- they're still in the histo... 
	//
	theProcess = new G4KaonMinusAbsorptionBertini();
     }
     else if ( part_name == "sigma-" )
     {
        theProcess = new G4SigmaMinusAbsorptionBertini();
     }
     else if ( part_name == "pi-" )
     {
        // for pi-, standard interface is available starting G4.9.5 
	// PreCo activation is hardcoded feature in it, 
	// i.e. no switch to turn it off - apparently, based on earlier validation results
	// that indicated great improvement if Bertini_PreCo used
	theProcess = new G4PiMinusAbsorptionBertini();
/*
        TestBertiniStopping* proc = new TestBertiniStopping();
        if ( gen_name == "BertiniPreCo" )
        {
           proc->UsePreCompound();
        }
        proc->InitTarget(mat);
        theProcess = proc;
*/
     }

  }
  else if ( gen_name == "FTF")
  {
     if ( part_name == "anti_proton" )
     {
        theProcess = new G4AntiProtonAbsorptionFritiof();
     }
  }
  else 
  {
    G4cout << gen_name
           << " generator is unkown - no hadron production" << G4endl;
  }
  
//  man->AddRestProcess(theProcess);
  theProcessMan->AddRestProcess(theProcess);

  G4cout <<  " Model <"
         << gen_name << "> is initialized"
         << G4endl;
  
  return theProcess;

}	

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
