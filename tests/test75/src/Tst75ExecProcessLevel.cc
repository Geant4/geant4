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

#include "Tst75ExecProcessLevel.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4ios.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleChange.hh"

#include "TstDiscreteProcessReader.hh"

//#include "G4Gamma.hh"
//#include "G4ForceCondition.hh"

#include "G4CascadeInterface.hh"
#include "G4PhotoNuclearProcess.hh" // hadronic process
// #include "G4GammaNuclearReaction.hh" // hadronic interaction (aka "model")


Tst75ExecProcessLevel::Tst75ExecProcessLevel( const TstReader* pset )
   : ExecProcessLevel(pset)
{

   InitProcess( pset );
   InitSetup(   pset );
   InitBeam(    pset );
   
   G4ParticleDefinition* partDef = ( G4ParticleTable::GetParticleTable())->FindParticle(pset->GetBeamParticle() );
   G4DynamicParticle dParticle( partDef, pset->GetDirection(), GetBeam()->GetBeamEnergy() );
   fXSecOnTarget = fProcWrapper->GetElementCrossSection( &dParticle, GetTarget()->GetCurrentMaterial()->GetElement(0) );

}

void Tst75ExecProcessLevel::InitProcess( const TstReader* pset )
{

   // G4String name = (dynamic_cast<const TstDiscreteProcessReader*>(pset))->GetProcessName();
   G4String NameModel = pset->GetPhysics();
   G4String NamePart  = pset->GetBeamParticle();

   // ProcessWrapper* proc = 0;
   G4HadronicProcess* proc = 0;
   // G4VDiscreteProcess* proc = 0;
   G4HadronicInteraction* model = 0;
    
   if ( NamePart =="gamma" )
   {
       proc = new G4PhotoNuclearProcess();
   }
    
   if ( NameModel == "Bertini" )
   {
       model = new G4CascadeInterface();
   }
//    else if ( NameModel == "CHIPS" )
//    {
//       model = new G4GammaNuclearReaction();
//    }
    
   // make sure the pointers are valid
   //
   assert(proc);
   assert(model);
    
   // configure tthe model 
   // (to a specific energy range)
   //
   G4double Mom = pset->GetBeamMomentum();
   model->SetMinEnergy(0.9*Mom);
   model->SetMaxEnergy(1.1*Mom);
    
   // now assign model to the process
   //
   proc->RegisterMe(model);
   
   fProcWrapper = proc;
    
   if (!fProcWrapper) 
   { 
      G4cout 
	     << " generator " << NameModel << " is unavailable"
	     << G4endl;
      exit(1);
   } 
    
   std::cout << " process = " << fProcWrapper->GetProcessName() << std::endl;
        
   // fProcWrapper->Compose();
   
   return;

}


