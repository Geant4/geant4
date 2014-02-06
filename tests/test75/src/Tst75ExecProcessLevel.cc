
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


