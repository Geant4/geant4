#include "ProcessWrapper.hh"


// -> #include "G4TheoFSGenerator.hh"

ProcessWrapper::~ProcessWrapper()
{

   
   // need to see how desctructors of other data members work !!!...
   // ... to avoid potential memory leaks...
   //
   if ( fInteractionModel ) delete fInteractionModel;

}


G4double ProcessWrapper::PostStepGetPhysicalInteractionLength(const G4Track&,
                        	     		                    G4double,
			                                            G4ForceCondition* condition)
{

  *condition = NotForced;
  G4double z = DBL_MAX;

  return z;

}

G4VParticleChange* ProcessWrapper::PostStepDoIt( const G4Track& track, const G4Step& )
{

  if ( !fInteractionModel ) return 0;  
  
  if ( !fTargetNucleus ) return 0;
  
  
  // fNSec = 0;
  
  G4HadProjectile proj(track);
  G4HadFinalState* result = fInteractionModel->ApplyYourself( proj, *fTargetNucleus );
  result->SetTrafoToLab( proj.GetTrafoToLab() );
  
  ClearNumberOfInteractionLengthLeft();

  fPartChange.Initialize(track);

  G4int NSec = result->GetNumberOfSecondaries();
  G4int nb = NSec;
  if( result->GetStatusChange() == isAlive ) nb++;
  
  fPartChange.ProposeTrackStatus(fStopAndKill);
  fPartChange.SetNumberOfSecondaries(nb);

  for(G4int i=0; i<NSec; i++) {
    G4Track* tr = new G4Track(result->GetSecondary(i)->GetParticle(),
                              track.GetGlobalTime(),
	                      track.GetPosition());
    fPartChange.AddSecondary(tr);
  }

  if(result->GetStatusChange() == isAlive) {
    G4DynamicParticle* dp = new G4DynamicParticle(*(track.GetDynamicParticle()));
    G4Track* tr = new G4Track(dp,track.GetGlobalTime(),track.GetPosition());
    tr->SetKineticEnergy(result->GetEnergyChange());
    tr->SetMomentumDirection(result->GetMomentumChange());
    fPartChange.AddSecondary(tr);
  }

  result->Clear();

  return &fPartChange;

}

/*
QGSPWrapper::QGSPWrapper( const G4String& name = "QGSPWrapper" )
   : ProcessWrapper(name)
{

   fInteractionModel = G4TheoFSGenerator(); // common
   
   fStringModel      = new G4QGSModel< G4QGSParticipants >; // specific
   
   fCascade          = new G4GeneratorPrecompoundInterface(); // common
   
   fCascade->SetDeExcitation( new G4PreCompoundModel( new G4ExcitationHandler() ) ); // specific
   fStringDecay      = G4ExcitedStringDecay( new G4QGSMFragmentation() ); // specific
   fStringModel->SetFragmentationModel( fStringDecay ); // common interface but different input

   fInteractionModel->SetQuasiElasticChannel( new G4QuasiElasticChannel() ); // specific
   fInteractionModel->SetTransport( fCascade ); // common
   fInteractionModel->SetHighEnergyGenerator( fStringModel ); // common

   fInteractionModel->SetMinEnergy(GeV);  // common
   fInteractionModel->SetMaxEnergy(100.*TeV); 

}
*/
