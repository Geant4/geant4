
#include "TstHistoSet.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"

#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"
#include "G4DecayProducts.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

void TstHistoSet::AccountForResDecay( G4VParticleChange* pchg )
{

   G4int NSec = pchg->GetNumberOfSecondaries();
        
   const G4Track*           trk = 0;
   const G4DynamicParticle* sec = 0;
       
   for (G4int i=0; i<NSec; i++) 
   {
      trk = pchg->GetSecondary(i);
      sec = trk->GetDynamicParticle();			
//      const G4String& pname = sec->GetDefinition()->GetParticleName();
      // G4 doesn't seem to produce an omega (783) meson (at all)
      // it's unclear if it even know such resonance
      // it only seems to know Omega baryons...
//       if ( pname == "eta" || pname == "eta_prime" || pname == "rho0"    ) 
//            || pname == "kaon0S" )
//           || pname == "kaon0L" )
      // if ( sec->GetDefinition()->IsShortLived() ) // somehow this flaf is ALWAYS FASLE !!!
      //
      // OK, I think the reasonable compromise would be to decay everything with life<0.01*ns
      // provided that it has a decay table !!! (because B's, D's, and several others do NOT)
      // (for example, K0S life~0.89ns, and Lambda life~0.026ns - those have been removed, 
      //  or at least attempted to remove from the NA49 samples)
      //
      if ( sec->GetDefinition()->GetPDGLifeTime() < 0.01*ns && sec->GetDefinition()->GetDecayTable() != NULL )
      {	 
	 G4DecayTable*    decTable    = sec->GetDefinition()->GetDecayTable();
         G4VDecayChannel* decChannel  = decTable->SelectADecayChannel();
         G4DecayProducts* decProducts = decChannel->DecayIt( sec->GetMass() );
         decProducts->Boost( sec->GetTotalEnergy(), sec->GetMomentumDirection() );
	 //
	 // there could also be (G4Decay::)DaughterPolarization( G4Track*, G4DecayProducts* )
	 // but in G4Decay this method is a) protected and b) empty !!!
	 //	       
         G4int ndec = decProducts->entries();
	 for ( G4int id=0; id<ndec; ++id )
	 {
	    G4DynamicParticle* dsec = (*decProducts)[id]; 
	    G4double time = trk->GetGlobalTime();
	    const G4ThreeVector& pos = trk->GetPosition();
	    pchg->AddSecondary( new G4Track( dsec, time, pos ) );
         }
      }        
   }   

   return;

}

/*
void TstHistoSet::AccountForPionDecay( const G4Track* trk )
{
   
   fPionDecay = 0;
   
   const G4DynamicParticle* sec = trk->GetDynamicParticle();
   
   const G4String& pname = sec->GetDefinition()->GetParticleName();
   if ( pname != "pi-" && pname != "pi+" ) return;
   
   G4DecayTable* decTable = sec->GetDefinition()->GetDecayTable();
   if ( decTable == NULL ) return;
   
   G4VDecayChannel* decChannel  = decTable->SelectADecayChannel();
   fPionDecay = decChannel->DecayIt( sec->GetMass() );
   fPionDecay->Boost( sec->GetTotalEnergy(), sec->GetMomentumDirection() );
   
   G4int ndec = fPionDecay->entries();
   for ( G4int id=0; id<ndec; ++id )
   {
      G4DynamicParticle* dsec = (*fPionDecay)[id]; 
      const G4String& dname = dsec->GetDefinition()->GetParticleName();
      if ( dname == "nu_mu" ) // pi+ --> mu+ nu_mu
      {
	 double nu_ekin = dsec->GetKineticEnergy() / GeV ;
	 if ( nu_ekin > 0. && nu_ekin <= 2. )
	 {
	    fNuERange = kR_0_2;
	 }
	 else if ( nu_ekin > 2. && nu_ekin <= 5. )
	 {
	    fNuERange = kR_2_5;
	 }
	 else if ( nu_ekin > 5. && nu_ekin <= 10. )
	 {
	    fNuERange = kR_5_10;
	 }
	 else if ( nu_ekin > 10. && nu_ekin <= 20. )
	 {
	    fNuERange = kR_10_20;
	 }
	 else if ( nu_ekin > 20. && nu_ekin <= 50. )
	 {
	    fNuERange = kR_20_50;
	 } 
	 break;
      }
      if ( dname == "anti_nu_mu" ) // pi- --> mu- anti_nu_mu
      {
	 double nubar_ekin = dsec->GetKineticEnergy() / GeV ;
	 if ( nubar_ekin > 0. && nubar_ekin <= 2. )
	 {
	    fNubarERange = kR_0_2;
	 }
	 else if ( nubar_ekin > 2. && nubar_ekin <= 5. )
	 {
	    fNubarERange = kR_2_5;
	 }
	 else if ( nubar_ekin > 5. && nubar_ekin <= 10. )
	 {
	    fNubarERange = kR_5_10;
	 }
	 else if ( nubar_ekin > 10. && nubar_ekin <= 20. )
	 {
	    fNubarERange = kR_10_20;
	 }
	 else if ( nubar_ekin > 20. && nubar_ekin <= 50. )
	 {
	    fNubarERange = kR_20_50;
	 } 
	 break;
      }
   }
         
   return;

}
*/
