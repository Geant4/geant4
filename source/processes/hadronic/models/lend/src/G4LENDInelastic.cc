#include "G4LENDInelastic.hh"

#include "G4Nucleus.hh"
#include "G4ParticleTable.hh"
  
G4HadFinalState * G4LENDInelastic::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTarg )
{

   G4ThreeVector proj_p = aTrack.Get4Momentum().vect();

   G4double temp = aTrack.GetMaterial()->GetTemperature();

   G4int iZ = int ( aTarg.GetZ() );
   G4int iA = int ( aTarg.GetN() );
   //G4cout << "target: Z = " << iZ << " N = " << iA << G4endl;

   G4double ke = aTrack.GetKineticEnergy();
   //G4cout << "projectile: KE = " << ke/MeV << " [MeV]" << G4endl;

   G4HadFinalState* theResult = &theParticleChange;
   theResult->Clear();

   G4GIDI_target* aTarget = usedTarget_map.find( endl_manager->GetNucleusEncoding( iZ , iA ) )->second->GetTarget();
   std::vector<G4GIDI_Product>* products = aTarget->getOthersFinalState( ke*MeV, temp, NULL, NULL );
   if ( products != NULL ) 
   {

      G4ThreeVector psum(0);

      int totN = 0;
      for ( G4int j = 0; j < int( products->size() ); j++ ) 
      {

         G4int jZ = (*products)[j].Z; 
         G4int jA = (*products)[j].A; 

         //G4cout << "ZA = " << 1000 * (*products)[j].Z + (*products)[j].A << "  EK = "
         //     << (*products)[j].kineticEnergy
         //     << " px  " <<  (*products)[j].px
         //     << " py  " <<  (*products)[j].py
         //     << " pz  " <<  (*products)[j].pz
         //     << G4endl;

         G4DynamicParticle* theSec = new G4DynamicParticle;

         if ( jA == 1 && jZ == 1 )
         {
            theSec->SetDefinition( G4Proton::Proton() );
            totN += 1;
         }
         else if ( jA == 1 && jZ == 0 )
         {
            theSec->SetDefinition( G4Neutron::Neutron() );
            totN += 1;
         } 
         else if ( jZ > 0 )
         {
            if ( jA != 0 )
            {
               theSec->SetDefinition( G4ParticleTable::GetParticleTable()->FindIon( jZ , jA , 0 , 0 ) );
               totN += jA;
            }
            else 
            {
               theSec->SetDefinition( G4ParticleTable::GetParticleTable()->FindIon( jZ , iA+1-totN , 0 , 0 ) );
            }
         } 
         else
         {
            theSec->SetDefinition( G4Gamma::Gamma() );
         } 

         G4ThreeVector p( (*products)[j].px*MeV , (*products)[j].py*MeV , (*products)[j].pz*MeV ); 
         psum += p; 
         if ( p.mag() == 0 ) p = proj_p - psum;

         theSec->SetMomentum( p );

         theResult->AddSecondary( theSec );
      } 
   }
   delete products;

   theResult->SetStatusChange( stopAndKill );

   return theResult; 

}
