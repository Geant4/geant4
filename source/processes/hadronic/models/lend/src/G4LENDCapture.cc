
#include "G4LENDCapture.hh"

#include "G4Nucleus.hh"
#include "G4ParticleTable.hh"
  
G4HadFinalState * G4LENDCapture::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTarg )
{

   G4double temp = aTrack.GetMaterial()->GetTemperature();

   //G4int iZ = int ( aTarg.GetZ() );
   //G4int iA = int ( aTarg.GetN() );
   //migrate to integer A and Z (GetN_asInt returns number of neutrons in the nucleus since this) 
   G4int iZ = aTarg.GetZ_asInt();
   G4int iA = aTarg.GetA_asInt();

   G4double ke = aTrack.GetKineticEnergy();

   G4HadFinalState* theResult = &theParticleChange;
   theResult->Clear();

   G4GIDI_target* aTarget = usedTarget_map.find( lend_manager->GetNucleusEncoding( iZ , iA ) )->second->GetTarget();
   std::vector<G4GIDI_Product>* products = aTarget->getCaptureFinalState( ke*MeV, temp, NULL, NULL );


   if ( products != NULL ) 
   {

      G4ThreeVector p(0,0,0);
      G4int totN = 0;

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

         G4ThreeVector dp((*products)[j].px,(*products)[j].py,(*products)[j].pz);
         p += dp;
          
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

         theSec->SetMomentum( G4ThreeVector( (*products)[j].px*MeV , (*products)[j].py*MeV , (*products)[j].pz*MeV ) );

         if ( dp.mag() == 0 ) 
         {
            theSec->SetMomentum( -p*MeV ); 
         }

         theResult->AddSecondary( theSec );
      } 
   }
   delete products;

   theResult->SetStatusChange( stopAndKill );

   return theResult; 

}
