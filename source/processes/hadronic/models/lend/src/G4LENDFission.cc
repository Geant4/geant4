#include "G4LENDFission.hh"

#include "G4Nucleus.hh"
#include "G4ParticleTable.hh"
  
G4HadFinalState * G4LENDFission::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTarg )
{

   G4double temp = aTrack.GetMaterial()->GetTemperature();

   G4int iZ = int ( aTarg.GetZ() );
   G4int iA = int ( aTarg.GetN() );

   G4double ke = aTrack.GetKineticEnergy();

   G4HadFinalState* theResult = &theParticleChange;
   theResult->Clear();

   GIDI4GEANT_target* aTarget = usedTarget_map.find( endl_manager->GetNucleusEncoding( iZ , iA ) )->second->GetTarget();
   vector<GIDI4GEANT_Product>* products = aTarget->getFissionFinalState( ke*MeV, temp, NULL, NULL );
   if ( products != NULL ) 
   {
      for ( G4int j = 0; j < int( products->size() ); j++ ) 
      {
         G4int iZ = (*products)[j].Z; 
         G4int iA = (*products)[j].A; 

         //G4cout << "Z = "    << (*products)[j].Z 
         //       << ", A = "  << (*products)[j].A 
         //       << ", EK = " << (*products)[j].kineticEnergy << " [MeV]" 
         //       << ", px = " << (*products)[j].px
         //       << ", py = " << (*products)[j].py
         //       << ", pz = " << (*products)[j].pz
         //       << G4endl;

         G4DynamicParticle* theSec = new G4DynamicParticle;

         if ( iZ > 0 )
         {
                                                                                        //Ex  j?
            theSec->SetDefinition( G4ParticleTable::GetParticleTable()->FindIon( iZ, iA , 0, 0 ) );
         } 
         else if ( iA == 1 && iZ == 0 )
         {
            theSec->SetDefinition( G4Neutron::Neutron() );
         } 
         else
         {
            theSec->SetDefinition( G4Gamma::Gamma() );
         } 

         theSec->SetMomentum( G4ThreeVector( (*products)[j].px*MeV , (*products)[j].py*MeV , (*products)[j].pz*MeV ) );
         //G4cout << theSec->GetDefinition()->GetParticleName() << G4endl;
         theResult->AddSecondary( theSec );
      } 
   }
   delete products;

   theResult->SetStatusChange( stopAndKill );

   return theResult; 

}
