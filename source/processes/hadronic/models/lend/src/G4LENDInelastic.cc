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
#include "G4LENDInelastic.hh"

#include "G4SystemOfUnits.hh"
#include "G4Nucleus.hh"
#include "G4IonTable.hh"
  
G4HadFinalState * G4LENDInelastic::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTarg )
{
   //return preco->ApplyYourself( aTrack, aTarg );

   G4ThreeVector proj_p = aTrack.Get4Momentum().vect();

   G4double temp = aTrack.GetMaterial()->GetTemperature();

   //G4int iZ = int ( aTarg.GetZ() );
   //G4int iA = int ( aTarg.GetN() );
   //migrate to integer A and Z (GetN_asInt returns number of neutrons in the nucleus since this) 
   G4int iZ = aTarg.GetZ_asInt();
   G4int iA = aTarg.GetA_asInt();
   //G4int iM = aTarg.GetM_asInt();
   G4int iM = 0;
   if ( aTarg.GetIsotope() != NULL ) {
      iM = aTarg.GetIsotope()->Getm();
   }
   //G4cout << "target: Z = " << iZ << " N = " << iA << G4endl;

   G4double ke = aTrack.GetKineticEnergy();
   //G4cout << "projectile: KE = " << ke/MeV << " [MeV]" << G4endl;

   G4HadFinalState* theResult = &theParticleChange;
   theResult->Clear();

   G4GIDI_target* aTarget = get_target_from_map( lend_manager->GetNucleusEncoding( iZ , iA , iM ) );
   if ( aTarget == NULL ) return returnUnchanged( aTrack , theResult );

   std::vector<G4GIDI_Product>* products; 
   for ( G4int i = 0 ; i != 1024 ; i++ ) {
      products = aTarget->getOthersFinalState( ke*MeV, temp, MyRNG, NULL );
      if ( products != NULL ) break;
   }
   //return preco->ApplyYourself( aTrack, aTarg );

   G4int iTotZ = iZ + aTrack.GetDefinition()->GetAtomicNumber();
   G4int iTotA = iA + aTrack.GetDefinition()->GetAtomicMass();

   if ( products != NULL ) 
   {
      //G4cout << "Using LENDModel" << G4endl;

      G4ThreeVector psum(0);
      G4bool needResidual = true;
      int totN = 0;
      int totZ = 0;
      for ( G4int j = 0; j < int( products->size() ); j++ ) 
      {

         G4int jZ = (*products)[j].Z; 
         G4int jA = (*products)[j].A; 
         G4int jm = (*products)[j].m; 
         //TK 
         //We need coordination LEND *products)[j].m and G4IonTable(Z,A,m) 
         //Excitation energy of isomer level is might (probably) different each other.
         //

         //G4cout << "ZA = " << 1000 * (*products)[j].Z + (*products)[j].A << "  EK = "
         //     << (*products)[j].kineticEnergy
         //     << " px  " <<  (*products)[j].px
         //     << " py  " <<  (*products)[j].py
         //     << " pz  " <<  (*products)[j].pz
         //     << G4endl;

         iTotZ -= jZ;
         iTotA -= jA;

         G4DynamicParticle* theSec = new G4DynamicParticle;

         if ( jA == 1 && jZ == 1 )
         {
            theSec->SetDefinition( G4Proton::Proton() );
            totN += 1;
            totZ += 1;
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
               theSec->SetDefinition( G4IonTable::GetIonTable()->GetIon( jZ , jA , jm ) );
               totN += jA;
               totZ += jZ;
            }
            else 
            {
               theSec->SetDefinition( G4IonTable::GetIonTable()->GetIon( jZ , iA+aTrack.GetDefinition()->GetAtomicMass()-totN , jm ) );
               iTotZ -= jZ;
               iTotA -= iA+aTrack.GetDefinition()->GetAtomicMass()-totN;
               needResidual=false;
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

      if ( !( iTotZ == 0 && iTotA == 0 ) ) {

         if ( iTotZ >= 0 && iTotA > 0 ) {
            if ( needResidual ) {
               G4DynamicParticle* residual = new G4DynamicParticle;
               if ( iTotZ > 0 ) {
                  residual->SetDefinition( G4IonTable::GetIonTable()->GetIon( iTotZ , iTotA ) );
               } else if ( iTotA == 1 ) {
                  residual->SetDefinition( G4Neutron::Neutron() );
               } else {
                  //G4cout << "Charge or Baryon Number Error #3 iTotZ = " << iTotZ << ", iTotA = " << iTotA << G4endl;
                  ;
               }
               residual->SetMomentum( proj_p - psum );
               theResult->AddSecondary( residual );
            } else { 
               //G4cout << "Charge or Baryon Number Error #1 iTotZ = " << iTotZ << ", iTotA = " << iTotA << G4endl;
               ;
            }
         } else {

            if ( needResidual ) {
               //G4cout << "Charge or Baryon Number Error #2 iTotZ = " << iTotZ << ", iTotA = " << iTotA << G4endl;
               ;
            }
         }

      }

   } 
   else {
      //G4cout << "Using PreCompoundModel" << G4endl;
      if ( aTrack.GetDefinition() == G4Proton::Proton() || 
           aTrack.GetDefinition() == G4Neutron::Neutron() ) {
         theResult = preco->ApplyYourself( aTrack, aTarg );
      } else {
         return theResult;
      }
   }
   delete products;

   theResult->SetStatusChange( stopAndKill );

   return theResult; 

}
