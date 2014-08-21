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

   G4GIDI_target* aTarget = usedTarget_map.find( lend_manager->GetNucleusEncoding( iZ , iA , iM ) )->second->GetTarget();
   std::vector<G4GIDI_Product>* products = aTarget->getOthersFinalState( ke*MeV, temp, NULL, NULL );
   if ( products != NULL ) 
   {

      G4ThreeVector psum(0);

      int totN = 0;
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
               theSec->SetDefinition( G4IonTable::GetIonTable()->GetIon( jZ , jA , jm ) );
               totN += jA;
            }
            else 
            {
               theSec->SetDefinition( G4IonTable::GetIonTable()->GetIon( jZ , iA+1-totN , jm ) );
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
