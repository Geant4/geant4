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

#include "G4LENDCapture.hh"
#include "G4Fragment.hh"
#include "G4PhotonEvaporation.hh"
#include "G4SystemOfUnits.hh"
#include "G4Nucleus.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
  
G4HadFinalState * G4LENDCapture::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTarg )
{

   G4double temp = aTrack.GetMaterial()->GetTemperature();

   //G4int iZ = int ( aTarg.GetZ() );
   //G4int iA = int ( aTarg.GetN() );
   //migrate to integer A and Z (GetN_asInt returns number of neutrons in the nucleus since this) 
   G4int iZ = aTarg.GetZ_asInt();
   G4int iA = aTarg.GetA_asInt();
   G4int iM = 0;
   if ( aTarg.GetIsotope() != NULL ) {
      iM = aTarg.GetIsotope()->Getm();
   }

   G4double ke = aTrack.GetKineticEnergy();

   G4HadFinalState* theResult = &theParticleChange;
   theResult->Clear();

   G4GIDI_target* aTarget = get_target_from_map( lend_manager->GetNucleusEncoding( iZ , iA , iM ) );
   if ( aTarget == NULL ) return returnUnchanged( aTrack , theResult );
   std::vector<G4GIDI_Product>* products = aTarget->getCaptureFinalState( ke*MeV, temp, MyRNG, NULL );

   G4int ipZ = aTrack.GetDefinition()->GetAtomicNumber();
   G4int ipA = aTrack.GetDefinition()->GetAtomicMass();

   G4bool needResidual=true;

   G4ThreeVector p(0,0,0);
   if ( products != NULL ) 
   {

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

         if ( jZ == iZ + ipZ && jA == iA + ipA ) needResidual = false;

         G4ThreeVector dp((*products)[j].px,(*products)[j].py,(*products)[j].pz);
         p += dp;
          
         G4DynamicParticle* theSec = new G4DynamicParticle;

         if ( jA == 1 && jZ == 1 ) {
            theSec->SetDefinition( G4Proton::Proton() );
            totN += 1;
         }
         else if ( jA == 1 && jZ == 0 )
         {
            theSec->SetDefinition( G4Neutron::Neutron() );
            totN += 1;
         } 
         else if ( jZ > 0 ) {
            if ( jA != 0 )
            {
               theSec->SetDefinition( G4IonTable::GetIonTable()->GetIon( jZ , jA , iM ) );
               totN += jA;
            }
            else 
            {
               theSec->SetDefinition( G4IonTable::GetIonTable()->GetIon( jZ , iA+1-totN , iM ) );
            }
         } 
         else {
            theSec->SetDefinition( G4Gamma::Gamma() );
         } 

         theSec->SetMomentum( G4ThreeVector( (*products)[j].px*MeV , (*products)[j].py*MeV , (*products)[j].pz*MeV ) );

/*
         if ( dp.mag() == 0 ) 
         {
            //theSec->SetMomentum( -p*MeV ); 
         }
*/

         theResult->AddSecondary( theSec );
      } 
   }
   else 
   {

      //For the case data does not provide final states 

      //G4cout << "products != NULL; iZ = " << iZ << ", iA = " << iA << G4endl; 

      // TK comment 
      // aTarg->ReturnTargetParticle()->Get4Momentum has trouble, thus we use following  
      G4Fragment nucleus( iA + ipA , iZ + ipZ ,  aTrack.Get4Momentum() + G4LorentzVector( G4ThreeVector(0,0,0) , G4IonTable::GetIonTable()->GetIon( iZ + ipZ , iA )->GetPDGMass() ) );
      G4PhotonEvaporation photonEvaporation;
      photonEvaporation.SetICM( TRUE );
      G4FragmentVector* products_from_PE = photonEvaporation.BreakItUp(nucleus);
      G4FragmentVector::iterator it;

      for ( it = products_from_PE->begin(); it != products_from_PE->end(); it++)
      {
         if ( (*it)->GetZ_asInt() == iZ + ipZ &&  (*it)->GetA_asInt() == iA + ipA )  needResidual = false;
         G4DynamicParticle* theSec = new G4DynamicParticle;
         if ( (*it)->GetParticleDefinition() != NULL ) {
            //G4cout << (*it)->GetParticleDefinition()->GetParticleName() << G4endl; 
            theSec->SetDefinition( (*it)->GetParticleDefinition() );
            theSec->Set4Momentum( (*it)->GetMomentum() );
         } else {
            //G4cout << (*it)->GetZ_asInt() << " " << (*it)->GetA_asInt() << G4endl;
            theSec->SetDefinition( G4IonTable::GetIonTable()->GetIon( (*it)->GetZ_asInt() , (*it)->GetA_asInt() ) );
            theSec->Set4Momentum( (*it)->GetMomentum() );
         }
         theResult->AddSecondary( theSec );
      }
   }
   
   //if necessary, generate residual nucleus
   if ( needResidual ) {
      G4DynamicParticle* residual = new G4DynamicParticle;
      residual->SetDefinition( G4IonTable::GetIonTable()->GetIon( iZ + ipZ , iA + ipA ) );
      residual->SetMomentum( -p*MeV ); 
      theResult->AddSecondary( residual );
   }

   delete products;

   theResult->SetStatusChange( stopAndKill );

   return theResult; 

}
