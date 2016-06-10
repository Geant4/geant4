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
// 080602 Fix memory leaks by T. Koi 
// 081120 Add deltaT in signature of CalKinematicsOfBinaryCollisions
//        Add several required updating of Mean Filed 
//        Modified handling of absorption case by T. Koi
// 090126 Fix in absorption case by T. Koi  
// 090331 Fix for gamma participant by T. Koi 
//
#include "G4QMDCollision.hh"
#include "G4Scatterer.hh"
#include "G4Pow.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

G4QMDCollision::G4QMDCollision()
: fdeltar ( 4.0 )
, fbcmax0 ( 1.323142 ) // NN maximum impact parameter
, fbcmax1 ( 2.523 )    // others maximum impact parameter
// , sig0 ( 55 )   // NN cross section
//110617 fix for gcc 4.6 compilation warnings 
//, sig1 ( 200 )  // others cross section
, fepse ( 0.0001 )
{
   //These two pointers will be set through SetMeanField method
   theSystem=NULL;
   theMeanField=NULL;
   theScatterer = new G4Scatterer();
}

/*
G4QMDCollision::G4QMDCollision( const G4QMDCollision& obj )
: fdeltar ( obj.fdeltar )
, fbcmax0 ( obj.fbcmax0 ) // NN maximum impact parameter
, fbcmax1 ( obj.fbcmax1 )    // others maximum impact parameter
, fepse ( obj.fepse )
{
   
   if ( obj.theSystem != NULL ) {
      theSystem = new G4QMDSystem;
      *theSystem = *obj.theSystem;
   } else {
      theSystem = NULL;
   }
   if ( obj.theMeanField != NULL ) {
      theMeanField = new G4QMDMeanField;
      *theMeanField = *obj.theMeanField;
   } else {
      theMeanField = NULL;
   }
   theScatterer = new G4Scatterer();
   *theScatterer = *obj.theScatterer;
}

G4QMDCollision & G4QMDCollision::operator= ( const G4QMDCollision& obj)
{
   fdeltar = obj.fdeltar;
   fbcmax0 = obj.fbcmax1;
   fepse = obj.fepse;

   if ( obj.theSystem != NULL ) {
      delete theSystem;
      theSystem = new G4QMDSystem;
      *theSystem = *obj.theSystem;
   } else {
      theSystem = NULL;
   }
   if ( obj.theMeanField != NULL ) {
      delete theMeanField;
      theMeanField = new G4QMDMeanField;
      *theMeanField = *obj.theMeanField;
   } else {
      theMeanField = NULL;
   }
   delete theScatterer;
   theScatterer = new G4Scatterer();
   *theScatterer = *obj.theScatterer;

   return *this;
}
*/


G4QMDCollision::~G4QMDCollision()
{
   //if ( theSystem != NULL ) delete theSystem;
   //if ( theMeanField != NULL ) delete theMeanField;
   delete theScatterer;
}


void G4QMDCollision::CalKinematicsOfBinaryCollisions( G4double dt )
{
   G4double deltaT = dt; 

   G4int n = theSystem->GetTotalNumberOfParticipant(); 
//081118
   //G4int nb = 0;
   for ( G4int i = 0 ; i < n ; i++ )
   {
      theSystem->GetParticipant( i )->UnsetHitMark();
      theSystem->GetParticipant( i )->UnsetHitMark();
      //nb += theSystem->GetParticipant( i )->GetBaryonNumber();
   }
   //G4cout << "nb = " << nb << " n = " << n << G4endl;


//071101
   for ( G4int i = 0 ; i < n ; i++ )
   {

      //std::cout << i << " " << theSystem->GetParticipant( i )->GetDefinition()->GetParticleName() << " " << theSystem->GetParticipant( i )->GetPosition() << std::endl;

      if ( theSystem->GetParticipant( i )->GetDefinition()->IsShortLived() )
      {

         G4bool decayed = false; 

         const G4ParticleDefinition* pd0 = theSystem->GetParticipant( i )->GetDefinition();
         G4ThreeVector p0 = theSystem->GetParticipant( i )->GetMomentum();
         G4ThreeVector r0 = theSystem->GetParticipant( i )->GetPosition();

         G4LorentzVector p40 = theSystem->GetParticipant( i )->Get4Momentum();

         G4double eini = theMeanField->GetTotalPotential() + p40.e();

         G4int n0 = theSystem->GetTotalNumberOfParticipant(); 
         G4int i0 = 0; 

         G4bool isThisEnergyOK = false;

         G4int maximumNumberOfTrial=4;
         for ( G4int ii = 0 ; ii < maximumNumberOfTrial ; ii++ )
         { 

            //G4LorentzVector p4 = theSystem->GetParticipant( i )->Get4Momentum();
            G4LorentzVector p400 = p40;

            p400 *= GeV;
            //G4KineticTrack kt( theSystem->GetParticipant( i )->GetDefinition() , 0.0 , (theSystem->GetParticipant( i )->GetPosition())*fermi , p4 );
            G4KineticTrack kt( pd0 , 0.0 , r0*fermi , p400 );
            //std::cout << "G4KineticTrack " << i << " " <<  kt.GetDefinition()->GetParticleName() <<  kt.GetPosition() << std::endl;
            G4KineticTrackVector* secs = NULL;
            secs = kt.Decay();
            G4int id = 0;
            G4double et = 0;
            if ( secs )
            {
               for ( G4KineticTrackVector::iterator it 
                     = secs->begin() ; it != secs->end() ; it++ )
               {
/*
                  G4cout << "G4KineticTrack" 
                  << " " << (*it)->GetDefinition()->GetParticleName()
                  << " " << (*it)->Get4Momentum()
                  << " " << (*it)->GetPosition()/fermi
                  << G4endl;
*/
                   if ( id == 0 ) 
                   {
                      theSystem->GetParticipant( i )->SetDefinition( (*it)->GetDefinition() );
                      theSystem->GetParticipant( i )->SetMomentum( (*it)->Get4Momentum().v()/GeV );
                      theSystem->GetParticipant( i )->SetPosition( (*it)->GetPosition()/fermi );
                      //theMeanField->Cal2BodyQuantities( i ); 
                      et += (*it)->Get4Momentum().e()/GeV;
                   }
                   if ( id > 0 )
                   {
                      // Append end;
                      theSystem->SetParticipant ( new G4QMDParticipant ( (*it)->GetDefinition() , (*it)->Get4Momentum().v()/GeV , (*it)->GetPosition()/fermi ) );
                      et += (*it)->Get4Momentum().e()/GeV;
                      if ( id > 1 )
                      {
                         //081118
                         //G4cout << "G4QMDCollision id >2; id= " << id  << G4endl; 
                      }
                   }
                   id++; // number of daughter particles

                   delete *it;
               }

               theMeanField->Update();
               i0 = id-1; // 0 enter to i

               delete secs;
            }

//          EnergyCheck  

            G4double efin = theMeanField->GetTotalPotential() + et; 
            //std::cout <<  std::abs ( eini - efin ) - fepse << std::endl; 
//            std::cout <<  std::abs ( eini - efin ) - fepse*10 << std::endl; 
//                                           *10 TK  
            if ( std::abs ( eini - efin ) < fepse*10 ) 
            {
               // Energy OK 
               isThisEnergyOK = true;
               break; 
            }
            else
            {

               theSystem->GetParticipant( i )->SetDefinition( pd0 );
               theSystem->GetParticipant( i )->SetPosition( r0 );
               theSystem->GetParticipant( i )->SetMomentum( p0 );

               //for ( G4int i0i = 0 ; i0i < id-1 ; i0i++ )
               //160210 deletion must be done in descending order
               for ( G4int i0i = id-2 ; 0 <= i0i ; i0i-- ) {
                  //081118
                  //std::cout << "Decay Energitically Blocked deleteing " << i0i+n0  << std::endl;
                  theSystem->DeleteParticipant( i0i+n0 );
               }
               //081103
               theMeanField->Update();
            }

         }


//       Pauli Check 
         if ( isThisEnergyOK == true )
         {
            if ( theMeanField->IsPauliBlocked ( i ) != true ) 
            {

               G4bool allOK = true; 
               for ( G4int i0i = 0 ; i0i < i0 ; i0i++ )
               {
                  if ( theMeanField->IsPauliBlocked ( i0i+n0 ) == true )
                  {
                     allOK = false;
                     break;
                  } 
               }

               if ( allOK ) 
               {
                  decayed = true; //Decay Succeeded
               }
            }

         }
//       

         if ( decayed )
         {
            //081119
            //G4cout << "Decay Suceeded! " << std::endl;
            theSystem->GetParticipant( i )->SetHitMark();
            for ( G4int i0i = 0 ; i0i < i0 ; i0i++ )
            {
                theSystem->GetParticipant( i0i+n0 )->SetHitMark();
            }

         }
         else
         {

//          Decay Blocked and re-enter orginal participant;

            if ( isThisEnergyOK == true )  // for false case already done
            {

               theSystem->GetParticipant( i )->SetDefinition( pd0 );
               theSystem->GetParticipant( i )->SetPosition( r0 );
               theSystem->GetParticipant( i )->SetMomentum( p0 );

               for ( G4int i0i = 0 ; i0i < i0 ; i0i++ )
               {
                  //081118
                  //std::cout << "Decay Blocked deleteing " << i0i+n0  << std::endl;
                  //160210 adding commnet: deletion must be done in descending order
                  theSystem->DeleteParticipant( i0+n0-i0i-1 );
               }
               //081103
               theMeanField->Update();
            }

         }

      }  //shortlive
   }  // go next participant 
//071101


   n = theSystem->GetTotalNumberOfParticipant(); 

//081118
   //for ( G4int i = 1 ; i < n ; i++ )
   for ( G4int i = 1 ; i < theSystem->GetTotalNumberOfParticipant() ; i++ )
   {

      //std::cout << "Collision i " << i << std::endl;
      if ( theSystem->GetParticipant( i )->IsThisHit() ) continue; 

      G4ThreeVector ri =  theSystem->GetParticipant( i )->GetPosition();
      G4LorentzVector p4i =  theSystem->GetParticipant( i )->Get4Momentum();
      G4double rmi =  theSystem->GetParticipant( i )->GetMass();
      const G4ParticleDefinition* pdi =  theSystem->GetParticipant( i )->GetDefinition();
//090331 gamma 
      if ( pdi->GetPDGMass() == 0.0 ) continue;

      //std::cout << " p4i00 " << p4i << std::endl;
      for ( G4int j = 0 ; j < i ; j++ )
      {


/*
         G4cout << "Collision " << i << " " << theSystem->GetParticipant( i )->IsThisProjectile() << G4endl;
         G4cout << "Collision " << j << " " << theSystem->GetParticipant( j )->IsThisProjectile() << G4endl;
         G4cout << "Collision " << i << " " << theSystem->GetParticipant( i )->IsThisTarget() << G4endl;
         G4cout << "Collision " << j << " " << theSystem->GetParticipant( j )->IsThisTarget() << G4endl;
*/

         // Only 1 Collision allowed for each particle in a time step. 
         //081119
         if ( theSystem->GetParticipant( i )->IsThisHit() ) continue; 
         if ( theSystem->GetParticipant( j )->IsThisHit() ) continue; 

         //std::cout << "Collision " << i << " " << j << std::endl;

         // Do not allow collision between nucleons in target/projectile til its first collision.
         if ( theSystem->GetParticipant( i )->IsThisProjectile() )
         {
            if ( theSystem->GetParticipant( j )->IsThisProjectile() ) continue;
         }
         else if ( theSystem->GetParticipant( i )->IsThisTarget() )
         {
            if ( theSystem->GetParticipant( j )->IsThisTarget() ) continue;
         }


         G4ThreeVector rj =  theSystem->GetParticipant( j )->GetPosition();
         G4LorentzVector p4j =  theSystem->GetParticipant( j )->Get4Momentum();
         G4double rmj =  theSystem->GetParticipant( j )->GetMass();
         const G4ParticleDefinition* pdj =  theSystem->GetParticipant( j )->GetDefinition();
//090331 gamma 
         if ( pdj->GetPDGMass() == 0.0 ) continue;

         G4double rr2 = theMeanField->GetRR2( i , j );

//       Here we assume elab (beam momentum less than 5 GeV/n )
         if ( rr2 > fdeltar*fdeltar ) continue;

         //G4double s = (p4i+p4j)*(p4i+p4j);
         //G4double srt = std::sqrt ( s );

         G4double srt = std::sqrt( (p4i+p4j)*(p4i+p4j) );

         G4double cutoff = 0.0;
         G4double fbcmax = 0.0;
         //110617 fix for gcc 4.6 compilation warnings 
         //G4double sig = 0.0;

         if ( rmi < 0.94 && rmj < 0.94 ) 
         {
//          nucleon or pion case
            cutoff = rmi + rmj + 0.02; 
            fbcmax = fbcmax0;
            //110617 fix for gcc 4.6 compilation warnings 
            //sig = sig0;
         }
         else
         {
            cutoff = rmi + rmj; 
            fbcmax = fbcmax1;
            //110617 fix for gcc compilation warnings 
            //sig = sig1;
         }

         //std::cout << "Collision cutoff " << i << " " << j << " " << cutoff << std::endl;
         if ( srt < cutoff ) continue; 
        
         G4ThreeVector dr = ri - rj;
         G4double rsq = dr*dr;

         G4double pij = p4i*p4j; 
         G4double pidr = p4i.vect()*dr;
         G4double pjdr = p4j.vect()*dr;

         G4double aij = 1.0 - ( rmi*rmj /pij ) * ( rmi*rmj /pij ); 
         G4double bij = pidr / rmi - pjdr*rmi/pij;
         G4double cij = rsq + ( pidr / rmi ) * ( pidr / rmi );
         G4double brel = std::sqrt ( std::abs ( cij - bij*bij/aij ) );
 
         if ( brel > fbcmax ) continue;
         //std::cout << "collisions3 " << std::endl;
     
         G4double bji = -pjdr/rmj + pidr * rmj /pij;
 
         G4double ti = ( pidr/rmi - bij / aij ) * p4i.e() / rmi;
         G4double tj = (-pjdr/rmj - bji / aij ) * p4j.e() / rmj;


/*
         G4cout << "collisions4  p4i " << p4i << G4endl;
         G4cout << "collisions4  ri " << ri << G4endl;
         G4cout << "collisions4  p4j " << p4j << G4endl;
         G4cout << "collisions4  rj " << rj << G4endl;
         G4cout << "collisions4  dr " << dr << G4endl;
         G4cout << "collisions4  pij " << pij << G4endl;
         G4cout << "collisions4  aij " << aij << G4endl;
         G4cout << "collisions4  bij bji " << bij << " " << bji << G4endl;
         G4cout << "collisions4  pidr pjdr " << pidr << " " << pjdr << G4endl;
         G4cout << "collisions4  p4i.e() p4j.e() " << p4i.e() << " " << p4j.e() << G4endl;
         G4cout << "collisions4  rmi rmj " << rmi << " " << rmj << G4endl;
         G4cout << "collisions4 " << ti << " " << tj << G4endl;
*/
         if ( std::abs ( ti + tj ) > deltaT ) continue;
         //std::cout << "collisions4 " << std::endl;

         G4ThreeVector beta = ( p4i + p4j ).boostVector();

         G4LorentzVector p = p4i;
         G4LorentzVector p4icm = p.boost( p.findBoostToCM ( p4j ) );
         G4ThreeVector pcm = p4icm.vect();
         
         G4double prcm = pcm.mag();

         if ( prcm <= 0.00001 ) continue; 
         //std::cout << "collisions5 " << std::endl;

         G4bool energetically_forbidden = !( CalFinalStateOfTheBinaryCollision ( i , j ) ); // Use Geant4 Collision Library
         //G4bool energetically_forbidden = !( CalFinalStateOfTheBinaryCollisionJQMD ( sig , cutoff , pcm , prcm , srt, beta , gamma , i , j ) ); // JQMD Elastic 

/*
         G4bool pauli_blocked = false;
         if ( energetically_forbidden == false ) // result true 
         { 
            if ( theMeanField->IsPauliBlocked ( i ) == true || theMeanField->IsPauliBlocked ( j ) == true ) 
            {
               pauli_blocked = true;
               //std::cout << "G4QMDRESULT Collsion Pauli Blocked " << std::endl;
            }
         }
         else
         {
            if ( theMeanField->IsPauliBlocked ( i ) == true || theMeanField->IsPauliBlocked ( j ) == true ) 
               pauli_blocked = false;
            //std::cout << "G4QMDRESULT Collsion Blocked " << std::endl;
         } 
*/

/*
            G4cout << "G4QMDRESULT Collsion initial p4 i and j " 
                      << p4i << " " << p4j
                      << G4endl;
*/
//       081118
         //if ( energetically_forbidden == true || pauli_blocked == true )
         if ( energetically_forbidden == true )
         {

            //G4cout << " energetically_forbidden  " << G4endl;
//          Collsion not allowed then re enter orginal participants 
//          Now only momentum, becasuse we only consider elastic scattering of nucleons

            theSystem->GetParticipant( i )->SetMomentum( p4i.vect() );
            theSystem->GetParticipant( i )->SetDefinition( pdi );
            theSystem->GetParticipant( i )->SetPosition( ri );

            theSystem->GetParticipant( j )->SetMomentum( p4j.vect() );
            theSystem->GetParticipant( j )->SetDefinition( pdj );
            theSystem->GetParticipant( j )->SetPosition( rj );

            theMeanField->Cal2BodyQuantities( i ); 
            theMeanField->Cal2BodyQuantities( j ); 

         }
         else 
         {

            
           G4bool absorption = false; 
           if ( n == theSystem->GetTotalNumberOfParticipant()+1 ) absorption = true;
           if ( absorption ) 
           {
              //G4cout << "Absorption happend " << G4endl; 
              i = i-1; 
              n = n-1;
           } 
              
//          Collsion allowed (really happened) 

            // Unset Projectile/Target flag
            theSystem->GetParticipant( i )->UnsetInitialMark();
            if ( !absorption ) theSystem->GetParticipant( j )->UnsetInitialMark();

            theSystem->GetParticipant( i )->SetHitMark();
            if ( !absorption ) theSystem->GetParticipant( j )->SetHitMark();

            theSystem->IncrementCollisionCounter();

/*
            G4cout << "G4QMDRESULT Collsion Really Happened between " 
                      << i << " and " << j 
                      << G4endl;
            G4cout << "G4QMDRESULT Collsion initial p4 i and j " 
                      << p4i << " " << p4j
                      << G4endl;
            G4cout << "G4QMDRESULT Collsion after p4 i and j " 
                      << theSystem->GetParticipant( i )->Get4Momentum()
                      << " " 
                      << theSystem->GetParticipant( j )->Get4Momentum()
                      << G4endl;
            G4cout << "G4QMDRESULT Collsion Diff " 
                      << p4i + p4j - theSystem->GetParticipant( i )->Get4Momentum() - theSystem->GetParticipant( j )->Get4Momentum()
                      << G4endl;
            G4cout << "G4QMDRESULT Collsion initial r i and j " 
                      << ri << " " << rj
                      << G4endl;
            G4cout << "G4QMDRESULT Collsion after r i and j " 
                      << theSystem->GetParticipant( i )->GetPosition()
                      << " " 
                      << theSystem->GetParticipant( j )->GetPosition()
                      << G4endl;
*/
             

         }

      }

   }


}



G4bool G4QMDCollision::CalFinalStateOfTheBinaryCollision( G4int i , G4int j )
{

//081103
   //G4cout << "CalFinalStateOfTheBinaryCollision " << i << " " << j << " " << theSystem->GetTotalNumberOfParticipant() << G4endl;

   G4bool result = false;
   G4bool energyOK = false; 
   G4bool pauliOK = false; 
   G4bool abs = false;
   G4QMDParticipant* absorbed = NULL;  

   G4LorentzVector p4i = theSystem->GetParticipant( i )->Get4Momentum();
   G4LorentzVector p4j = theSystem->GetParticipant( j )->Get4Momentum();

//071031

   G4double epot = theMeanField->GetTotalPotential();

   G4double eini = epot + p4i.e() + p4j.e();

//071031
   // will use KineticTrack
   const G4ParticleDefinition* pdi0 =theSystem->GetParticipant( i )->GetDefinition();
   const G4ParticleDefinition* pdj0 =theSystem->GetParticipant( j )->GetDefinition();
   G4LorentzVector p4i0 = p4i*GeV;
   G4LorentzVector p4j0 = p4j*GeV;
   G4ThreeVector ri0 = ( theSystem->GetParticipant( i )->GetPosition() )*fermi;
   G4ThreeVector rj0 = ( theSystem->GetParticipant( j )->GetPosition() )*fermi;

   for ( G4int iitry = 0 ; iitry < 4 ; iitry++ )
   {

      abs = false;

      G4KineticTrack kt1( pdi0 , 0.0 , ri0 , p4i0 );
      G4KineticTrack kt2( pdj0 , 0.0 , rj0 , p4j0 );

      G4LorentzVector p4ix_new; 
      G4LorentzVector p4jx_new; 
      G4KineticTrackVector* secs = NULL;
      secs = theScatterer->Scatter( kt1 , kt2 );

      //std::cout << "G4QMDSCATTERER BEFORE " << kt1.GetDefinition()->GetParticleName() << " " << kt1.Get4Momentum()/GeV << " " << kt1.GetPosition()/fermi << std::endl;
      //std::cout << "G4QMDSCATTERER BEFORE " << kt2.GetDefinition()->GetParticleName() << " " << kt2.Get4Momentum()/GeV << " " << kt2.GetPosition()/fermi << std::endl;
      //std::cout << "THESCATTERER " << theScatterer->GetCrossSection ( kt1 , kt2 )/millibarn << " " << elastic << " " << sig << std::endl;


      if ( secs )
      {
         G4int iti = 0;
         if (  secs->size() == 2 )
         {
            for ( G4KineticTrackVector::iterator it 
                = secs->begin() ; it != secs->end() ; it++ )
            {
               if ( iti == 0 ) 
               {
                  theSystem->GetParticipant( i )->SetDefinition( (*it)->GetDefinition() );
                  p4ix_new = (*it)->Get4Momentum()/GeV;
                  //std::cout << "THESCATTERER " << (*it)->GetDefinition()->GetParticleName() << std::endl;
                  theSystem->GetParticipant( i )->SetMomentum( p4ix_new.v() );
               } 
               if ( iti == 1 )
               {
                  theSystem->GetParticipant( j )->SetDefinition( (*it)->GetDefinition() );
                  p4jx_new = (*it)->Get4Momentum()/GeV;
                  //std::cout << "THESCATTERER " << p4jx_new.e()-p4jx_new.m() << std::endl;
                  theSystem->GetParticipant( j )->SetMomentum( p4jx_new.v() );
               }
               //std::cout << "G4QMDSCATTERER AFTER " << (*it)->GetDefinition()->GetParticleName() << " " << (*it)->Get4Momentum()/GeV << std::endl;
               iti++;
            }
         }
         else if ( secs->size() == 1 )
         {
//081118
            abs = true;
            //G4cout << "G4QMDCollision pion absrorption " << secs->front()->GetDefinition()->GetParticleName() << G4endl;
            //secs->front()->Decay();
            theSystem->GetParticipant( i )->SetDefinition( secs->front()->GetDefinition() );
            p4ix_new = secs->front()->Get4Momentum()/GeV;
            theSystem->GetParticipant( i )->SetMomentum( p4ix_new.v() );

         } 

//081118
         if ( secs->size() > 2 ) 
         {

            G4cout << "G4QMDCollision secs size > 2;  " << secs->size() << G4endl;

            for ( G4KineticTrackVector::iterator it 
                = secs->begin() ; it != secs->end() ; it++ )
            {
               G4cout << "G4QMDSCATTERER AFTER " << (*it)->GetDefinition()->GetParticleName() << " " << (*it)->Get4Momentum()/GeV << G4endl;
            }

         }

         // deleteing KineticTrack
         for ( G4KineticTrackVector::iterator it 
               = secs->begin() ; it != secs->end() ; it++ )
         {  
            delete *it;
         }

         delete secs;
      }
//071031

      if ( !abs )
      { 
         theMeanField->Cal2BodyQuantities( i ); 
         theMeanField->Cal2BodyQuantities( j ); 
      } 
      else
      {
         absorbed = theSystem->EraseParticipant( j ); 
         theMeanField->Update();
      }

      epot = theMeanField->GetTotalPotential();

      G4double efin = epot + p4ix_new.e() + p4jx_new.e(); 

      //std::cout << "Collision NEW epot " << i << " " << j << " " << epot << " " << std::abs ( eini - efin ) - fepse << std::endl;

/*
      G4cout << "Collision efin " << i << " " << j << " " << efin << G4endl;
      G4cout << "Collision " << i << " " << j << " " << std::abs ( eini - efin ) << " " << fepse << G4endl;
      G4cout << "Collision " << std::abs ( eini - efin ) << " " << fepse << G4endl;
*/

//071031
      if ( std::abs ( eini - efin ) < fepse ) 
      {
         // Collison OK 
         //std::cout << "collisions6" << std::endl;
         //std::cout << "collisions before " << p4i << " " << p4j << std::endl;
         //std::cout << "collisions after " << theSystem->GetParticipant( i )->Get4Momentum() << " " << theSystem->GetParticipant( j )->Get4Momentum() << std::endl;
         //std::cout << "collisions dif " << ( p4i + p4j ) - ( theSystem->GetParticipant( i )->Get4Momentum() + theSystem->GetParticipant( j )->Get4Momentum() ) << std::endl;
         //std::cout << "collisions before " << ri0/fermi << " " << rj0/fermi << std::endl;
         //std::cout << "collisions after " << theSystem->GetParticipant( i )->GetPosition() << " " << theSystem->GetParticipant( j )->GetPosition() << std::endl;
         energyOK = true;
         break;
      }
      else
      {
         //G4cout << "Energy Not OK " << G4endl;
         if ( abs )
         {
            //G4cout << "TKDB reinsert j " << G4endl;
            theSystem->InsertParticipant( absorbed , j );   
            theMeanField->Update();
         }
         // do not need reinsert in no absroption case 
      }
//071031
   }

// Energetically forbidden collision

   if ( energyOK )
   {
      // Pauli Check 
      //G4cout << "Pauli Checking " << theSystem->GetTotalNumberOfParticipant() << G4endl;
      if ( !abs ) 
      {
         if ( !( theMeanField->IsPauliBlocked ( i ) == true || theMeanField->IsPauliBlocked ( j ) == true ) ) 
         {
            //G4cout << "Binary Collision Happen " << theSystem->GetTotalNumberOfParticipant() << G4endl;
            pauliOK = true;
         }
      }
      else 
      {
         //if ( theMeanField->IsPauliBlocked ( i ) == false ) 
         //090126                            i-1 cause jth is erased
         if ( theMeanField->IsPauliBlocked ( i-1 ) == false ) 
         { 
            //G4cout << "Absorption Happen " << theSystem->GetTotalNumberOfParticipant() << G4endl;
            delete absorbed;
            pauliOK = true;
         }
      }
      

      if ( pauliOK ) 
      {
         result = true;
      }
      else
      {
         //G4cout << "Pauli Blocked" << G4endl;
         if ( abs )
         {
            //G4cout << "TKDB reinsert j pauli block" << G4endl;
            theSystem->InsertParticipant( absorbed , j );   
            theMeanField->Update(); 
         }
      }
   }

   return result;

} 



G4bool G4QMDCollision::CalFinalStateOfTheBinaryCollisionJQMD( G4double sig , G4double cutoff , G4ThreeVector pcm , G4double prcm , G4double srt , G4ThreeVector beta , G4double gamma , G4int i , G4int j )
{

   //G4cout << "CalFinalStateOfTheBinaryCollisionJQMD" << G4endl;

   G4bool result = true;

   G4LorentzVector p4i =  theSystem->GetParticipant( i )->Get4Momentum();
   G4double rmi =  theSystem->GetParticipant( i )->GetMass();
   G4int zi =  theSystem->GetParticipant( i )->GetChargeInUnitOfEplus();

   G4LorentzVector p4j =  theSystem->GetParticipant( j )->Get4Momentum();
   G4double rmj =  theSystem->GetParticipant( j )->GetMass();
   G4int zj =  theSystem->GetParticipant( j )->GetChargeInUnitOfEplus();

   G4double pr = prcm;

   G4double c2  = pcm.z()/pr;
    
   G4double csrt = srt - cutoff;

   //G4double pri = prcm;
   //G4double prf = sqrt ( 0.25 * srt*srt -rm2 );
    
   G4double asrt = srt - rmi - rmj;
   G4double pra = prcm;
   


   G4double elastic = 0.0; 

   if ( zi == zj )
   {
      if ( csrt < 0.4286 )
      {
         elastic = 35.0 / ( 1. + csrt * 100.0 )  +  20.0;
      }
      else
      {
         elastic = ( - std::atan( ( csrt - 0.4286 ) * 1.5 - 0.8 )
                 *   2. / pi + 1.0 ) * 9.65 + 7.0;
      }         
   }
   else
   {
      if ( csrt < 0.4286 )
      {
         elastic = 28.0 / ( 1. + csrt * 100.0 )  +  27.0;
      }
      else
      {
         elastic = ( - std::atan( ( csrt - 0.4286 ) * 1.5 - 0.8 )
                 *   2. / pi + 1.0 ) * 12.34 + 10.0;
      }         
   }

//   std::cout << "Collision csrt " << i << " " << j << " " << csrt << std::endl;
//   std::cout << "Collision elstic " << i << " " << j << " " << elastic << std::endl;


//   std::cout << "Collision sig " << i << " " << j  << " " << sig << std::endl;
   if ( G4UniformRand() > elastic / sig ) 
   { 
      //std::cout << "Inelastic " << std::endl; 
      //std::cout << "elastic/sig " << elastic/sig << std::endl; 
      return result; 
   }
   else
   {
      //std::cout << "Elastic " << std::endl; 
   } 
//   std::cout << "Collision ELSTIC " << i << " " << j << std::endl;

   
   G4double as = G4Pow::GetInstance()->powN ( 3.65 * asrt , 6 );
   G4double a = 6.0 * as / (1.0 + as);
   G4double ta = -2.0 * pra*pra;
   G4double x = G4UniformRand(); 
   G4double t1 = G4Log( (1-x) * G4Exp(2.*a*ta) + x )  /  a;
   G4double c1 = 1.0 - t1/ta;
 
   if( std::abs(c1) > 1.0 ) c1 = 2.0 * x - 1.0;

/*
   G4cout << "Collision as " << i << " " << j << " " << as << G4endl;
   G4cout << "Collision a " << i << " " << j << " " << a << G4endl;
   G4cout << "Collision ta " << i << " " << j << " " << ta << G4endl;
   G4cout << "Collision x " << i << " " << j << " " << x << G4endl;
   G4cout << "Collision t1 " << i << " " << j << " " << t1 << G4endl;
   G4cout << "Collision c1 " << i << " " << j << " " << c1 << G4endl;
*/
   t1 = 2.0*pi*G4UniformRand(); 
//   std::cout << "Collision t1 " << i << " " << j << " " << t1 << std::endl;
   G4double t2 = 0.0;
   if ( pcm.x() == 0.0 && pcm.y() == 0 )
   {
      t2 = 0.0;
   }
   else 
   {
      t2 = std::atan2( pcm.y() , pcm.x() );
   }
//      std::cout << "Collision t2 " << i << " " << j << " " << t2 << std::endl;

   G4double s1 = std::sqrt ( 1.0 - c1*c1 );
   G4double s2 = std::sqrt ( 1.0 - c2*c2 );
 
   G4double ct1 = std::cos(t1);      
   G4double st1 = std::sin(t1);      

   G4double ct2 = std::cos(t2);      
   G4double st2 = std::sin(t2);      

   G4double ss = c2*s1*ct1 + s2*c1;

   pcm.setX( pr * ( ss*ct2 - s1*st1*st2) );
   pcm.setY( pr * ( ss*st2 + s1*st1*ct2) );
   pcm.setZ( pr * ( c1*c2 - s1*s2*ct1) );

// std::cout << "Collision pcm " << i << " " << j << " " << pcm << std::endl;

   G4double epot = theMeanField->GetTotalPotential();

   G4double eini = epot + p4i.e() + p4j.e();
   G4double etwo = p4i.e() + p4j.e();

/*
   G4cout << "Collision epot " << i << " " << j << " " << epot << G4endl;
   G4cout << "Collision eini " << i << " " << j << " " << eini << G4endl;
   G4cout << "Collision etwo " << i << " " << j << " " << etwo << G4endl;
*/


   for ( G4int itry = 0 ; itry < 4 ; itry++ )
   {

      G4double eicm = std::sqrt ( rmi*rmi + pcm*pcm );
      G4double pibeta = pcm*beta;

      G4double trans = gamma * ( gamma * pibeta / ( gamma + 1 ) + eicm );
   
      G4ThreeVector pi_new = beta*trans + pcm;
   
      G4double ejcm = std::sqrt ( rmj*rmj + pcm*pcm );
      trans = gamma * ( gamma * pibeta / ( gamma + 1 ) + ejcm );

      G4ThreeVector pj_new = beta*trans - pcm;

//
// Delete old 
// Add new Particitipants
//
// Now only change momentum ( Beacuse we only have elastic sctter of nucleon
// In future Definition also will be change 
//

      theSystem->GetParticipant( i )->SetMomentum( pi_new );
      theSystem->GetParticipant( j )->SetMomentum( pj_new );

      G4double pi_new_e = (theSystem->GetParticipant( i )->Get4Momentum()).e();
      G4double pj_new_e = (theSystem->GetParticipant( j )->Get4Momentum()).e();

      theMeanField->Cal2BodyQuantities( i ); 
      theMeanField->Cal2BodyQuantities( j ); 

      epot = theMeanField->GetTotalPotential();

      G4double efin = epot + pi_new_e + pj_new_e ; 

      //std::cout << "Collision NEW epot " << i << " " << j << " " << epot << " " << std::abs ( eini - efin ) - fepse << std::endl;
/*
      G4cout << "Collision efin " << i << " " << j << " " << efin << G4endl;
      G4cout << "Collision " << i << " " << j << " " << std::abs ( eini - efin ) << " " << fepse << G4endl;
      G4cout << "Collision " << std::abs ( eini - efin ) << " " << fepse << G4endl;
*/

//071031
      if ( std::abs ( eini - efin ) < fepse ) 
      {
	 // Collison OK 
         //std::cout << "collisions6" << std::endl;
         //std::cout << "collisions before " << p4i << " " << p4j << std::endl;
         //std::cout << "collisions after " << theSystem->GetParticipant( i )->Get4Momentum() << " " << theSystem->GetParticipant( j )->Get4Momentum() << std::endl;
         //std::cout << "collisions dif " << ( p4i + p4j ) - ( theSystem->GetParticipant( i )->Get4Momentum() + theSystem->GetParticipant( j )->Get4Momentum() ) << std::endl;
         //std::cout << "collisions before " << rix/fermi << " " << rjx/fermi << std::endl;
         //std::cout << "collisions after " << theSystem->GetParticipant( i )->GetPosition() << " " << theSystem->GetParticipant( j )->GetPosition() << std::endl;
      }
//071031

         if ( std::abs ( eini - efin ) < fepse ) return result;  // Collison OK 
      
         G4double cona = ( eini - efin + etwo ) / gamma;
         G4double fac2 = 1.0 / ( 4.0 * cona*cona * pr*pr ) *
                       ( ( cona*cona - ( rmi*rmi + rmj*rmj ) )*( cona*cona - ( rmi*rmi + rmj*rmj ) )
                       - 4.0 * rmi*rmi * rmj*rmj );

         if ( fac2 > 0 ) 
         {
            G4double fact = std::sqrt ( fac2 ); 
            pcm = fact*pcm;
         }


   }

// Energetically forbidden collision
   result = false;

   return result;

} 
