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
#include "G4QMDCollision.hh"
#include "G4ParticleDefinition.hh"
#include "G4Scatterer.hh"
#include "Randomize.hh"
#include "G4PionZero.hh"

G4QMDCollision::G4QMDCollision()
: deltar ( 4 )
, bcmax0 ( 1.323142 ) // NN maximum impact parameter
, bcmax1 ( 2.523 )    // others maximum impact parameter
, sig0 ( 55 )   // NN cross section
, sig1 ( 200 )  // others cross section
, epse ( 0.0001 )
{
   theScatterer = new G4Scatterer();
}



G4QMDCollision::~G4QMDCollision()
{
   delete theScatterer;
}


void G4QMDCollision::CalKinematicsOfBinaryCollisions()
{


   G4int n = theSystem->GetTotalNumberOfParticipant(); 

//071101
   for ( G4int i = 0 ; i < n ; i++ )
   {
      //std::cout << i << " " << theSystem->GetParticipant( i )->GetDefinition()->GetParticleName() << " " << theSystem->GetParticipant( i )->GetPosition() << std::endl;
      if ( theSystem->GetParticipant( i )->GetDefinition()->IsShortLived() )
      {
         G4ParticleDefinition* pd0 = theSystem->GetParticipant( i )->GetDefinition();
         G4ThreeVector p0 = theSystem->GetParticipant( i )->GetMomentum();
         G4ThreeVector r0 = theSystem->GetParticipant( i )->GetPosition();

         G4LorentzVector p40 = theSystem->GetParticipant( i )->Get4Momentum();

         G4double epot = theMeanField->GetTotalPotential();
         G4double eini = epot + p40.e();

         G4int n0 = theSystem->GetTotalNumberOfParticipant(); 
         G4int i0 = 0; 
G4bool isThisEnergyOK = false;
         for ( G4int ii = 0 ; ii < 4 ; ii++ )
{ 

         //G4LorentzVector p4 = theSystem->GetParticipant( i )->Get4Momentum();
         G4LorentzVector p400 = p40;

         p400 *= GeV;
         //G4KineticTrack kt( theSystem->GetParticipant( i )->GetDefinition() , 0.0 , (theSystem->GetParticipant( i )->GetPosition())*fermi , p4 );
         G4KineticTrack kt( pd0 , 0.0 , r0*fermi , p400 );
//         std::cout << "G4KineticTrack " << i << " " <<  kt.GetDefinition()->GetParticleName() <<  kt.GetPosition() << std::endl;
         G4KineticTrackVector* secs = NULL;
         secs = kt.Decay();
         G4int id = 0;
         G4double et = 0;
         if ( secs )
         {
            for ( G4KineticTrackVector::iterator it 
                  = secs->begin() ; it != secs->end() ; it++ )
            {
//              std::cout << "G4KineticTrack" 
//                 << " " << (*it)->GetDefinition()->GetParticleName()
//                 << " " << (*it)->Get4Momentum()
//                 << " " << (*it)->GetPosition()/fermi
//                         << std::endl;
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
 //                     std::cout << "NAGISA id >2; id= " << id  << std::endl; 
                   }
                }
                id++; 

                delete *it;
            }
            theMeanField->SetSystem ( theSystem );
            i0 = id-1; // 0 enter to i
         }

//       EnergyCheck  

         G4double epot = theMeanField->GetTotalPotential();
         G4double efin = epot + et; 
         //std::cout <<  std::abs ( eini - efin ) - epse << std::endl; 
//         std::cout <<  std::abs ( eini - efin ) - epse*10 << std::endl; 
//071031
//                                        *10 TK  
         if ( std::abs ( eini - efin ) < epse*10 ) 
         {
            // Energy OK 
//            std::cout << "Decay Succeeded Energy OK" << std::endl;
            isThisEnergyOK = true;
            break; 
         }
         else
         {
            for ( G4int i0i = 0 ; i0i < id-1 ; i0i++ )
            {
//               std::cout << "Decay Energitically Blocked deleteing " << i0i+n0  << std::endl;
               theSystem->DeleteParticipant( i0i+n0 );
            }
         }
}

//       Pauli Check 
         if ( isThisEnergyOK == true )
         {
//          if ( theMeanField->IsPauliBlocked ( i ) != true ) 
            {
               bool allOK = true; 
               for ( G4int i0i = 0 ; i0i < i0 ; i0i++ )
               {
                  if ( theMeanField->IsPauliBlocked ( i0i+n0 ) == true )
                  {
                     allOK = false;
                     break;
                  } 
               }

//               if ( allOK ) std::cout << "Decay Succeeded" << std::endl;
               if ( allOK ) continue; //Do not Pauli Blocked
            }
         }
//       

//         std::cout << "Decay Blocked" << std::endl;
         theSystem->GetParticipant( i )->SetDefinition( pd0 );
         theSystem->GetParticipant( i )->SetPosition( r0 );
         theSystem->GetParticipant( i )->SetMomentum( p0 );

         if ( isThisEnergyOK == true )
         {
         for ( G4int i0i = 0 ; i0i < i0 ; i0i++ )
         {
//            std::cout << "Decay Blocked deleteing " << i0i+n0  << std::endl;
            theSystem->DeleteParticipant( i0i+n0 );
         }
         }
         
      }
   }
//071101


   n = theSystem->GetTotalNumberOfParticipant(); 
   //std::cout << "Collision n " << n << std::endl;

   std::vector< G4bool > isCollided ( n , false );

   for ( G4int i = 1 ; i < n ; i++ )
   {

      //std::cout << "Collision i " << i << std::endl;

      G4ThreeVector ri =  theSystem->GetParticipant( i )->GetPosition();
      G4LorentzVector p4i =  theSystem->GetParticipant( i )->Get4Momentum();
      G4double rmi =  theSystem->GetParticipant( i )->GetMass();
      G4ParticleDefinition* pdi =  theSystem->GetParticipant( i )->GetDefinition();

      //std::cout << " p4i00 " << p4i << std::endl;
      for ( G4int j = 0 ; j < i ; j++ )
      {
//         std::cout << "Collision " << i << " " << j << std::endl;

/*
         std::cout << "Collision " << i << " " << theSystem->GetParticipant( i )->IsThisProjectile() << std::endl;
         std::cout << "Collision " << j << " " << theSystem->GetParticipant( j )->IsThisProjectile() << std::endl;
         std::cout << "Collision " << i << " " << theSystem->GetParticipant( i )->IsThisTarget() << std::endl;
         std::cout << "Collision " << j << " " << theSystem->GetParticipant( j )->IsThisTarget() << std::endl;
*/

         // Only 1 Collision allowed for each particle in a time step. 
         if ( isCollided[ i ] == true ) continue;
         if ( isCollided[ j ] == true ) continue;

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
         G4ParticleDefinition* pdj =  theSystem->GetParticipant( j )->GetDefinition();

         G4double rr2 = theMeanField->GetRR2( i , j );

//       Here we assume elab (beam momentum less than 5 GeV/n )
         if ( rr2 > deltar*deltar ) continue;

         G4double s = (p4i+p4j)*(p4i+p4j);

         G4double srt = std::sqrt ( s );

         G4double cutoff = 0.0;
         G4double bcmax = 0.0;
         G4double sig = 0.0;

         if ( rmi < 0.94 && rmj < 0.94 ) 
         {
//          nucleon or pion case
            cutoff = rmi + rmj + 0.02; 
            bcmax = bcmax0;
            sig = sig0;
         }
         else
         {
            cutoff = rmi + rmj; 
            bcmax = bcmax1;
            sig = sig1;
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
 
         if ( brel > bcmax ) continue;
         //std::cout << "collisions3 " << std::endl;
     
         G4double bji = -pjdr/rmj + pidr * rmj /pij;
 
         G4double ti = ( pidr/rmi - bij / aij ) * p4i.e() / rmi;
         G4double tj = (-pjdr/rmj - bji / aij ) * p4j.e() / rmj;

         G4double deltaT = 0.0;
         deltaT = 1.0; // TK  

/*
         std::cout << "collisions4  p4i " << p4i << std::endl;
         std::cout << "collisions4  ri " << ri << std::endl;
         std::cout << "collisions4  p4j " << p4j << std::endl;
         std::cout << "collisions4  rj " << rj << std::endl;
         std::cout << "collisions4  dr " << dr << std::endl;
         std::cout << "collisions4  pij " << pij << std::endl;
         std::cout << "collisions4  aij " << aij << std::endl;
         std::cout << "collisions4  bij bji " << bij << " " << bji << std::endl;
         std::cout << "collisions4  pidr pjdr " << pidr << " " << pjdr << std::endl;
         std::cout << "collisions4  p4i.e() p4j.e() " << p4i.e() << " " << p4j.e() << std::endl;
         std::cout << "collisions4  rmi rmj " << rmi << " " << rmj << std::endl;
         std::cout << "collisions4 " << ti << " " << tj << std::endl;
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

         G4bool pauli_blocked = false;
         if ( energetically_forbidden != true ) 
         { 
            if ( theMeanField->IsPauliBlocked ( i ) == true || theMeanField->IsPauliBlocked ( j ) == true ) 
            {
               pauli_blocked = true;
               //std::cout << "G4QMDRESULT Collsion Pauli Blocked " << std::endl;
            }
         }
         else
         {
            //std::cout << "G4QMDRESULT Collsion Blocked " << std::endl;
         } 

/*
            std::cout << "G4QMDRESULT Collsion initial p4 i and j " 
                      << p4i << " " << p4j
                      << std::endl;
*/

 
         if ( energetically_forbidden == true ||  pauli_blocked == true )
         {
//       Collsion not allowed then re enter orginal participants 
//       Now only momentum, becasuse we only consider elastic scattering of nucleons

            theSystem->GetParticipant( i )->SetMomentum( p4i.vect() );
            theSystem->GetParticipant( i )->SetDefinition( pdi );
            theSystem->GetParticipant( i )->SetPosition( ri );
            theSystem->GetParticipant( j )->SetMomentum( p4j.vect() );
            theSystem->GetParticipant( j )->SetDefinition( pdj );
            theSystem->GetParticipant( j )->SetPosition( rj );

         }
         else 
         {
//       Collsion allowed (really happened) 

            // Unset Projectile/Target flag
            theSystem->GetParticipant( i )->UnsetInitialMark();
            theSystem->GetParticipant( j )->UnsetInitialMark();

            isCollided[ i ] = true; 
            isCollided[ j ] = true; 

            theSystem->IncrementCollisionCounter();

/*
            std::cout << "G4QMDRESULT Collsion Really Happened between " 
                      << i << " and " << j 
                      << std::endl;
            std::cout << "G4QMDRESULT Collsion initial p4 i and j " 
                      << p4i << " " << p4j
                      << std::endl;
            std::cout << "G4QMDRESULT Collsion after p4 i and j " 
                      << theSystem->GetParticipant( i )->Get4Momentum()
                      << " " 
                      << theSystem->GetParticipant( j )->Get4Momentum()
                      << std::endl;
            std::cout << "G4QMDRESULT Collsion Diff " 
                      << p4i + p4j - theSystem->GetParticipant( i )->Get4Momentum() - theSystem->GetParticipant( j )->Get4Momentum()
                      << std::endl;
            std::cout << "G4QMDRESULT Collsion initial r i and j " 
                      << ri << " " << rj
                      << std::endl;
            std::cout << "G4QMDRESULT Collsion after r i and j " 
                      << theSystem->GetParticipant( i )->GetPosition()
                      << " " 
                      << theSystem->GetParticipant( j )->GetPosition()
                      << std::endl;
*/

         }

//         theMeanField

      }
   }

//071106
   n = theSystem->GetTotalNumberOfParticipant(); 
   G4bool isThisModefied = false;
   for ( G4int i = 0 ; i < n ; i++ )
   {
      if ( theSystem->GetParticipant( i )->GetDefinition() == G4PionZero::PionZero() )
      {
         if ( theSystem->GetParticipant( i )->GetPosition().mag() > 1.0e9 )
         {
//            std::cout << "Deleting " << i << " " << theSystem->GetParticipant( i )->GetPosition().mag() << std::endl;
            theSystem->DeleteParticipant( i );
            isThisModefied = true;
         }
      } 
   }
   if ( isThisModefied == true ) theMeanField->SetSystem ( theSystem );
//071106

}



G4bool G4QMDCollision::CalFinalStateOfTheBinaryCollision( G4int i , G4int j )
{

   //G4cout << "CalFinalStateOfTheBinaryCollision " << G4endl;

   G4bool result = true;

   G4LorentzVector p4i =  theSystem->GetParticipant( i )->Get4Momentum();
   G4LorentzVector p4j =  theSystem->GetParticipant( j )->Get4Momentum();

//071031
   // will use KineticTrack
   G4LorentzVector p4ix = p4i*GeV;
   G4LorentzVector p4jx = p4j*GeV;
   G4ThreeVector rix = (theSystem->GetParticipant( i )->GetPosition())*fermi; 
   G4ThreeVector rjx = (theSystem->GetParticipant( j )->GetPosition())*fermi; 
//071031

   G4double epot = theMeanField->GetTotalPotential();

   G4double eini = epot + p4i.e() + p4j.e();


//071031
   G4ParticleDefinition* pdi0 =theSystem->GetParticipant( i )->GetDefinition();
   G4ParticleDefinition* pdj0 =theSystem->GetParticipant( j )->GetDefinition();
   G4ThreeVector ri0 =(theSystem->GetParticipant( i )->GetPosition())*fermi;
   G4ThreeVector rj0 =(theSystem->GetParticipant( j )->GetPosition())*fermi;

   for ( G4int iitry = 0 ; iitry < 4 ; iitry++ )
   {

      G4KineticTrack kt1( pdi0 , 0.0 , ri0 , p4ix );
      G4KineticTrack kt2( pdj0 , 0.0 , rj0 , p4jx );
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
         else
         {
            //std::cout << "NAGISA pion absrorption " << secs->front()->GetDefinition()->GetParticleName() << std::endl;
            //secs->front()->Decay();
            theSystem->GetParticipant( i )->SetDefinition( secs->front()->GetDefinition() );
            p4ix_new = secs->front()->Get4Momentum()/GeV;
            theSystem->GetParticipant( i )->SetMomentum( p4ix_new.v() );
            
              //std::cout << "THESCATTERER " << (*it)->GetDefinition()->GetParticleName() << std::endl;
              p4jx_new( 0 ); 
              //theSystem->GetParticipant( j )->SetDefinition( G4Gamma::Gamma() );
              //theSystem->GetParticipant( j )->SetDefinition( G4Neutron::Neutron() );
              theSystem->GetParticipant( j )->SetDefinition( G4PionZero::PionZero() );
              theSystem->GetParticipant( j )->SetMomentum( G4ThreeVector( G4UniformRand() )*eV );
              theSystem->GetParticipant( j )->SetPosition( G4ThreeVector( 1000, 1000, 1000 )*km );

         } 

         if ( secs->size() > 2 ) std::cout << "NAGISA secs size > 2;  " << secs->size() << std::endl;

         // deleteing KineticTrack
         for ( G4KineticTrackVector::iterator it 
               = secs->begin() ; it != secs->end() ; it++ )
         {  
            delete *it;
         }
      }
//071031

      theMeanField->Cal2BodyQuantities( i ); 
      theMeanField->Cal2BodyQuantities( j ); 

      epot = theMeanField->GetTotalPotential();

      G4double efin = epot + p4ix_new.e() + p4jx_new.e(); 

      //std::cout << "Collision NEW epot " << i << " " << j << " " << epot << " " << std::abs ( eini - efin ) - epse << std::endl;

/*
      std::cout << "Collision efin " << i << " " << j << " " << efin << std::endl;
      std::cout << "Collision " << i << " " << j << " " << std::abs ( eini - efin ) << " " << epse << std::endl;
      std::cout << "Collision " << std::abs ( eini - efin ) << " " << epse << std::endl;
*/

//071031
      if ( std::abs ( eini - efin ) < epse ) 
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

      if ( std::abs ( eini - efin ) < epse ) return result;  // Collison OK 

   }

//  Energetically forbidden collision
    result = false;

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

   
   G4double as = std::pow ( 3.65 * asrt , 6 );
   G4double a = 6.0 * as / (1.0 + as);
   G4double ta = -2.0 * pra*pra;
   G4double x = G4UniformRand(); 
   G4double t1 = std::log( (1-x) * std::exp(2.*a*ta) + x )  /  a;
   G4double c1 = 1.0 - t1/ta;
 
   if( std::abs(c1) > 1.0 ) c1 = 2.0 * x - 1.0;

/*
   std::cout << "Collision as " << i << " " << j << " " << as << std::endl;
   std::cout << "Collision a " << i << " " << j << " " << a << std::endl;
   std::cout << "Collision ta " << i << " " << j << " " << ta << std::endl;
   std::cout << "Collision x " << i << " " << j << " " << x << std::endl;
   std::cout << "Collision t1 " << i << " " << j << " " << t1 << std::endl;
   std::cout << "Collision c1 " << i << " " << j << " " << c1 << std::endl;
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
   std::cout << "Collision epot " << i << " " << j << " " << epot << std::endl;
   std::cout << "Collision eini " << i << " " << j << " " << eini << std::endl;
   std::cout << "Collision etwo " << i << " " << j << " " << etwo << std::endl;
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

      //std::cout << "Collision NEW epot " << i << " " << j << " " << epot << " " << std::abs ( eini - efin ) - epse << std::endl;
/*
      std::cout << "Collision efin " << i << " " << j << " " << efin << std::endl;
      std::cout << "Collision " << i << " " << j << " " << std::abs ( eini - efin ) << " " << epse << std::endl;
      std::cout << "Collision " << std::abs ( eini - efin ) << " " << epse << std::endl;
*/

//071031
      if ( std::abs ( eini - efin ) < epse ) 
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

         if ( std::abs ( eini - efin ) < epse ) return result;  // Collison OK 
      
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
