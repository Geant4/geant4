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
// 081120 Add Update by T. Koi
//

#include <map>
#include <algorithm>
#include <numeric>

#include <cmath>
#include <CLHEP/Random/Stat.h>

#include "G4QMDMeanField.hh"
#include "G4QMDParameters.hh"
#include "G4Exp.hh"
#include "G4Pow.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

G4QMDMeanField::G4QMDMeanField()
: rclds ( 4.0 )    // distance for cluster judgement        
, epsx ( -20.0 )   // gauss term      
, epscl ( 0.0001 ) // coulomb term     
, irelcr ( 1 )     
{

   G4QMDParameters* parameters = G4QMDParameters::GetInstance(); 
   wl = parameters->Get_wl();
   cl = parameters->Get_cl();
   rho0 = parameters->Get_rho0(); 
   hbc = parameters->Get_hbc();
   gamm = parameters->Get_gamm();

   cpw = parameters->Get_cpw();
   cph = parameters->Get_cph();
   cpc = parameters->Get_cpc();

   c0 = parameters->Get_c0();
   c3 = parameters->Get_c3();
   cs = parameters->Get_cs();

// distance
   c0w = 1.0/4.0/wl;
   //c3w = 1.0/4.0/wl; //no need
   c0sw = std::sqrt( c0w );
   clw = 2.0 / std::sqrt ( 4.0 * pi * wl );

// graduate
   c0g = - c0 / ( 2.0 * wl );
   c3g = - c3 / ( 4.0 * wl ) * gamm;
   csg = - cs / ( 2.0 * wl );
   pag = gamm - 1;

   system = NULL; // will be set through SetSystem method
}



G4QMDMeanField::~G4QMDMeanField()
{
   ;
}



void G4QMDMeanField::SetSystem ( G4QMDSystem* aSystem ) 
{ 

   //std::cout << "QMDMeanField SetSystem" << std::endl;

   system = aSystem; 

   G4int n = system->GetTotalNumberOfParticipant();
  
   pp2.clear();
   rr2.clear();
   rbij.clear();
   rha.clear();
   rhe.clear();
   rhc.clear();

   rr2.resize( n );
   pp2.resize( n );
   rbij.resize( n );
   rha.resize( n );
   rhe.resize( n );
   rhc.resize( n );

   for ( int i = 0 ; i < n ; i++ )
   {
      rr2[i].resize( n );
      pp2[i].resize( n );
      rbij[i].resize( n );
      rha[i].resize( n );
      rhe[i].resize( n );
      rhc[i].resize( n );
   }


   ffr.clear();
   ffp.clear();
   rh3d.clear();

   ffr.resize( n );
   ffp.resize( n );
   rh3d.resize( n );

   Cal2BodyQuantities();

}

void G4QMDMeanField::SetNucleus ( G4QMDNucleus* aNucleus ) 
{

   //std::cout << "QMDMeanField SetNucleus" << std::endl;

   SetSystem( aNucleus );

   G4double totalPotential = GetTotalPotential(); 
   aNucleus->SetTotalPotential( totalPotential );

   aNucleus->CalEnergyAndAngularMomentumInCM();

}



void G4QMDMeanField::Cal2BodyQuantities()
{

   if ( system->GetTotalNumberOfParticipant() < 2 ) return;

   for ( G4int j = 1 ; j < system->GetTotalNumberOfParticipant() ; j++ )
   {

      G4ThreeVector rj = system->GetParticipant( j )->GetPosition();
      G4LorentzVector p4j = system->GetParticipant( j )->Get4Momentum();

      for ( G4int i = 0 ; i < j ; i++ )
      {

         G4ThreeVector ri = system->GetParticipant( i )->GetPosition();
         G4LorentzVector p4i = system->GetParticipant( i )->Get4Momentum();

         G4ThreeVector rij = ri - rj;
         G4ThreeVector pij = (p4i - p4j).v();
         G4LorentzVector p4ij = p4i - p4j;
         G4ThreeVector bij = ( p4i + p4j ).boostVector();
         G4double gammaij = ( p4i + p4j ).gamma();

         G4double eij = ( p4i + p4j ).e();

         G4double rbrb = rij*bij;
//         G4double bij2 = bij*bij;
         G4double rij2 = rij*rij;
         G4double pij2 = pij*pij;

         rbrb = irelcr * rbrb;
         G4double  gamma2_ij = gammaij*gammaij;


         rr2[i][j] = rij2 + gamma2_ij * rbrb*rbrb;
         rr2[j][i] = rr2[i][j];

         rbij[i][j] = gamma2_ij * rbrb;
         rbij[j][i] = - rbij[i][j];

         pp2[i][j] = pij2
                   + irelcr * ( - G4Pow::GetInstance()->powN ( p4i.e() - p4j.e() , 2 )
                   + gamma2_ij * G4Pow::GetInstance()->powN ( ( ( p4i.m2() - p4j.m2() ) / eij ) , 2 ) );


         pp2[j][i] = pp2[i][j];

//       Gauss term

         G4double expa1 = - rr2[i][j] * c0w;

         G4double rh1;
         if ( expa1 > epsx )
         {
            rh1 = G4Exp( expa1 );
         }
         else
         {
            rh1 = 0.0;
         }

         G4int ibry = system->GetParticipant(i)->GetBaryonNumber();
         G4int jbry = system->GetParticipant(j)->GetBaryonNumber();


         rha[i][j] = ibry*jbry*rh1;
         rha[j][i] = rha[i][j];

//       Coulomb terms

         G4double rrs2 = rr2[i][j] + epscl;
         G4double rrs = std::sqrt ( rrs2 );

         G4int icharge = system->GetParticipant(i)->GetChargeInUnitOfEplus();
         G4int jcharge = system->GetParticipant(j)->GetChargeInUnitOfEplus();

         G4double xerf = 0.0;
         // T. K. add this protection. 5.8 is good enough for double
         if ( rrs*c0sw < 5.8 ) {
            //erf = G4RandStat::erf ( rrs*c0sw );
            //Restore to CLHEP for avoiding compilation error in MT
            //erf = CLHEP::HepStat::erf ( rrs*c0sw );
            //Use cmath 
#if defined WIN32-VC
            xerf = CLHEP::HepStat::erf ( rrs*c0sw );
#else
            xerf = erf ( rrs*c0sw );
#endif
         } else {
            xerf = 1.0;
         }

         G4double erfij = xerf/rrs;


         rhe[i][j] = icharge*jcharge * erfij;

         rhe[j][i] = rhe[i][j];

         rhc[i][j] = icharge*jcharge * ( - erfij + clw * rh1 ) / rrs2;

         rhc[j][i] = rhc[i][j];

      }  // i
   }  // j
}



void G4QMDMeanField::Cal2BodyQuantities( G4int i )
{

   //std::cout << "Cal2BodyQuantities " << i << std::endl;

   G4ThreeVector ri = system->GetParticipant( i )->GetPosition();  
   G4LorentzVector p4i = system->GetParticipant( i )->Get4Momentum();  


   for ( G4int j = 0 ; j < system->GetTotalNumberOfParticipant() ; j ++ )
   {
      if ( j == i ) continue; 

      G4ThreeVector rj = system->GetParticipant( j )->GetPosition();  
      G4LorentzVector p4j = system->GetParticipant( j )->Get4Momentum();  

         G4ThreeVector rij = ri - rj;
         G4ThreeVector pij = (p4i - p4j).v();
         G4LorentzVector p4ij = p4i - p4j;
         G4ThreeVector bij = ( p4i + p4j ).boostVector();
         G4double gammaij = ( p4i + p4j ).gamma();

         G4double eij = ( p4i + p4j ).e();

         G4double rbrb = rij*bij;
//         G4double bij2 = bij*bij;
         G4double rij2 = rij*rij;
         G4double pij2 = pij*pij;

         rbrb = irelcr * rbrb;
         G4double  gamma2_ij = gammaij*gammaij;

/*
      G4double rbrb = 0.0;
      G4double beta2_ij = 0.0;
      G4double rij2 = 0.0;
      G4double pij2 = 0.0;

//    
      G4LorentzVector p4ip4j = p4i + p4j;
      G4double eij = p4ip4j.e();

      G4ThreeVector r = ri - rj;  
      G4LorentzVector p4 = p4i - p4j;  

      rbrb = r.x()*p4ip4j.x()/eij
           +  r.y()*p4ip4j.y()/eij
           +  r.z()*p4ip4j.z()/eij;

      beta2_ij = ( p4ip4j.x()*p4ip4j.x() + p4ip4j.y()*p4ip4j.y() + p4ip4j.z()*p4ip4j.z() ) / ( eij*eij ); 
      rij2 = r*r;
      pij2 = p4.v()*p4.v();

      rbrb = irelcr * rbrb;

      G4double gamma2_ij = 1 / ( 1 - beta2_ij ); 
*/

      rr2[i][j] = rij2 + gamma2_ij * rbrb*rbrb;
      rr2[j][i] = rr2[i][j];

      rbij[i][j] = gamma2_ij * rbrb;
      rbij[j][i] = - rbij[i][j];
      
      pp2[i][j] = pij2
                + irelcr * ( - G4Pow::GetInstance()->powN ( p4i.e() - p4j.e() , 2 )
                + gamma2_ij * G4Pow::GetInstance()->powN ( ( ( p4i.m2() - p4j.m2() ) / eij ) , 2 ) );

      pp2[j][i] = pp2[i][j];

//    Gauss term

      G4double expa1 = - rr2[i][j] * c0w;  

      G4double rh1;
      if ( expa1 > epsx ) 
      {
         rh1 = G4Exp( expa1 );
      }
      else 
      {
         rh1 = 0.0;  
      } 

      G4int ibry = system->GetParticipant(i)->GetBaryonNumber();  
      G4int jbry = system->GetParticipant(j)->GetBaryonNumber();  


      rha[i][j] = ibry*jbry*rh1;  
      rha[j][i] = rha[i][j]; 

//    Coulomb terms

      G4double rrs2 = rr2[i][j] + epscl; 
      G4double rrs = std::sqrt ( rrs2 ); 

      G4int icharge = system->GetParticipant(i)->GetChargeInUnitOfEplus();
      G4int jcharge = system->GetParticipant(j)->GetChargeInUnitOfEplus();

      G4double xerf = 0.0;
      // T. K. add this protection. 5.8 is good enough for double
      if ( rrs*c0sw < 5.8 ) {
         //xerf = G4RandStat::erf ( rrs*c0sw );
         //Use cmath 
#if defined WIN32-VC
         xerf = CLHEP::HepStat::erf ( rrs*c0sw );
#else
         xerf = erf ( rrs*c0sw );
#endif
      } else {
         xerf = 1.0;
      }

      G4double erfij = xerf/rrs; 
      

      rhe[i][j] = icharge*jcharge * erfij;

      rhe[j][i] = rhe[i][j]; 

//    G4double clw;

      rhc[i][j] = icharge*jcharge * ( - erfij + clw * rh1 ) / rrs2;

      rhc[j][i] = rhc[i][j]; 

   }

}



void G4QMDMeanField::CalGraduate()
{

   ffr.resize( system->GetTotalNumberOfParticipant() );
   ffp.resize( system->GetTotalNumberOfParticipant() );
   rh3d.resize( system->GetTotalNumberOfParticipant() );

   for ( G4int i = 0 ; i < system->GetTotalNumberOfParticipant() ; i ++ )
   {
      G4double rho3 = 0.0;
      for ( G4int j = 0 ; j < system->GetTotalNumberOfParticipant() ; j ++ )
      {
         rho3 += rha[j][i]; 
      }
      rh3d[i] = G4Pow::GetInstance()->powA ( rho3 , pag ); 
   }


   for ( G4int i = 0 ; i < system->GetTotalNumberOfParticipant() ; i ++ )
   {

      G4ThreeVector ri = system->GetParticipant( i )->GetPosition();  
      G4LorentzVector p4i = system->GetParticipant( i )->Get4Momentum();  

      G4ThreeVector betai = p4i.v()/p4i.e();
      
//    R-JQMD
      G4double Vi = GetPotential( i );
      G4double p_zero = std::sqrt( p4i.e()*p4i.e() + 2*p4i.m()*Vi);
      G4ThreeVector betai_R = p4i.v()/p_zero;
      G4double mi_R = p4i.m()/p_zero;
//
      ffr[i] = betai_R;
      ffp[i] = G4ThreeVector( 0.0 );

      if ( false ) 
      {
         ffr[i] = betai;
         mi_R = 1.0;
      }

      for ( G4int j = 0 ; j < system->GetTotalNumberOfParticipant() ; j ++ )
      {

         G4ThreeVector rj = system->GetParticipant( j )->GetPosition();  
         G4LorentzVector p4j = system->GetParticipant( j )->Get4Momentum();  

         G4double eij = p4i.e() + p4j.e(); 

         G4int icharge = system->GetParticipant(i)->GetChargeInUnitOfEplus();
         G4int jcharge = system->GetParticipant(j)->GetChargeInUnitOfEplus();

         G4int inuc = system->GetParticipant(i)->GetNuc();
         G4int jnuc = system->GetParticipant(j)->GetNuc();

         G4double ccpp = c0g * rha[j][i]
                       + c3g * rha[j][i] * ( rh3d[j] + rh3d[i] )
                       + csg * rha[j][i] * jnuc * inuc
                           * ( 1. - 2. * std::abs( jcharge - icharge ) )
                       + cl * rhc[j][i];
         ccpp *= mi_R;

/*
           G4cout << c0g << " " << c3g << " " << csg << " " << cl << G4endl;
           G4cout << "ccpp " << i << " " << j << " " << ccpp << G4endl;
           G4cout << "rha[j][i] " << rha[j][i] << G4endl;
           G4cout << "rh3d " << rh3d[j] << " " << rh3d[i] << G4endl;
           G4cout << "rhc[j][i] " << rhc[j][i] << G4endl;
*/

         G4double grbb = - rbij[j][i];
         G4double ccrr = grbb * ccpp / eij;

/*
           G4cout << "ccrr " << ccrr << G4endl;
           G4cout << "grbb " << grbb << G4endl;
*/


         G4ThreeVector rij = ri - rj;   
         G4ThreeVector betaij =  ( p4i + p4j ).v()/eij;   

         G4ThreeVector cij = betaij - betai;   

         ffr[i] = ffr[i] + 2*ccrr* ( rij + grbb*cij );
            
         ffp[i] = ffp[i] - 2*ccpp* ( rij + grbb*betaij );

      }
   }

   //std::cout << "gradu 0 " << ffr[0] << " " << ffp[0] << std::endl;
   //std::cout << "gradu 1 " << ffr[1] << " " << ffp[1] << std::endl;

}



G4double G4QMDMeanField::GetPotential( G4int i )
{
   G4int n = system->GetTotalNumberOfParticipant();

   G4double rhoa = 0.0;
   G4double rho3 = 0.0;
   G4double rhos = 0.0;
   G4double rhoc = 0.0;


   G4int icharge = system->GetParticipant(i)->GetChargeInUnitOfEplus();
   G4int inuc = system->GetParticipant(i)->GetNuc();

   for ( G4int j = 0 ; j < n ; j ++ )
   {
      G4int jcharge = system->GetParticipant(j)->GetChargeInUnitOfEplus();
      G4int jnuc = system->GetParticipant(j)->GetNuc();

      rhoa += rha[j][i];
      rhoc += rhe[j][i];
      rhos += rha[j][i] * jnuc * inuc
                * ( 1 - 2 * std::abs ( jcharge - icharge ) );
   }

   rho3 = G4Pow::GetInstance()->powA ( rhoa , gamm );

   G4double potential = c0 * rhoa
                      + c3 * rho3
                      + cs * rhos
                      + cl * rhoc;

   return potential;
}



G4double G4QMDMeanField::GetTotalPotential()
{

   G4int n = system->GetTotalNumberOfParticipant();

   std::vector < G4double > rhoa ( n , 0.0 ); 
   std::vector < G4double > rho3 ( n , 0.0 ); 
   std::vector < G4double > rhos ( n , 0.0 ); 
   std::vector < G4double > rhoc ( n , 0.0 ); 
   

   for ( G4int i = 0 ; i < n ; i ++ )
   {
      G4int icharge = system->GetParticipant(i)->GetChargeInUnitOfEplus();
      G4int inuc = system->GetParticipant(i)->GetNuc();

      for ( G4int j = 0 ; j < n ; j ++ )
      {
         G4int jcharge = system->GetParticipant(j)->GetChargeInUnitOfEplus();
         G4int jnuc = system->GetParticipant(j)->GetNuc();

         rhoa[i] += rha[j][i];
         rhoc[i] += rhe[j][i];
         rhos[i] += rha[j][i] * jnuc * inuc 
                   * ( 1 - 2 * std::abs ( jcharge - icharge ) );
      }

      rho3[i] = G4Pow::GetInstance()->powA ( rhoa[i] , gamm );
   }

   G4double potential = c0 * std::accumulate( rhoa.begin() , rhoa.end() , 0.0 ) 
                      + c3 * std::accumulate( rho3.begin() , rho3.end() , 0.0 ) 
                      + cs * std::accumulate( rhos.begin() , rhos.end() , 0.0 ) 
                      + cl * std::accumulate( rhoc.begin() , rhoc.end() , 0.0 );

   return potential;

}



G4double G4QMDMeanField::calPauliBlockingFactor( G4int i )
{

   G4double pf = 0.0;
// i is supposed beyond total number of Participant()
   G4int icharge = system->GetParticipant(i)->GetChargeInUnitOfEplus();

   for ( G4int j = 0 ; j < system->GetTotalNumberOfParticipant() ; j++ )
   {

      G4int jcharge = system->GetParticipant(j)->GetChargeInUnitOfEplus();
      G4int jnuc = system->GetParticipant(j)->GetNuc();

      if ( jcharge == icharge && jnuc == 1 )
      {

/*
   G4cout << "Pauli i j " << i << " " << j << G4endl;
   G4cout << "Pauli icharge " << icharge << G4endl;
   G4cout << "Pauli jcharge " << jcharge << G4endl;
*/
         G4double expa = -rr2[i][j]*cpw;   


         if ( expa > epsx ) 
         {
            expa = expa - pp2[i][j]*cph;
/*
   G4cout << "Pauli cph " << cph << G4endl;
   G4cout << "Pauli pp2 " << pp2[i][j] << G4endl;
   G4cout << "Pauli expa " <<  expa  << G4endl;
   G4cout << "Pauli epsx " <<  epsx  << G4endl;
*/
            if ( expa > epsx ) 
            { 
//   std::cout << "Pauli phase " <<  pf  << std::endl;
               pf = pf + G4Exp ( expa );
            }
         }
      }
      
   }


   pf  = ( pf - 1.0 ) * cpc;

   //std::cout << "Pauli pf " << pf << std::endl;

   return pf;
    
}



G4bool G4QMDMeanField::IsPauliBlocked( G4int i )
{
    G4bool result = false; 
    
    if ( system->GetParticipant( i )->GetNuc() == 1 )
    {
       G4double pf = calPauliBlockingFactor( i );
       G4double rand = G4UniformRand(); 
       if ( pf > rand ) result = true;
    }

    return result; 
}



void G4QMDMeanField::DoPropagation( G4double dt )
{
 
   G4double cc2 = 1.0; 
   G4double cc1 = 1.0 - cc2; 
   G4double cc3 = 1.0 / 2.0 / cc2; 

   G4double dt3 = dt * cc3;
   G4double dt1 = dt * ( cc1 - cc3 );
   G4double dt2 = dt * cc2;

   CalGraduate(); 

   G4int n = system->GetTotalNumberOfParticipant();

// 1st Step 

   std::vector< G4ThreeVector > f0r, f0p;
   f0r.resize( n );
   f0p.resize( n );

   for ( G4int i = 0 ; i < n ; i++ )
   {
      G4ThreeVector ri = system->GetParticipant( i )->GetPosition();  
      G4ThreeVector p3i = system->GetParticipant( i )->GetMomentum();  

      ri += dt3* ffr[i];
      p3i += dt3* ffp[i];

      f0r[i] = ffr[i];
      f0p[i] = ffp[i];
      
      system->GetParticipant( i )->SetPosition( ri );  
      system->GetParticipant( i )->SetMomentum( p3i );  

//    we do not need set total momentum by ourselvs
   }

// 2nd Step
   Cal2BodyQuantities(); 
   CalGraduate(); 

   for ( G4int i = 0 ; i < n ; i++ )
   {
      G4ThreeVector ri = system->GetParticipant( i )->GetPosition();  
      G4ThreeVector p3i = system->GetParticipant( i )->GetMomentum();  

      ri += dt1* f0r[i] + dt2* ffr[i];
      p3i += dt1* f0p[i] + dt2* ffp[i];

      system->GetParticipant( i )->SetPosition( ri );  
      system->GetParticipant( i )->SetMomentum( p3i );  

//    we do not need set total momentum by ourselvs
   }

   Cal2BodyQuantities(); 

}



std::vector< G4QMDNucleus* > G4QMDMeanField::DoClusterJudgment()
{

   //std::cout << "MeanField DoClusterJudgemnt" << std::endl;

   Cal2BodyQuantities(); 

   G4double cpf2 = G4Pow::GetInstance()->A23 ( 1.5 * pi*pi * G4Pow::GetInstance()->powA ( 4.0 * pi * wl , -1.5 ) ) * hbc * hbc;
   G4double rcc2 = rclds*rclds;

   G4int n = system->GetTotalNumberOfParticipant();
   std::vector < G4double > rhoa;
   rhoa.resize ( n );

   for ( G4int i = 0 ; i < n ; i++ )
   {
      rhoa[i] = 0.0;

      if ( system->GetParticipant( i )->GetBaryonNumber() == 1 )  
      {
      for ( G4int j = 0 ; j < n ; j++ )
      {
         if ( system->GetParticipant( j )->GetBaryonNumber() == 1 )  
         rhoa[i] += rha[i][j];
      }
      }

      rhoa[i] = G4Pow::GetInstance()->A13 ( rhoa[i] + 1 );

   }

// identification of the cluster

   std::map < G4int , std::vector < G4int > > cluster_map;
   std::vector < G4bool > is_already_belong_some_cluster;

   //         cluster_id   participant_id
   std::multimap < G4int , G4int > comb_map; 
   std::multimap < G4int , G4int > assign_map; 
   assign_map.clear(); 

   std::vector < G4int > mascl;
   std::vector < G4int > num;
   mascl.resize ( n );
   num.resize ( n );
   is_already_belong_some_cluster.resize ( n );

   std::vector < G4int > is_assigned_to ( n , -1 );
   std::multimap < G4int , G4int > clusters;

   for ( G4int i = 0 ; i < n ; i++ )
   {
      mascl[i] = 1;
      num[i] = 1;

      is_already_belong_some_cluster[i] = false;
   }


   G4int ichek = 1;


   G4int id = 0;
   G4int cluster_id = -1;  
   for ( G4int i = 0 ; i < n-1 ; i++ )
   { 

      G4bool hasThisCompany = false;
//    Check only for bryons?
//      std::cout << "Check Baryon " << i << std::endl;

      if ( system->GetParticipant( i )->GetBaryonNumber() == 1 )  
      {

//      if ( is_already_belong_some_cluster[i] != true )  
//      {
      //G4int j1 = ichek + 1;
      G4int j1 = i + 1;
      for ( G4int j = j1 ; j < n ; j++ )
      {

         std::vector < G4int > cluster_participants;
         if ( system->GetParticipant( j )->GetBaryonNumber() == 1 )  
         {
         G4double rdist2 = rr2[ i ][ j ];
         G4double pdist2 = pp2[ i ][ j ];
         //G4double rdist2 = rr2[ num[i] ][ num[j] ];
         //G4double pdist2 = pp2[ num[i] ][ num[j] ];
         G4double pcc2 = cpf2 
                       * ( rhoa[ i ] + rhoa[ j ] )
                       * ( rhoa[ i ] + rhoa[ j ] );

//       Check phase space: close enough?
         if ( rdist2 < rcc2 && pdist2 < pcc2 ) 
         {

/*
            G4cout << "G4QMDRESULT "
                   << i << " " << j << " " << id << " " 
                   << is_assigned_to [ i ] << " " << is_assigned_to [ j ] 
                   << G4endl;
*/

            if ( is_assigned_to [ j ] == -1 )
            {
               if ( is_assigned_to [ i ] == -1 )
               {
                  if ( clusters.size() != 0 )
                  {
                     id = clusters.rbegin()->first + 1;  
                     //std::cout << "id is increare " << id << std::endl;
                  }
                  else
                  {
                     id = 0;
                  }
                  clusters.insert ( std::multimap<G4int,G4int>::value_type ( id , i ) );
                  is_assigned_to [ i ] = id; 
                  clusters.insert ( std::multimap<G4int,G4int>::value_type ( id , j ) );
                  is_assigned_to [ j ] = id; 
               }
               else
               {
                  clusters.insert ( std::multimap<G4int,G4int>::value_type ( is_assigned_to [ i ] , j ) );
                  is_assigned_to [ j ] = is_assigned_to [ i ]; 
               }
            } 
            else
            {
//             j is already belong to some cluester 
               if ( is_assigned_to [ i ] == -1 )
               {
                  clusters.insert ( std::multimap<G4int,G4int>::value_type ( is_assigned_to [ j ] , i ) );
                  is_assigned_to [ i ] = is_assigned_to [ j ]; 
               }
               else
               {
//             i has companion  
                  if ( is_assigned_to [ i ] != is_assigned_to [ j ] )
                  {
//                   move companions to the cluster    
//                 
                     //std::cout <<  "combine " << is_assigned_to [ i ] << " to " << is_assigned_to [ j ] << std::endl;
                     std::multimap< G4int , G4int > clusters_tmp;
                     G4int target_cluster_id; 
                     if ( is_assigned_to [ i ] > is_assigned_to [ j ] )
                        target_cluster_id = is_assigned_to [ i ];
                     else
                        target_cluster_id = is_assigned_to [ j ];
                      
                     for ( std::multimap< G4int , G4int >::iterator it
                         = clusters.begin() ; it != clusters.end() ; it++ ) 
                     {

                        //std::cout << it->first << " " << it->second  << " " << target_cluster_id << std::endl;
                        if ( it->first == target_cluster_id )
                        { 
                           //std::cout << "move " <<  it->first << " " << it->second << std::endl;
                           is_assigned_to [ it->second ] = is_assigned_to [ j ]; 
                           clusters_tmp.insert ( std::multimap<G4int,G4int>::value_type (  is_assigned_to [ j ] , it->second ) );
                        }
                        else
                        {
                           clusters_tmp.insert ( std::multimap<G4int,G4int>::value_type ( it->first , it->second ) );
                        }
                     }

                     clusters = clusters_tmp;
                     //id = clusters.rbegin()->first;
                     //id = target_cluster_id;
                     //std::cout <<  "id " << id << std::endl;
                  }
               }
            } 

            //std::cout << "combination " << i << " " << j << std::endl;
            comb_map.insert( std::multimap<G4int,G4int>::value_type ( i , j ) ); 
            cluster_participants.push_back ( j );



            if ( assign_map.find( cluster_id ) == assign_map.end() ) 
            {  
               is_already_belong_some_cluster[i] = true;
               assign_map.insert ( std::multimap<G4int,G4int>::value_type ( cluster_id , i ) );
               hasThisCompany = true; 
            } 
            assign_map.insert ( std::multimap<G4int,G4int>::value_type ( cluster_id , j ) );  
            is_already_belong_some_cluster[j] = true;

         } 

         if ( ichek == i ) 
         {
            ichek++;
         }
         }

         if ( cluster_participants.size() > 0 ) 
         {
//                                         cluster , participant 
            cluster_map.insert ( std::pair < G4int , std::vector < G4int > > ( i , cluster_participants ) );
         }
      }
//      }
      } 
      if ( hasThisCompany == true ) cluster_id++;
   } 

   //std::cout << " id " << id << std::endl;

// sort
// Heavy cluster comes first
//                size    cluster_id 
   std::multimap< G4int , G4int > sorted_cluster_map; 
   for ( G4int i = 0 ; i <= id ; i++ )  // << "<=" because id is highest cluster nubmer.
   {
 
      //std::cout << i << " cluster has " << clusters.count( i )  << " nucleons." << std::endl;
      sorted_cluster_map.insert ( std::multimap<G4int,G4int>::value_type ( (G4int) clusters.count( i ) , i ) );

   }


// create nucleus from devided clusters
   std::vector < G4QMDNucleus* > result;
   for ( std::multimap < G4int , G4int >::reverse_iterator it 
       = sorted_cluster_map.rbegin() ; it != sorted_cluster_map.rend()  ; it ++) 
   {

      //G4cout << "Add Participants to cluseter " << it->second << G4endl;

      if ( it->first != 0 ) 
      {
         G4QMDNucleus* nucleus = new G4QMDNucleus();
         for ( std::multimap < G4int , G4int >::iterator itt
             = clusters.begin() ; itt != clusters.end()  ; itt ++)
         {

            if ( it->second == itt->first )
            {
               nucleus->SetParticipant( system->GetParticipant ( itt->second ) );
               //G4cout << "Add Participants " << itt->second << " "  << system->GetParticipant ( itt->second )->GetPosition() << G4endl;
             }

         }
         result.push_back( nucleus );
      } 

   }

// delete participants from current system

   for ( std::vector < G4QMDNucleus* > ::iterator it 
       = result.begin() ; it != result.end() ; it++ ) 
   {
      system->SubtractSystem ( *it ); 
   }
   
   return result;
   
}



void G4QMDMeanField::Update() 
{ 
   SetSystem( system ); 
}
