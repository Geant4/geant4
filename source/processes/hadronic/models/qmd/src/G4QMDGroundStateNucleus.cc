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
// 081024 G4NucleiPropertiesTable:: to G4NucleiProperties::
//
#include "G4QMDGroundStateNucleus.hh"

#include "G4NucleiProperties.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Exp.hh"
#include "G4Pow.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

G4QMDGroundStateNucleus::G4QMDGroundStateNucleus( G4int z , G4int a )
: maxTrial ( 1000 )
, r00 ( 1.124 )  // radius parameter for Woods-Saxon [fm] 
, r01 ( 0.5 )    // radius parameter for Woods-Saxon 
, saa ( 0.2 )    // diffuse parameter for initial Woods-Saxon shape
, rada ( 0.9 )   // cutoff parameter
, radb ( 0.3 )   // cutoff parameter
, dsam ( 1.5 )   // minimum distance for same particle [fm]
, ddif ( 1.0 )   // minimum distance for different particle
, edepth ( 0.0 )
, epse ( 0.000001 )  // torelance for energy in [GeV]
, meanfield ( NULL ) 
{

   //std::cout << " G4QMDGroundStateNucleus( G4int z , G4int a ) Begin " << z << " " << a << std::endl;

   dsam2 = dsam*dsam;
   ddif2 = ddif*ddif;

   G4QMDParameters* parameters = G4QMDParameters::GetInstance();

   hbc = parameters->Get_hbc();
   gamm = parameters->Get_gamm();
   cpw = parameters->Get_cpw();
   cph = parameters->Get_cph();
   epsx = parameters->Get_epsx();
   cpc = parameters->Get_cpc();

   cdp = parameters->Get_cdp();
   c0p = parameters->Get_c0p();
   c3p = parameters->Get_c3p();
   csp = parameters->Get_csp();
   clp = parameters->Get_clp();

   //edepth = 0.0; 

   for ( int i = 0 ; i < a ; i++ )
   {

      G4ParticleDefinition* pd; 

      if ( i < z )
      { 
         pd = G4Proton::Proton();
      }
      else
      {
         pd = G4Neutron::Neutron();
      }
         
      G4ThreeVector p( 0.0 );
      G4ThreeVector r( 0.0 );
      G4QMDParticipant* aParticipant = new G4QMDParticipant( pd , p , r );
      SetParticipant( aParticipant );

   }

   G4double radious = r00 * G4Pow::GetInstance()->A13( double ( GetMassNumber() ) ); 

   rt00 = radious - r01; 
   radm = radious - rada * ( gamm - 1.0 ) + radb;
   rmax = 1.0 / ( 1.0 + G4Exp ( -rt00/saa ) );

   //maxTrial = 1000;
   
   //Nucleon primary or target case;
   if ( z == 1 && a == 1 ) {  // Hydrogen  Case or proton primary 
      SetParticipant( new G4QMDParticipant( G4Proton::Proton() , G4ThreeVector( 0.0 ) , G4ThreeVector( 0.0 ) ) );
      ebini = 0.0; 
      return;
   }
   else if ( z == 0 && a == 1 ) { // Neutron primary 
      SetParticipant( new G4QMDParticipant( G4Neutron::Neutron() , G4ThreeVector( 0.0 ) , G4ThreeVector( 0.0 ) ) );
      ebini = 0.0; 
      return;
   }

   
   meanfield = new G4QMDMeanField();
   meanfield->SetSystem( this );

   //std::cout << "G4QMDGroundStateNucleus( G4int z , G4int a ) packNucleons Begin ( z , a ) ( " << z << ", " << a << " )" << std::endl;
   packNucleons();
   //std::cout << "G4QMDGroundStateNucleus( G4int z , G4int a ) packNucleons End" << std::endl;

   delete meanfield;

}



void G4QMDGroundStateNucleus::packNucleons()
{

   //std::cout << "G4QMDGroundStateNucleus::packNucleons" << std::endl;

   ebini = - G4NucleiProperties::GetBindingEnergy( GetMassNumber() , GetAtomicNumber() ) / GetMassNumber();

   G4double ebin00 = ebini * 0.001;

   G4double ebin0 = 0.0;  
   G4double ebin1 = 0.0;  

   if ( GetMassNumber() != 4  )
   {
      ebin0 = ( ebini - 0.5 ) * 0.001;
      ebin1 = ( ebini + 0.5 ) * 0.001;
   }
   else
   {
      ebin0 = ( ebini - 1.5 ) * 0.001;
      ebin1 = ( ebini + 1.5 ) * 0.001;
   }

   G4int n0Try = 0; 
   G4bool isThisOK = false;
   while ( n0Try < maxTrial ) // Loop checking, 11.03.2015, T. Koi
   {
      n0Try++;
      //std::cout << "TKDB packNucleons n0Try " << n0Try << std::endl;

//    Sampling Position

      //std::cout << "TKDB Sampling Position " << std::endl;

      G4bool areThesePsOK = false;
      G4int npTry = 0;
      while ( npTry < maxTrial ) // Loop checking, 11.03.2015, T. Koi
      {
      //std::cout << "TKDB Sampling Position npTry " << npTry << std::endl;
         npTry++; 
         G4int i = 0; 
         if ( samplingPosition( i ) ) 
         {
            //std::cout << "packNucleons samplingPosition 0 succeed " << std::endl;
            for ( i = 1 ; i < GetMassNumber() ; i++ )
            {
               //std::cout << "packNucleons samplingPosition " << i  << " trying " << std::endl;
               if ( !( samplingPosition( i ) ) ) 
               {
                  //std::cout << "packNucleons samplingPosition " << i << " failed" << std::endl;
                  break;
               }
            }
            if ( i == GetMassNumber() ) 
            {
               //std::cout << "packNucleons samplingPosition all scucceed " << std::endl;
               areThesePsOK = true;
               break; 
            }
         }
      }
      if ( areThesePsOK == false ) continue;

      //std::cout << "TKDB Sampling Position End" << std::endl;

//    Calculate Two-body quantities

      meanfield->Cal2BodyQuantities(); 
      std::vector< G4double > rho_a ( GetMassNumber() , 0.0 );
      std::vector< G4double > rho_s ( GetMassNumber() , 0.0 );
      std::vector< G4double > rho_c ( GetMassNumber() , 0.0 );

      rho_l.resize ( GetMassNumber() , 0.0 );
      d_pot.resize ( GetMassNumber() , 0.0 );

      for ( G4int i = 0 ; i < GetMassNumber() ; i++ )
      {
         for ( G4int j = 0 ; j < GetMassNumber() ; j++ )
         {

            rho_a[ i ] += meanfield->GetRHA( j , i ); 
            G4int k = 0; 
            if ( participants[j]->GetDefinition() != participants[i]->GetDefinition() )
            {
               k = 1;
            } 

            rho_s[ i ] += meanfield->GetRHA( j , i )*( 1.0 - 2.0 * k ); // OK?  

            rho_c[ i ] += meanfield->GetRHE( j , i ); 
         }

      }

      for ( G4int i = 0 ; i < GetMassNumber() ; i++ )
      {
         rho_l[i] = cdp * ( rho_a[i] + 1.0 );     
         d_pot[i] = c0p * rho_a[i]
                  + c3p * G4Pow::GetInstance()->powA ( rho_a[i] , gamm )
                  + csp * rho_s[i]
                  + clp * rho_c[i];

         //std::cout << "d_pot[i] " << i << " " << d_pot[i] << std::endl; 
      }


// Sampling Momentum

      //std::cout << "TKDB Sampling Momentum " << std::endl;

      phase_g.clear();
      phase_g.resize( GetMassNumber() , 0.0 );
      
      //std::cout << "TKDB Sampling Momentum 1st " << std::endl;
      G4bool isThis1stMOK = false;
      G4int nmTry = 0; 
      while ( nmTry < maxTrial ) // Loop checking, 11.03.2015, T. Koi
      {
         nmTry++;
         G4int i = 0; 
         if ( samplingMomentum( i ) ) 
         {
           isThis1stMOK = true;
           break; 
         }
      }
      if ( isThis1stMOK == false ) continue;

      //std::cout << "TKDB Sampling Momentum 2nd so on" << std::endl;

      G4bool areTheseMsOK = true;
      nmTry = 0; 
      while ( nmTry < maxTrial ) // Loop checking, 11.03.2015, T. Koi
      {
         nmTry++;
            G4int i = 0; 
            for ( i = 1 ; i < GetMassNumber() ; i++ )
            {
               //std::cout << "TKDB packNucleons samplingMomentum try " << i << std::endl;
               if ( !( samplingMomentum( i ) ) ) 
               {
                  //std::cout << "TKDB packNucleons samplingMomentum " << i << " failed" << std::endl;
                  areTheseMsOK = false;
                  break;
               }
            }
            if ( i == GetMassNumber() ) 
            {
               areTheseMsOK = true;
            }

            break;
      }
      if ( areTheseMsOK == false ) continue;
     
// Kill Angluar Momentum

      //std::cout << "TKDB Sampling Kill Angluar Momentum " << std::endl;

      killCMMotionAndAngularM();    


// Check Binding Energy

      //std::cout << "packNucleons Check Binding Energy Begin " << std::endl;

      G4double ekinal = 0.0;
      for ( int i = 0 ; i < GetMassNumber() ; i++ )
      {
         ekinal += participants[i]->GetKineticEnergy();
      }

      meanfield->Cal2BodyQuantities();

      G4double totalPotentialE = meanfield->GetTotalPotential(); 
      G4double ebinal = ( totalPotentialE + ekinal ) / double ( GetMassNumber() );

      //std::cout << "packNucleons totalPotentialE " << totalPotentialE << std::endl;
      //std::cout << "packNucleons ebinal " << ebinal << std::endl;
      //std::cout << "packNucleons ekinal " << ekinal << std::endl;

      if ( ebinal < ebin0 || ebinal > ebin1 ) 
      {
         //std::cout << "packNucleons ebin0 " << ebin0 << std::endl;
         //std::cout << "packNucleons ebin1 " << ebin1 << std::endl;
         //std::cout << "packNucleons ebinal " << ebinal << std::endl;
      //std::cout << "packNucleons Check Binding Energy Failed " << std::endl;
         continue;
      }

      //std::cout << "packNucleons Check Binding Energy End = OK " << std::endl;


// Energy Adujstment

      G4double dtc = 1.0;
      G4double frg = -0.1;
      G4double rdf0 = 0.5;
      
      G4double edif0 = ebinal - ebin00;

      G4double cfrc = 0.0;
      if ( 0 < edif0 ) 
      {
         cfrc = frg;
      }
      else
      {
         cfrc = -frg;
      }

      G4int ifrc = 1;

      G4int neaTry = 0;

      G4bool isThisEAOK = false;
      while ( neaTry < maxTrial )  // Loop checking, 11.03.2015, T. Koi
      {
         neaTry++;

         G4double  edif = ebinal - ebin00; 

         //std::cout << "TKDB edif " << edif << std::endl; 
         if ( std::abs ( edif ) < epse )
         {
   
            isThisEAOK = true;
            //std::cout << "isThisEAOK " << isThisEAOK << std::endl; 
            break;
         }

         G4int jfrc = 0;
         if ( edif < 0.0 ) 
         {
            jfrc = 1;
         }
         else
         {
            jfrc = -1;
         }

         if ( jfrc != ifrc ) 
         {
            cfrc = -rdf0 * cfrc;
            dtc = rdf0 * dtc;
         }

         if ( jfrc == ifrc && std::abs( edif0 ) < std::abs( edif ) )
         {
            cfrc = -rdf0 * cfrc;
            dtc = rdf0 * dtc;
         }

         ifrc = jfrc;
         edif0 = edif;

         meanfield->CalGraduate();

         for ( int i = 0 ; i < GetMassNumber() ; i++ )
         {
            G4ThreeVector ri = participants[i]->GetPosition(); 
            G4ThreeVector p_i = participants[i]->GetMomentum(); 

            ri += dtc * ( meanfield->GetFFr(i) - cfrc * ( meanfield->GetFFp(i) ) );
            p_i += dtc * ( meanfield->GetFFp(i) + cfrc * ( meanfield->GetFFr(i) ) );

            participants[i]->SetPosition( ri ); 
            participants[i]->SetMomentum( p_i ); 
         }

         ekinal = 0.0;     

         for ( int i = 0 ; i < GetMassNumber() ; i++ )
         {
            ekinal += participants[i]->GetKineticEnergy(); 
         }

         meanfield->Cal2BodyQuantities(); 
         totalPotentialE = meanfield->GetTotalPotential(); 

         ebinal = ( totalPotentialE + ekinal ) / double ( GetMassNumber() );

      }
      //std::cout << "isThisEAOK " << isThisEAOK << std::endl; 
      if ( isThisEAOK == false ) continue;
   
      isThisOK = true;
      //std::cout << "isThisOK " << isThisOK << std::endl; 
      break; 

   }

   if ( isThisOK == false )
   {
      G4cout << "GroundStateNucleus state cannot be created. Try again with another parameters." << G4endl;
   } 

   //std::cout << "packNucleons End" << std::endl;
   return;
}


G4bool G4QMDGroundStateNucleus::samplingPosition( G4int i )
{

   G4bool result = false;

   G4int nTry = 0; 
   
   while ( nTry < maxTrial )  // Loop checking, 11.03.2015, T. Koi
   {

      //std::cout << "samplingPosition i th particle, nTtry " << i << " " << nTry << std::endl;  

      G4double rwod = -1.0;        
      G4double rrr = 0.0; 

      G4double rx = 0.0;
      G4double ry = 0.0;
      G4double rz = 0.0;

      G4int icounter = 0;
      G4int icounter_max = 1024;
      while ( G4UniformRand() * rmax > rwod ) // Loop checking, 11.03.2015, T. Koi
      {
         icounter++;
         if ( icounter > icounter_max ) {
	    G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
            break;
         }

         G4double rsqr = 10.0; 
         G4int jcounter = 0;
         G4int jcounter_max = 1024;
         while ( rsqr > 1.0 ) // Loop checking, 11.03.2015, T. Koi
         {
            jcounter++;
            if ( jcounter > jcounter_max ) {
	       G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
               break;
            }
            rx = 1.0 - 2.0 * G4UniformRand();
            ry = 1.0 - 2.0 * G4UniformRand();
            rz = 1.0 - 2.0 * G4UniformRand();
            rsqr = rx*rx + ry*ry + rz*rz; 
         }
         rrr = radm * std::sqrt ( rsqr );
         rwod = 1.0 / ( 1.0 + G4Exp ( ( rrr - rt00 ) / saa ) );

      } 

      participants[i]->SetPosition( G4ThreeVector( rx , ry , rz )*radm );

      if ( i == 0 )   
      {
         result = true; 
         return result;
      }

//    i > 1 ( Second Particle or later )
//    Check Distance to others 

      G4bool isThisOK = true;
      for ( G4int j = 0 ; j < i ; j++ )
      {

         G4double r2 =  participants[j]->GetPosition().diff2( participants[i]->GetPosition() );   
         G4double dmin2 = 0.0;

         if ( participants[j]->GetDefinition() == participants[i]->GetDefinition() )
         {
            dmin2 = dsam2;
         }
         else
         {
            dmin2 = ddif2;
         }

         //std::cout << "distance between j and i " << j << " " << i << " " << r2 << " " << dmin2 << std::endl;
         if ( r2 < dmin2 )
         {
            isThisOK = false;
            break; 
         }

      }

      if ( isThisOK == true )
      {
         result = true;
         return result; 
      }

      nTry++; 

   }

// Here return "false" 
   return result;
}



G4bool G4QMDGroundStateNucleus::samplingMomentum( G4int i )
{

   //std::cout << "TKDB samplingMomentum for " << i << std::endl;
   
   G4bool result = false;

   G4double pfm = hbc * G4Pow::GetInstance()->A13 ( ( 3.0 / 2.0 * pi*pi * rho_l[i] ) );

   if ( 10 < GetMassNumber() &&  -5.5 < ebini ) 
   {
      pfm = pfm * ( 1.0 + 0.2 * std::sqrt( std::abs( 8.0 + ebini ) / 8.0 ) );
   }

   //std::cout << "TKDB samplingMomentum pfm " << pfm << std::endl;

   std::vector< G4double > phase; 
   phase.resize( i+1 ); // i start from 0

   G4int ntry = 0;
// 710 
   while ( ntry < maxTrial )  // Loop checking, 11.03.2015, T. Koi
   {

      //std::cout << " TKDB ntry " << ntry << std::endl;
      ntry++;

      G4double ke = DBL_MAX;

      G4int tkdb_i =0;
// 700
      G4int icounter = 0;
      G4int icounter_max = 1024;
      while ( ke + d_pot [i] > edepth ) // Loop checking, 11.03.2015, T. Koi
      {
         icounter++;
         if ( icounter > icounter_max ) {
	    G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
            break;
         }
      
         G4double psqr = 10.0;
         G4double px = 0.0;
         G4double py = 0.0;
         G4double pz = 0.0;

         G4int jcounter = 0;
         G4int jcounter_max = 1024;
         while ( psqr > 1.0 ) // Loop checking, 11.03.2015, T. Koi
         {
            jcounter++;
            if ( jcounter > jcounter_max ) {
	       G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
               break;
            }
            px = 1.0 - 2.0*G4UniformRand();
            py = 1.0 - 2.0*G4UniformRand();
            pz = 1.0 - 2.0*G4UniformRand();

            psqr = px*px + py*py + pz*pz;
         }

         G4ThreeVector p ( px , py , pz ); 
         p = pfm * p;
         participants[i]->SetMomentum( p );
         G4LorentzVector p4 = participants[i]->Get4Momentum();
         //ke = p4.e() - p4.restMass();
         ke = participants[i]->GetKineticEnergy();
   

         tkdb_i++;  
         if ( tkdb_i > maxTrial ) return result; // return false

      }

      //std::cout << "TKDB ke d_pot[i] " << ke << " " << d_pot[i] << std::endl;


      if ( i == 0 )   
      {
         result = true; 
         return result;
      }

      G4bool isThisOK = true;

      // Check Pauli principle

      phase[ i ] = 0.0; 

      //std::cout << "TKDB Check Pauli Principle " << i << std::endl;

      for ( G4int j = 0 ; j < i ; j++ )
      {
         phase[ j ] = 0.0;
         //std::cout << "TKDB Check Pauli Principle  i , j " << i << " , " << j << std::endl;
         G4double expa = 0.0;
         if ( participants[j]->GetDefinition() ==  participants[i]->GetDefinition() )
         {

            expa = - meanfield->GetRR2(i,j) * cpw;

            if ( expa > epsx ) 
            {
               G4ThreeVector p_i = participants[i]->GetMomentum();  
               G4ThreeVector pj = participants[j]->GetMomentum();  
               G4double dist2_p = p_i.diff2( pj ); 

               dist2_p = dist2_p*cph;
               expa = expa - dist2_p; 

               if ( expa > epsx ) 
               {

                  phase[j] = G4Exp ( expa );

                  if ( phase[j] * cpc > 0.2 ) 
                  { 
/*
         G4cout << "TKDB Check Pauli Principle A i , j " << i << " , " << j << G4endl;
         G4cout << "TKDB Check Pauli Principle phase[j] " << phase[j] << G4endl;
         G4cout << "TKDB Check Pauli Principle phase[j]*cpc > 0.2 " << phase[j]*cpc << G4endl;
*/
                     isThisOK = false;
                     break;
                  }
                  if ( ( phase_g[j] + phase[j] ) * cpc > 0.5 ) 
                  { 
/*
         G4cout << "TKDB Check Pauli Principle B i , j " << i << " , " << j << G4endl;
         G4cout << "TKDB Check Pauli Principle B phase_g[j] " << phase_g[j] << G4endl;
         G4cout << "TKDB Check Pauli Principle B phase[j] " << phase[j] << G4endl;
         G4cout << "TKDB Check Pauli Principle B phase_g[j] + phase[j] ) * cpc > 0.5  " <<  ( phase_g[j] + phase[j] ) * cpc << G4endl;
*/
                     isThisOK = false;
                     break;
                  }

                  phase[i] += phase[j];
                  if ( phase[i] * cpc > 0.3 ) 
                  { 
/*
         G4cout << "TKDB Check Pauli Principle C i , j " << i << " , " << j << G4endl;
         G4cout << "TKDB Check Pauli Principle C phase[i] " << phase[i] << G4endl;
         G4cout << "TKDB Check Pauli Principle C phase[i] * cpc > 0.3 " <<  phase[i] * cpc << G4endl;
*/
                     isThisOK = false;
                     break;
                  }

                  //std::cout << "TKDB Check Pauli Principle OK i , j " << i << " , " << j << std::endl;

               }
               else
               {
                  //std::cout << "TKDB Check Pauli Principle OK i , j " << i << " , " << j << std::endl;
               }

            }
            else
            {
               //std::cout << "TKDB Check Pauli Principle OK i , j " << i << " , " << j << std::endl;
            }

         }
         else
         {
            //std::cout << "TKDB Check Pauli Principle OK i , j " << i << " , " << j << std::endl;
         }

      }

      if ( isThisOK == true )
      {

         phase_g[i] = phase[i];

         for ( int j = 0 ; j < i ; j++ )
         {
            phase_g[j] += phase[j]; 
         }

         result = true; 
         return result;
      }

   }

   return result;

}



void G4QMDGroundStateNucleus::killCMMotionAndAngularM()
{

//   CalEnergyAndAngularMomentumInCM();

   //std::vector< G4ThreeVector > p ( GetMassNumber() , 0.0 );
   //std::vector< G4ThreeVector > r ( GetMassNumber() , 0.0 );

// Move to cm system

   G4ThreeVector pcm_tmp ( 0.0 );
   G4ThreeVector rcm_tmp ( 0.0 );
   G4double sumMass = 0.0;

   for ( G4int i = 0 ; i < GetMassNumber() ; i++ )
   {
      pcm_tmp += participants[i]->GetMomentum();
      rcm_tmp += participants[i]->GetPosition() * participants[i]->GetMass();
      sumMass += participants[i]->GetMass();
   }

   pcm_tmp = pcm_tmp/GetMassNumber();
   rcm_tmp = rcm_tmp/sumMass;

   for ( G4int i = 0 ; i < GetMassNumber() ; i++ )
   {
      participants[i]->SetMomentum( participants[i]->GetMomentum() - pcm_tmp );
      participants[i]->SetPosition( participants[i]->GetPosition() - rcm_tmp );
   }

// kill the angular momentum

   G4ThreeVector ll ( 0.0 );
   for ( G4int i = 0 ; i < GetMassNumber() ; i++ )
   {
      ll += participants[i]->GetPosition().cross ( participants[i]->GetMomentum() );
   }

   G4double rr[3][3];
   G4double ss[3][3];
   for ( G4int i = 0 ; i < 3 ; i++ )
   {
      for ( G4int j = 0 ; j < 3 ; j++ )
      {
         rr [i][j] = 0.0;

         if ( i == j ) 
         {
            ss [i][j] = 1.0;
         }
         else
         {
            ss [i][j] = 0.0;
         }
      } 
   }

   for ( G4int i = 0 ; i < GetMassNumber() ; i++ )
   {
      G4ThreeVector r = participants[i]->GetPosition();
      rr[0][0] += r.y() * r.y() + r.z() * r.z(); 
      rr[1][0] += - r.y() * r.x();
      rr[2][0] += - r.z() * r.x();
      rr[0][1] += - r.x() * r.y();
      rr[1][1] += r.z() * r.z() + r.x() * r.x(); 
      rr[2][1] += - r.z() * r.y();
      rr[2][0] += - r.x() * r.z();
      rr[2][1] += - r.y() * r.z();
      rr[2][2] += r.x() * r.x() + r.y() * r.y(); 
   }

   for ( G4int i = 0 ; i < 3 ; i++ )
   {
      G4double x = rr [i][i];
      for ( G4int j = 0 ; j < 3 ; j++ )
      {
         rr[i][j] = rr[i][j] / x;
         ss[i][j] = ss[i][j] / x;
      }
      for ( G4int j = 0 ; j < 3 ; j++ )
      {
         if ( j != i ) 
         {
            G4double y = rr [j][i];
            for ( G4int k = 0 ; k < 3 ; k++ )
            {
               rr[j][k] += -y * rr[i][k];
               ss[j][k] += -y * ss[i][k];
            }
         }
      }
   }

   G4double opl[3];
   G4double rll[3];

   rll[0] = ll.x();
   rll[1] = ll.y();
   rll[2] = ll.z();
   
   for ( G4int i = 0 ; i < 3 ; i++ )
   {
      opl[i] = 0.0;

      for ( G4int j = 0 ; j < 3 ; j++ )
      {
	 opl[i] += ss[i][j]*rll[j];
      }
   }

   for ( G4int i = 0 ; i < GetMassNumber() ; i++ )
   {
      G4ThreeVector p_i = participants[i]->GetMomentum() ;
      G4ThreeVector ri = participants[i]->GetPosition() ;
      G4ThreeVector opl_v ( opl[0] , opl[1] , opl[2] );  

      p_i += ri.cross(opl_v);

      participants[i]->SetMomentum( p_i );
   }

}
