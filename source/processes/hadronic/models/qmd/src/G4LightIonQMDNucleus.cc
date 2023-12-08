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
// 230309 Skyrme-QMD parameters added by Y-H. Sato and A. Haga
// 230309 Total energy evaluated by Lorentz covariant version by Y-H. Sato and A. Haga

#include <numeric>

#include "G4LightIonQMDNucleus.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4NucleiProperties.hh"
#include "G4HadronicException.hh"

#include "G4LightIonQMDParameters.hh"    // 20230309
#include "G4PhysicalConstants.hh"   // 20230309
#include <cmath>                    // 20230309
#include <CLHEP/Random/Stat.h>      // 20230309

G4LightIonQMDNucleus::G4LightIonQMDNucleus()
{
   G4LightIonQMDParameters* parameters = G4LightIonQMDParameters::GetInstance();
   hbc = parameters->Get_hbc();
   
   jj = 0; // will be calcualted in CalEnergyAndAngularMomentumInCM;
   potentialEnergy = 0.0; // will be set through set method 
   excitationEnergy = 0.0;

   // Following Parameters are added (20230309)
   wl = parameters->Get_wl();
   cl = parameters->Get_cl();
   rho0 = parameters->Get_rho0();
   gamm = parameters->Get_gamm();
   eta = parameters->Get_eta(); // Skyrme-QMD
   kappas = parameters->Get_kappas(); // Skyrme-QMD
    
   cpw = parameters->Get_cpw();
   cph = parameters->Get_cph();
   cpc = parameters->Get_cpc();
    
   c0 = parameters->Get_c0();
   c3 = parameters->Get_c3();
   cs = parameters->Get_cs();
   g0 = parameters->Get_g0(); // Skyrme-QMD
   g0iso = parameters->Get_g0iso(); // Skyrme-QMD
   gtau0 = parameters->Get_gtau0(); // Skyrme-QMD
    
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
   pag_tau = eta - 1; // Skyrme-QMD
   cg0 = - g0 / ( 2.0 * wl );  // Skyrme-QMD
   cgtau0 = - gtau0 / ( 4.0 * wl ) * eta;  // Skyrme-QMD
}



//G4LightIonQMDNucleus::~G4LightIonQMDNucleus()
//{
//   ;
//}


G4LorentzVector G4LightIonQMDNucleus::Get4Momentum()
{
   G4LorentzVector p( 0 );
   std::vector< G4QMDParticipant* >::iterator it; 
   for ( it = participants.begin() ; it != participants.end() ; it++ ) 
      p += (*it)->Get4Momentum();   

   return p;
}



G4int G4LightIonQMDNucleus::GetMassNumber()
{

   G4int A = 0; 
   std::vector< G4QMDParticipant* >::iterator it; 
   for ( it = participants.begin() ; it != participants.end() ; it++ ) 
   {
      if ( (*it)->GetDefinition() == G4Proton::Proton() 
        || (*it)->GetDefinition() == G4Neutron::Neutron() ) 
         A++; 
   }

   if ( A == 0 ) {
      throw G4HadronicException(__FILE__, __LINE__, "G4LightIonQMDNucleus has the mass number of 0!");
   }

   return A;
}



G4int G4LightIonQMDNucleus::GetAtomicNumber()
{
   G4int Z = 0; 
   std::vector< G4QMDParticipant* >::iterator it; 
   for ( it = participants.begin() ; it != participants.end() ; it++ ) 
   {
      if ( (*it)->GetDefinition() == G4Proton::Proton() ) 
         Z++; 
   }
   return Z;
}



G4double G4LightIonQMDNucleus::GetNuclearMass()
{

   G4double mass = G4NucleiProperties::GetNuclearMass( GetMassNumber() , GetAtomicNumber() );
   
   if ( mass == 0.0 )
   {

      G4int Z = GetAtomicNumber();
      G4int A = GetMassNumber();
      G4int N = A - Z;

// Weizsacker-Bethe 

      G4double Av = 16*MeV; 
      G4double As = 17*MeV; 
      G4double Ac = 0.7*MeV; 
      G4double Asym = 23*MeV; 

      G4double BE = Av * A 
                  - As * G4Pow::GetInstance()->A23 ( G4double ( A ) ) 
                  - Ac * Z*Z/G4Pow::GetInstance()->A13 ( G4double ( A ) )
                  - Asym * ( N - Z )* ( N - Z ) / A; 

      mass = Z * G4Proton::Proton()->GetPDGMass() 
           + N * G4Neutron::Neutron()->GetPDGMass()
           - BE;

   }

   return mass;
}



void G4LightIonQMDNucleus::CalEnergyAndAngularMomentumInCM()
{

   //G4cout << "CalEnergyAndAngularMomentumInCM " << this->GetAtomicNumber() << " " << GetMassNumber() << G4endl;

   G4double gamma = Get4Momentum().gamma();
   G4ThreeVector beta = Get4Momentum().v()/ Get4Momentum().e();

   G4ThreeVector pcm0( 0.0 ) ;

   G4int n = GetTotalNumberOfParticipant();
   pcm.resize( n );

   for ( G4int i= 0; i < n ; i++ ) 
   {
      G4ThreeVector p_i = GetParticipant( i )->GetMomentum();

      G4double trans = gamma / ( gamma + 1.0 ) * p_i * beta; 
      pcm[i] = p_i - trans*beta;

      pcm0 += pcm[i];
   }

   pcm0 = pcm0 / double ( n );

   //G4cout << "pcm0 " << pcm0 << G4endl;

   for ( G4int i= 0; i < n ; i++ ) 
   {
      pcm[i] += -pcm0;
      //G4cout << "pcm " << i << " " << pcm[i] << G4endl;
   }


   G4double tmass = 0;
   G4ThreeVector rcm0( 0.0 ) ;
   rcm.resize( n );
   es.resize( n );

   // binding energy should be evaluated with a relativistic version: 20230308 by Y-H. Sato and A. Haga
   for ( G4int i= 0; i < n ; i++ )
   {
      G4ThreeVector ri = GetParticipant( i )->GetPosition();
      G4double trans = gamma / ( gamma + 1.0 ) * ri * beta;
      G4double nucpote = GetNuclPotential( i );

      es[i] = std::sqrt ( G4Pow::GetInstance()->powN ( GetParticipant( i )->GetMass() , 2 ) + pcm[i]*pcm[i] + 2.0*GetParticipant( i )->GetMass()*nucpote) - GetParticipant( i )->GetMass(); //R-JQMD

      rcm[i] = ri + trans*beta;

      rcm0 += rcm[i]*es[i];

      tmass += es[i];
   }

   rcm0 = rcm0/tmass;

   for ( G4int i= 0; i < n ; i++ ) 
   {
      rcm[i] += -rcm0;
      //G4cout << "rcm " << i << " " << rcm[i] << G4endl;
   }

// Angular momentum

   G4ThreeVector rl ( 0.0 ); 
   for ( G4int i= 0; i < n ; i++ ) 
   {
      rl += rcm[i].cross ( pcm[i] );
   }

// DHW: move hbc outside of sqrt to get correct units
//  jj = int ( std::sqrt ( rl*rl / hbc ) + 0.5 );

   jj = int (std::sqrt(rl*rl)/hbc + 0.5);

// kinetic energy per nucleon in CM

    /*
   G4double totalMass = 0.0;
   for ( G4int i= 0; i < n ; i++ ) 
   {
      // following two lines are equivalent
      //totalMass += GetParticipant( i )->GetDefinition()->GetPDGMass()/GeV;
      totalMass += GetParticipant( i )->GetMass();
   }
   */

   //G4double kineticEnergyPerNucleon = ( std::accumulate ( es.begin() , es.end() , 0.0 ) - totalMass )/n;

// Total (not per nucleion ) Binding Energy
   // relativistic version Y-H. Sato and A. Haga 20230309
   G4double bindingEnergy =  ( std::accumulate ( es.begin() , es.end() , 0.0 ) );

   //G4cout << "n " << n << "totalpote " << totalpote << " " << potentialEnergy << " " << bindingEnergy << G4endl;
   //G4cout << "KineticEnergyPerNucleon in GeV " << kineticEnergyPerNucleon << G4endl;
   //G4cout << "KineticEnergySum in GeV " << std::accumulate ( es.begin() , es.end() , 0.0 ) - totalMass << G4endl;
   //G4cout << "PotentialEnergy in GeV " << potentialEnergy << G4endl;
   //G4cout << "BindingEnergy in GeV " << bindingEnergy << G4endl;
   //G4cout << "G4BindingEnergy in GeV " << G4NucleiProperties::GetBindingEnergy( GetAtomicNumber() , GetMassNumber() )/GeV << G4endl;

   excitationEnergy = bindingEnergy + G4NucleiProperties::GetBindingEnergy( GetMassNumber() , GetAtomicNumber() )/GeV;
   if ( excitationEnergy < 0 ) excitationEnergy = 0.0;

 }

// Get potential with a relativistic version added by Y-H. Sato and A. Haga 20230309
G4double G4LightIonQMDNucleus::GetNuclPotential( G4int i )
{
    epsx = -20.0;
    epscl = 0.0001; // coulomb term
    irelcr = 1;
    G4int n = GetTotalNumberOfParticipant();

    G4double rhoa = 0.0;
    G4double rho3 = 0.0;
    G4double fsij_rhoa = 0.0; // Skyrme-QMD
    //    G4double fsij_rhos = 0.0; // Skyrme-QMD
    G4double rho3_tau = 0.0; // Skyrme-QMD
    G4double rhos = 0.0;
    G4double rhoc = 0.0;
    
    
    G4int icharge = GetParticipant(i)->GetChargeInUnitOfEplus();
    G4int inuc = GetParticipant(i)->GetNuc();
    G4int ibry = GetParticipant(i)->GetBaryonNumber();
    
    G4ThreeVector ri = GetParticipant( i )->GetPosition();
    G4LorentzVector p4i = GetParticipant( i )->Get4Momentum();

    for ( G4int j = 0 ; j < n ; j ++ )
    {
        G4double cef = 1.0;
        if (i == j)
        {
            cef = 0.0;
        }

        
        G4int jcharge = GetParticipant(j)->GetChargeInUnitOfEplus();
        G4int jnuc = GetParticipant(j)->GetNuc();
        G4int jbry = GetParticipant(j)->GetBaryonNumber();
        
        G4ThreeVector rj = GetParticipant( j )->GetPosition();
        G4LorentzVector p4j = GetParticipant( j )->Get4Momentum();
        
        G4ThreeVector rij = ri - rj;
        G4ThreeVector pij = (p4i - p4j).v();
        G4LorentzVector p4ij = p4i - p4j;
        G4ThreeVector bij = ( p4i + p4j ).boostVector();
        G4double gammaij = ( p4i + p4j ).gamma();
        
        //G4double eij = ( p4i + p4j ).e();
        
        G4double rbrb = rij*bij;
        //         G4double bij2 = bij*bij;
        G4double rij2 = rij*rij;
        //G4double pij2 = pij*pij;
        
        
        rbrb = irelcr * rbrb;
        G4double  gamma2_ij = gammaij*gammaij;
        
        G4double rr2 = rij2 + gamma2_ij * rbrb*rbrb;
        
        G4double expa1 = - (rij2 + gamma2_ij * rbrb*rbrb) * c0w;
        G4double rh1;
        if ( expa1 > epsx )
        {
            rh1 = G4Exp( expa1 );
        }
        else
        {
            rh1 = 0.0;
        }
        
        G4double rrs2 = (rij2 + gamma2_ij * rbrb*rbrb) + epscl;
        G4double rrs = std::sqrt ( rrs2 );
        
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
            xerf = std::erf ( rrs*c0sw );
#endif
        } else {
            xerf = 1.0;
        }
        
        G4double erfij = xerf/rrs;
        
        G4double fsij = 3.0/(2*wl) - rr2/(2*wl)/(2*wl); // Add for Skyrme-QMD
        
        rhoa += ibry*jbry*rh1*cef;
        fsij_rhoa += fsij * ibry*jbry*rh1*cef; // Skyrme-QMD
        rhoc += icharge*jcharge * erfij * cef;
        rhos += ibry*jbry*rh1 * jnuc * inuc * cef
        * ( 1 - 2 * std::abs ( jcharge - icharge ) )
        * (1. - kappas * fsij);
        
        //G4cout << i << " " << j << " " << ( - erfij ) << " " << clw << G4endl;
        

    }
    
    rho3 = G4Pow::GetInstance()->powA ( rhoa , gamm );
    rho3_tau = G4Pow::GetInstance()->powA ( rhoa , eta );
    
    G4double potential = c0 * rhoa
    + c3 * rho3
    + g0 * fsij_rhoa  // Skyrme-QMD
    // + g0iso * fsij_rhos  // Skyrme-QMD
    + gtau0 * rho3_tau  // Skyrme-QMD
    + cs * rhos
    + cl * rhoc;
    
    //G4cout << "n " << n << " " << rho3 << G4endl;
    return potential;
}
