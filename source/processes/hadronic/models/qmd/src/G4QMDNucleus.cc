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
#include <numeric>

#include "G4QMDNucleus.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4NucleiProperties.hh"
#include "G4HadronicException.hh"

G4QMDNucleus::G4QMDNucleus()
{
   G4QMDParameters* parameters = G4QMDParameters::GetInstance();
   hbc = parameters->Get_hbc();
   
   jj = 0; // will be calcualted in CalEnergyAndAngularMomentumInCM;
   potentialEnergy = 0.0; // will be set through set method 
   excitationEnergy = 0.0;
}



//G4QMDNucleus::~G4QMDNucleus()
//{
//   ;
//}


G4LorentzVector G4QMDNucleus::Get4Momentum()
{
   G4LorentzVector p( 0 );
   std::vector< G4QMDParticipant* >::iterator it; 
   for ( it = participants.begin() ; it != participants.end() ; it++ ) 
      p += (*it)->Get4Momentum();   

   return p;
}



G4int G4QMDNucleus::GetMassNumber()
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
      throw G4HadronicException(__FILE__, __LINE__, "G4QMDNucleus has the mass number of 0!");
   }

   return A;
}



G4int G4QMDNucleus::GetAtomicNumber()
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



G4double G4QMDNucleus::GetNuclearMass()
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



void G4QMDNucleus::CalEnergyAndAngularMomentumInCM()
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

   for ( G4int i= 0; i < n ; i++ ) 
   {
      G4ThreeVector ri = GetParticipant( i )->GetPosition();
      G4double trans = gamma / ( gamma + 1.0 ) * ri * beta; 

      es[i] = std::sqrt ( G4Pow::GetInstance()->powN ( GetParticipant( i )->GetMass() , 2 ) + pcm[i]*pcm[i] );

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

// Angluar momentum

   G4ThreeVector rl ( 0.0 ); 
   for ( G4int i= 0; i < n ; i++ ) 
   {
      rl += rcm[i].cross ( pcm[i] );
   }

   jj = int ( std::sqrt ( rl*rl / hbc ) + 0.5 );


// kinetic energy per nucleon in CM

   G4double totalMass = 0.0;
   for ( G4int i= 0; i < n ; i++ ) 
   {
      // following two lines are equivalent
      //totalMass += GetParticipant( i )->GetDefinition()->GetPDGMass()/GeV;
      totalMass += GetParticipant( i )->GetMass();
   }

   //G4double kineticEnergyPerNucleon = ( std::accumulate ( es.begin() , es.end() , 0.0 ) - totalMass )/n;

// Total (not per nucleion ) Binding Energy 
   G4double bindingEnergy =  ( std::accumulate ( es.begin() , es.end() , 0.0 ) -totalMass ) + potentialEnergy;

   //G4cout << "KineticEnergyPerNucleon in GeV " << kineticEnergyPerNucleon << G4endl;
   //G4cout << "KineticEnergySum in GeV " << std::accumulate ( es.begin() , es.end() , 0.0 ) - totalMass << G4endl;
   //G4cout << "PotentialEnergy in GeV " << potentialEnergy << G4endl;
   //G4cout << "BindingEnergy in GeV " << bindingEnergy << G4endl;
   //G4cout << "G4BindingEnergy in GeV " << G4NucleiProperties::GetBindingEnergy( GetAtomicNumber() , GetMassNumber() )/GeV << G4endl;

   excitationEnergy = bindingEnergy + G4NucleiProperties::GetBindingEnergy( GetMassNumber() , GetAtomicNumber() )/GeV;
   //G4cout << "excitationEnergy in GeV " << excitationEnergy << G4endl;
   if ( excitationEnergy < 0 ) excitationEnergy = 0.0; 

}
