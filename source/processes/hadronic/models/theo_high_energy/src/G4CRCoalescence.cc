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
//
//---------------------------------------------------------------------------
//
// ClassName:    G4CRCoalescence   ("CR" stands for "Cosmic Ray")
//
// Author:       2020 Alberto Ribon , based on code written by
//               Diego Mauricio Gomez Coral for the GAPS Collaboration
//
// Description:  This class can be optionally used in the method:
//
//                 G4TheoFSGenerator::ApplyYourself
//
//               to coalesce pairs of proton-neutron and antiproton-antineutron
//               into deuterons and antideuterons, respectively, from the list
//               of secondaries produced by a string model.
//               This class can be useful in particular for Cosmic Ray (CR)
//               applications.
//               By default, this class is not used.
//               However, it can be enabled via the UI command:
//
//                 /process/had/enableCRCoalescence true
//
//               It is assumed that the candidate proton-neutron and
//               antiproton-antideuteron pairs originate from the same
//               spatial position, so the condition for coalescence takes
//               into account only their closeness in momentum space.
//
//               This class is based entirely on code written by
//               Diego Mauricio Gomez Coral for the GAPS Collaboration.
//               The main application of this work is for cosmic ray physics.
//
//               Notes:
//               -  In its current version, coalescence can occur only for
//                  proton projectile (because the coalescence parameters
//                  for deuteron and antideuteron are set to non-null values
//                  only for the case of proton projectile).
//               -  This class is not meant be used for secondaries produces
//                  by intranuclear cascade models - such as BERT, BIC and
//                  INCL - which should have already a coalescence phase.
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include "G4DynamicParticle.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Deuteron.hh"
#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"
#include "G4CRCoalescence.hh"
#include "G4ReactionProduct.hh"
#include "G4IonTable.hh"
#include "G4PhysicsModelCatalog.hh"


G4CRCoalescence::G4CRCoalescence() : G4HadronicInteraction("G4CRCoalescence" ),
  fP0_d( 0.0 ), fP0_dbar( 0.0 ), secID( -1 ) {
  secID = G4PhysicsModelCatalog::GetModelID( "model_G4CRCoalescence" );
}


G4CRCoalescence::~G4CRCoalescence() {}


void G4CRCoalescence::SetP0Coalescence( const G4HadProjectile &thePrimary, G4String /* model */ ) {
  // Note by A.R. : in the present version, the coalescence parameters are set only for
  //                proton projectile. If we want to extend this coalescence algorithm
  //                for other applications, besides cosmic rays, we need to set these
  //                coalescence parameters also for all projectiles.
  //                (Note that the method "GenerateDeuterons", instead, can be already used
  //                 as it is for all projectiles.)
  fP0_dbar = 0.0;
  fP0_d    = 0.0;
  if ( thePrimary.GetDefinition()->GetPDGEncoding() == 2212 ) {  // proton
    G4double mproj = thePrimary.GetDefinition()->GetPDGMass();
    G4double pz = thePrimary.Get4Momentum().z();
    G4double ekin = std::sqrt( pz*pz + mproj*mproj ) - mproj;
    if ( ekin > 10.0 ) {
      fP0_dbar = 130.0 / ( 1.0 + std::exp( 21.6 - std::log( 0.001*ekin )/0.089 ) );  // set p0 for antideuteron 
      fP0_d    = 118.1 * ( 1.0 + std::exp( 5.53 - std::log( 0.001*ekin )/0.43 ) );   // set p0 for deuteron
    }
  }
  //G4cout << "Coalescence parameter p0 deuteron / antideuteron: " << fP0_d << " / " << fP0_dbar << G4endl;
}


void G4CRCoalescence::GenerateDeuterons( G4ReactionProductVector* result ) {
  // Deuteron clusters are made with the first nucleon pair that fulfills
  // the coalescence conditions, starting with the protons.
  // A deuteron is a pair (i,j) where i is the proton and j the neutron in current event
  // with the relative momentum less than p0 (i.e. within a sphere of radius p0).
  // The same applies for antideuteron clusters, with antiprotons and antineutrons,
  // instead of protons and neutrons, respectively.
  
  // Vectors of index-position and 3-momentum pairs for, respectively:
  // protons, neutrons, antiprotons and antineutrons
  std::vector< std::pair< G4int, G4ThreeVector > > proton;
  std::vector< std::pair< G4int, G4ThreeVector > > neutron;
  std::vector< std::pair< G4int, G4ThreeVector > > antiproton;
  std::vector< std::pair< G4int, G4ThreeVector > > antineutron;
  for ( unsigned int i = 0; i < result->size(); ++i ) {
    G4int pdgid = result->operator[](i)->GetDefinition()->GetPDGEncoding();
    if ( pdgid == 2212 ) {  // proton
      proton.push_back( std::make_pair( i, result->operator[](i)->GetMomentum() ) );
      result->erase( result->begin() + i );
    }
  }
  for ( unsigned int i = 0; i < result->size(); ++i ) {
    G4int pdgid = result->operator[](i)->GetDefinition()->GetPDGEncoding();
    if ( pdgid == 2112 ) {  // neutron
      neutron.push_back( std::make_pair( i, result->operator[](i)->GetMomentum() ) );
      result->erase( result->begin() + i );
    }
  }
  for ( unsigned int i = 0; i < result->size(); ++i ) {
    G4int pdgid = result->operator[](i)->GetDefinition()->GetPDGEncoding();
    if ( pdgid == -2212 ) {  // antiproton
      antiproton.push_back( std::make_pair( i, result->operator[](i)->GetMomentum() ) );
      result->erase( result->begin() + i );
    }
  }
  for ( unsigned int i = 0; i < result->size(); ++i ) {
    G4int pdgid = result->operator[](i)->GetDefinition()->GetPDGEncoding();
    if ( pdgid == -2112 ) {  // antineutron
      antineutron.push_back( std::make_pair( i, result->operator[](i)->GetMomentum() ) );
      result->erase( result->begin() + i );
    }
  }

  for ( unsigned int i = 0; i < proton.size(); ++i ) {  // loop over protons 
    if ( proton.at(i).first == -1 ) continue;
    G4ThreeVector p1 = proton.at(i).second;
    int partner1 = FindPartner( p1, G4Proton::Proton()->GetPDGMass(), neutron,
				G4Neutron::Neutron()->GetPDGMass(), 1 );
    if ( partner1 == -1 ) {  // if no partner found, then the proton is a final-state secondary
    	G4ParticleDefinition* prt = G4ParticleTable::GetParticleTable()->FindParticle( "proton" );
    	G4ReactionProduct* finalp = new G4ReactionProduct;
    	finalp->SetDefinition( prt );
    	G4double massp = prt->GetPDGMass();
    	G4double totalEnergy = std::sqrt( p1.mag()*p1.mag() + massp*massp );	
    	finalp->SetMomentum( p1 );
    	finalp->SetTotalEnergy( totalEnergy );
    	finalp->SetMass( massp );	
    	result->push_back( finalp );		
    	continue;
    }
    G4ThreeVector p2 = neutron.at(partner1).second;
    PushDeuteron( p1, p2, 1, result );
    neutron.at(partner1).first = -1;  // tag the bound neutron
  }

  for ( unsigned int i = 0; i < neutron.size(); ++i ) {  // loop over neutrons 
    if ( neutron.at(i).first == -1 ) continue;  // Skip already bound neutron, else it is a final-state secondary
    G4ParticleDefinition* nrt = G4ParticleTable::GetParticleTable()->FindParticle( "neutron" );
    G4ReactionProduct* finaln = new G4ReactionProduct;
    finaln->SetDefinition( nrt );
    G4ThreeVector p2 = neutron.at(i).second;
    G4double massn = nrt->GetPDGMass();
    G4double totalEnergy = std::sqrt( p2.mag()*p2.mag() + massn*massn );	
    finaln->SetMomentum( p2 );
    finaln->SetTotalEnergy( totalEnergy );
    finaln->SetMass( massn );	
    result->push_back( finaln );			
  }

  for ( unsigned int i = 0; i < antiproton.size(); ++i ) {  // loop over antiprotons
    if ( antiproton.at(i).first == -1 ) continue;
    G4ThreeVector p1 = antiproton.at(i).second;
    int partner1 = FindPartner( p1, G4Proton::Proton()->GetPDGMass(), antineutron,
				G4Neutron::Neutron()->GetPDGMass(), -1 );
    if ( partner1 == -1 ) {  // if no partner found, then the antiproton is a final-state secondary
    	G4ParticleDefinition* pbar = G4ParticleTable::GetParticleTable()->FindAntiParticle( "proton" );
    	G4ReactionProduct* finalpbar = new G4ReactionProduct;
    	finalpbar->SetDefinition( pbar );
    	G4double massp = pbar->GetPDGMass();
    	G4double totalEnergy = std::sqrt( p1.mag()*p1.mag() + massp*massp );
    	finalpbar->SetMomentum( p1 );
    	finalpbar->SetTotalEnergy( totalEnergy );
    	finalpbar->SetMass( massp );	
    	result->push_back( finalpbar );		
    	continue;
    }
    G4ThreeVector p2 = antineutron.at(partner1).second;
    PushDeuteron( p1, p2, -1, result );
    antineutron.at(partner1).first = -1;  // tag the bound antineutron
  }

  for ( unsigned int i = 0; i < antineutron.size(); ++i ) {  // loop over antineutrons 
    if ( antineutron.at(i).first == -1 ) continue;  // Skip already bound antineutron, else it is a final-state secondary
    G4ParticleDefinition* nbar = G4ParticleTable::GetParticleTable()->FindAntiParticle( "neutron" );
    G4ReactionProduct* finalnbar = new G4ReactionProduct;
    finalnbar->SetDefinition( nbar );
    G4ThreeVector p2 = antineutron.at(i).second;
    G4double massn = nbar->GetPDGMass();
    G4double totalEnergy = std::sqrt( p2.mag()*p2.mag() + massn*massn );	
    finalnbar->SetMomentum( p2 );
    finalnbar->SetTotalEnergy( totalEnergy );
    finalnbar->SetMass( massn );	
    result->push_back( finalnbar );			
  }
}


void G4CRCoalescence::PushDeuteron( const G4ThreeVector &p1, const G4ThreeVector &p2, G4int charge,  // input
				    G4ReactionProductVector* result ) {                              // output
  // Create a deuteron or antideuteron (depending on "charge") object (of type G4ReactionProduct)
  // from the two input momenta "p1" and "p2", and push it to the vector "result".
  if ( charge > 0 ) {
    G4ParticleDefinition* deuteron = G4ParticleTable::GetParticleTable()->FindParticle( "deuteron" );		
    G4ReactionProduct* finaldeut = new G4ReactionProduct;
    finaldeut->SetDefinition( deuteron );
    G4ThreeVector psum = p1 + p2;
    G4double massd = deuteron->GetPDGMass();
    G4double totalEnergy = std::sqrt( psum.mag()*psum.mag() + massd*massd );	
    finaldeut->SetMomentum( psum );
    finaldeut->SetTotalEnergy( totalEnergy );
    finaldeut->SetMass( massd );
    finaldeut->SetCreatorModelID( secID );
    result->push_back( finaldeut );
  } else {
    G4ParticleDefinition* antideuteron = G4ParticleTable::GetParticleTable()->FindAntiParticle( "deuteron" );	
    G4ReactionProduct* finalantideut = new G4ReactionProduct;
    finalantideut->SetDefinition( antideuteron );
    G4ThreeVector psum = p1 + p2;
    G4double massd = antideuteron->GetPDGMass();
    G4double totalEnergy = std::sqrt( psum.mag()*psum.mag() + massd*massd );	
    finalantideut->SetMomentum( psum );
    finalantideut->SetTotalEnergy( totalEnergy );
    finalantideut->SetMass( massd );
    finalantideut->SetCreatorModelID( secID );
    result->push_back( finalantideut );
  }
}


G4int G4CRCoalescence::FindPartner( const G4ThreeVector &p1, G4double m1,
				    std::vector< std::pair< G4int, G4ThreeVector > > &neutron,
				    G4double m2, G4int charge ) {
  // Find a nucleon/antinucleon (depending on "charge") partner, from the input list "neutron"
  // (which is a vector of either neutron or antineutron particles depending on "charge")
  // within a sphere of radius p0 centered at the input momentum "p1"; exclude already bound
  // particles (neutrons or antineutrons depending on "charge") of "neutron". 
  for ( unsigned int j = 0; j < neutron.size(); ++j ) {
    if ( neutron.at(j).first == -1 ) continue;  // skip already bound particle
    G4ThreeVector p2 = neutron.at(j).second;
    if ( ! Coalescence( p1, m1, p2, m2, charge ) ) continue;
    return j;
  }
  return -1;  // no partner found
}


G4bool G4CRCoalescence::Coalescence( const G4ThreeVector &p1, G4double m1,
				     const G4ThreeVector &p2, G4double m2, G4int charge ) {
  // Returns true if the momenta of the two nucleons/antinucleons (depending on "charge") are
  // inside of an sphere of radius p0 (assuming that the two particles are in the same spatial place).
  return Coalescence( p1.x(), p1.y(), p1.z(), m1, p2.x(), p2.y(), p2.z(), m2, charge );
}


G4bool G4CRCoalescence::Coalescence( G4double p1x, G4double p1y, G4double p1z, G4double m1,
                                     G4double p2x, G4double p2y, G4double p2z, G4double m2,
				     G4int charge ) {
  // Returns true if the momenta of the two nucleons/antinucleons (depending on "charge") are
  // inside of a sphere of radius p0 (assuming that the two particles are in the same spatial place).
  G4double deltaP = GetPcm( p1x, p1y, p1z, m1, p2x, p2y, p2z, m2 );
  if ( charge > 0 ) return ( deltaP < fP0_d );
  else              return ( deltaP < fP0_dbar );
}


G4double G4CRCoalescence::GetPcm( const G4ThreeVector& p1, G4double m1,
				  const G4ThreeVector& p2, G4double m2 ) {
  // Momentum in the center-of-mass frame of two particles from LAB values.
  return GetPcm( p1.x(), p1.y(), p1.z(), m1, p2.x(), p2.y(), p2.z(), m2 );
}


G4double G4CRCoalescence::GetPcm( G4double p1x, G4double p1y, G4double p1z, G4double m1,
				  G4double p2x, G4double p2y, G4double p2z, G4double m2 ) {
  // Momentum in the center-of-mass frame of two particles from LAB values.
  G4double scm = GetS( p1x, p1y, p1z, m1, p2x, p2y, p2z, m2 );
  return std::sqrt( (scm - (m1-m2)*(m1-m2))*(scm - (m1+m2)*(m1+m2)) ) / (2.0*std::sqrt( scm ));
}


G4double G4CRCoalescence::GetS( G4double p1x, G4double p1y, G4double p1z, G4double m1,
				G4double p2x, G4double p2y, G4double p2z, G4double m2 ) {
  // Square of center-of-mass energy of two particles from LAB values.
  G4double E1 = std::sqrt( p1x*p1x + p1y*p1y + p1z*p1z + m1*m1 );
  G4double E2 = std::sqrt( p2x*p2x + p2y*p2y + p2z*p2z + m2*m2 );	
  return (E1+E2)*(E1+E2) - (p1x+p2x)*(p1x+p2x) - (p1y+p2y)*(p1y+p2y) - (p1z+p2z)*(p1z+p2z);
}
