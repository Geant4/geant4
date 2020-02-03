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
//

// ------------------------------------------------------------
//      GEANT 4 class implemetation file
//
//      ---------------- G4ElasticHNScattering --------------
//                   by V. Uzhinsky, March 2008.
//             elastic scattering used by Fritiof model
//                 Take a projectile and a target
//                 scatter the projectile and target
// ---------------------------------------------------------------------

#include "globals.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4ElasticHNScattering.hh"
#include "G4LorentzRotation.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4VSplitableHadron.hh"
#include "G4ExcitedString.hh"
#include "G4FTFParameters.hh"

#include "G4SampleResonance.hh"

#include "G4Exp.hh"
#include "G4Log.hh"

//============================================================================

G4ElasticHNScattering::G4ElasticHNScattering() {}


//============================================================================

G4bool G4ElasticHNScattering::ElasticScattering( G4VSplitableHadron* projectile, 
                                                 G4VSplitableHadron* target,
                                                 G4FTFParameters* theParameters ) const {
  projectile->IncrementCollisionCount( 1 );
  target->IncrementCollisionCount( 1 );

  if ( projectile->Get4Momentum().z() < 0.0 ) return false;  //Uzhi Aug.2019

  // Projectile parameters
  G4LorentzVector Pprojectile = projectile->Get4Momentum();
  G4double M0projectile = Pprojectile.mag();
  G4double M0projectile2 = M0projectile * M0projectile;

  // Target parameters
  G4LorentzVector Ptarget = target->Get4Momentum();
  G4double M0target = Ptarget.mag();
  G4double M0target2 = M0target * M0target;

  G4double AveragePt2 = theParameters->GetAvaragePt2ofElasticScattering();

  // Transform momenta to cms and then rotate parallel to z axis;
  G4LorentzVector Psum;
  Psum = Pprojectile + Ptarget;
  G4LorentzRotation toCms( -1*Psum.boostVector() );
  G4LorentzVector Ptmp = toCms*Pprojectile;
  if ( Ptmp.pz() <= 0.0 ) return false;                                 
  //"String" moving backwards in  CMS, abort collision !
  //G4cout << " abort Collision! " << G4endl;
  toCms.rotateZ( -1*Ptmp.phi() );
  toCms.rotateY( -1*Ptmp.theta() );
  G4LorentzRotation toLab( toCms.inverse() );
  Pprojectile.transform( toCms );
  Ptarget.transform( toCms );

  G4double PZcms2, PZcms;                                          
  G4double S = Psum.mag2();                                          
  G4double SqrtS = std::sqrt( S );
  if ( SqrtS < M0projectile + M0target ) return false;

  PZcms2 = ( S*S + sqr( M0projectile2 ) + sqr( M0target2 )
             - 2*S*M0projectile2 - 2*S*M0target2 - 2*M0projectile2*M0target2 ) / 4.0 / S;

  PZcms = ( PZcms2 > 0.0 ? std::sqrt( PZcms2 ) : 0.0 );

  G4double maxPtSquare = PZcms2;

  // Now we can calculate the transferred Pt
  G4double Pt2;                                                    
  G4double ProjMassT2, ProjMassT;                                  
  G4double TargMassT2, TargMassT;
  G4LorentzVector Qmomentum;

  const G4int maxNumberOfLoops = 1000;
  G4int loopCounter = 0;
  do {
    Qmomentum = G4LorentzVector( GaussianPt( AveragePt2, maxPtSquare ), 0.0 );
    Pt2 = G4ThreeVector( Qmomentum.vect() ).mag2();                  
    ProjMassT2 = M0projectile2 + Pt2;                           
    ProjMassT = std::sqrt( ProjMassT2 );                            
    TargMassT2 = M0target2 + Pt2;                               
    TargMassT = std::sqrt( TargMassT2 );                            
  } while ( ( SqrtS < ProjMassT + TargMassT ) && 
            ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */
  if ( loopCounter >= maxNumberOfLoops ) {
    return false;
  }

  PZcms2 = ( S*S + sqr( ProjMassT2 ) + sqr( TargMassT2 )                           
             - 2.0*S*ProjMassT2 - 2.0*S*TargMassT2 - 2.0*ProjMassT2*TargMassT2 ) / 4.0 / S;

  if ( PZcms2 < 0.0 ) { PZcms2 = 0.0; };  // to avoid the exactness problem
  PZcms = std::sqrt( PZcms2 );                                    
  Pprojectile.setPz( PZcms );  
  Ptarget.setPz( -PZcms ); 
  Pprojectile += Qmomentum;
  Ptarget     -= Qmomentum;

  // Transform back and update SplitableHadron Participant.
  Pprojectile.transform( toLab );
  Ptarget.transform( toLab );

  // Calculation of the creation time
  projectile->SetTimeOfCreation( target->GetTimeOfCreation() );
  projectile->SetPosition( target->GetPosition() );

  // Creation time and position of target nucleon were determined at
  // ReggeonCascade() of G4FTFModel

  projectile->Set4Momentum( Pprojectile );
  target->Set4Momentum( Ptarget );

  //projectile->IncrementCollisionCount( 1 );
  //target->IncrementCollisionCount( 1 );

  return true;
}


//============================================================================

G4ThreeVector G4ElasticHNScattering::GaussianPt( G4double AveragePt2, 
                                                 G4double maxPtSquare ) const {
  // @@ this method is used in FTFModel as well. Should go somewhere common!
  G4double Pt2( 0.0 );
  if ( AveragePt2 <= 0.0 ) {
    Pt2 = 0.0;
  } else {
    Pt2 = -AveragePt2 * G4Log( 1.0 + G4UniformRand() * ( G4Exp( -maxPtSquare/AveragePt2 ) -1.0 ) ); 
  }
  G4double Pt = ( Pt2 > 0.0 ? std::sqrt( Pt2 ) : 0.0 );
  G4double phi = G4UniformRand() * twopi;
  return G4ThreeVector( Pt * std::cos( phi ), Pt * std::sin( phi ), 0.0 );    
}


//============================================================================

G4ElasticHNScattering::G4ElasticHNScattering( const G4ElasticHNScattering& ) {
  throw G4HadronicException( __FILE__, __LINE__, 
                             "G4ElasticHNScattering copy constructor not meant to be called" );
}


//============================================================================

G4ElasticHNScattering::~G4ElasticHNScattering() {}


//============================================================================

const G4ElasticHNScattering & G4ElasticHNScattering::operator=( const G4ElasticHNScattering& ) {
  throw G4HadronicException( __FILE__, __LINE__, 
                             "G4ElasticHNScattering = operator not meant to be called" );
}


//============================================================================

G4bool G4ElasticHNScattering::operator==( const G4ElasticHNScattering& ) const {
 throw G4HadronicException( __FILE__, __LINE__, 
                            "G4ElasticHNScattering == operator not meant to be called" );
}


//============================================================================

G4bool G4ElasticHNScattering::operator!=( const G4ElasticHNScattering& ) const {
  throw G4HadronicException( __FILE__, __LINE__, 
                            "G4ElasticHNScattering != operator not meant to be called" );
}

