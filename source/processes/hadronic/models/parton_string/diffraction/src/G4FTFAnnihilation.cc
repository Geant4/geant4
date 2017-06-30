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
// $Id: G4FTFAnnihilation.cc 102029 2016-12-16 14:53:08Z gcosmo $
//

// ------------------------------------------------------------
//      GEANT 4 class implemetation file
//
//      ---------------- G4FTFAnnihilation --------------
//                by V. Uzhinsky, Spring 2011.
//                Take a projectile and a target
//        make annihilation or re-orangement of quarks and anti-quarks.
//     Ideas of Quark-Gluon-String model my A. Capella and A.B. Kaidalov
//                       are implemented.
// ---------------------------------------------------------------------

#include "globals.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4DiffractiveSplitableHadron.hh"
#include "G4DiffractiveExcitation.hh"
#include "G4FTFParameters.hh"
#include "G4ElasticHNScattering.hh"
#include "G4FTFAnnihilation.hh"

#include "G4LorentzRotation.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh" 
#include "G4VSplitableHadron.hh"
#include "G4ExcitedString.hh"
#include "G4ParticleTable.hh"
#include "G4Neutron.hh"
#include "G4ParticleDefinition.hh"

#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

//#include "G4ios.hh"
//#include "UZHI_diffraction.hh"

#include "G4ParticleTable.hh"     // Uzhi March 2016
//============================================================================

//#define debugFTFannih


//============================================================================

G4FTFAnnihilation::G4FTFAnnihilation() {}


//============================================================================

G4FTFAnnihilation::~G4FTFAnnihilation() {}


//============================================================================

G4bool G4FTFAnnihilation::Annihilate( G4VSplitableHadron* projectile, 
                                      G4VSplitableHadron* target,
                                      G4VSplitableHadron*& AdditionalString,
                                      G4FTFParameters* theParameters ) const  {

//theParameters->SetProbabilityOfAnnihilation( 0.0 );  // Uzhi March 2016 ??? for other Anti_bar annih.

  #ifdef debugFTFannih 
  G4cout << "---------------------------- Annihilation----------------" << G4endl;
  #endif

  // Projectile parameters
  G4LorentzVector Pprojectile = projectile->Get4Momentum();
  G4int ProjectilePDGcode = projectile->GetDefinition()->GetPDGEncoding();
  if ( ProjectilePDGcode > 0 ) {
    target->SetStatus( 3 );                                      // 2->3   Uzhi Oct 2014
    return false;
  } 
  //G4double M0projectile = Pprojectile.mag();  
  //G4double M0projectile2 = projectile->GetDefinition()->GetPDGMass() *
  //                         projectile->GetDefinition()->GetPDGMass(); 
  G4double M0projectile2 = Pprojectile.mag2();

  // Target parameters
  G4int TargetPDGcode = target->GetDefinition()->GetPDGEncoding();
  G4LorentzVector Ptarget = target->Get4Momentum();
  //G4double M0target = Ptarget.mag();
  //G4double M0target2 = target->GetDefinition()->GetPDGMass() *
  //                     target->GetDefinition()->GetPDGMass();
  G4double M0target2 = Ptarget.mag2();

  #ifdef debugFTFannih
  G4cout << "PDG codes " << ProjectilePDGcode << " " << TargetPDGcode << G4endl
         << "Pprojec " << Pprojectile << " " << Pprojectile.mag() << G4endl
         << "Ptarget " << Ptarget    << " " << Ptarget.mag() << G4endl
         << "M0 proj target " << std::sqrt( M0projectile2 ) 
         << " " << std::sqrt( M0target2 ) << G4endl;
  #endif

  G4double AveragePt2 = theParameters->GetAveragePt2();

  // Kinematical properties of the interactions
  G4LorentzVector Psum;  // 4-momentum in CMS
  Psum = Pprojectile + Ptarget;
  G4double S = Psum.mag2(); 

  #ifdef debugFTFannih
  G4cout << "Psum SqrtS S " << Psum << " " << std::sqrt( S ) << " " << S << G4endl;
  #endif

  // Transform momenta to cms and then rotate parallel to z axis
  G4LorentzRotation toCms( -1*Psum.boostVector() );
  G4LorentzVector Ptmp = toCms*Pprojectile;
  toCms.rotateZ( -1*Ptmp.phi() );
  toCms.rotateY( -1*Ptmp.theta() );
  G4LorentzRotation toLab( toCms.inverse() );

  G4double SqrtS = std::sqrt( S );
  G4double maxPtSquare;
  G4double X_a( 0.0 ), X_b( 0.0 ), X_c( 0.0 ), X_d( 0.0 );

  G4double MesonProdThreshold = projectile->GetDefinition()->GetPDGMass() +
                                target->GetDefinition()->GetPDGMass() +
                                ( 2.0*140.0 + 16.0 )*MeV;  // 2 Mpi + DeltaE
  G4double Prel2 = S*S + M0projectile2*M0projectile2 + M0target2*M0target2 -
                   2.0*S*M0projectile2 - 2.0*S*M0target2 - 2.0*M0projectile2*M0target2;
  Prel2 /= S;
  //G4cout << "Prel2 " << Prel2 << G4endl;
  if ( Prel2 <= 0.0 ) {  // *MeV*MeV 1600.
    // Annihilation at rest! Values are copied from Parameters
    X_a =  625.1;    // mb  // 3-shirt diagram
    X_b =    0.0;    // 9.780 12 Dec. 2012;  // mb  // anti-quark-quark annihilation
    X_c =   49.989;  // mb
    X_d =    6.614;  // mb

    #ifdef debugFTFannih 
    G4cout << "Annih at Rest X a b c d " << X_a << " " << X_b << " " << X_c << " " << X_d 
           << G4endl;
    #endif

  } else { // Annihilation in flight!
    G4double FlowF = 1.0 / std::sqrt( Prel2 )*GeV;

    // Process cross sections
    X_a = 25.0*FlowF;  // mb 3-shirt diagram
    if ( SqrtS < MesonProdThreshold ) {
      X_b = 3.13 + 140.0*G4Pow::GetInstance()->powA( ( MesonProdThreshold - SqrtS )/GeV, 2.5 ); 
    } else {
      X_b = 6.8*GeV / SqrtS;  // mb anti-quark-quark annihilation
    }
    if ( projectile->GetDefinition()->GetPDGMass() + target->GetDefinition()->GetPDGMass()
         > SqrtS ) {
      X_b = 0.0;
    }
    // This can be in an interaction of low energy anti-baryon with off-shell nuclear nucleon
    X_c = 2.0 * FlowF * sqr( projectile->GetDefinition()->GetPDGMass() +
                             target->GetDefinition()->GetPDGMass() ) / S; // mb re-arrangement of
                                                                          // 2 quarks and 2 anti-quarks
    X_d = 23.3*GeV*GeV / S; // mb anti-quark-quark string creation

    #ifdef debugFTFannih 
    G4cout << "Annih in Flight X a b c d " << X_a << " " << X_b << " " << X_c << " " << X_d
           << G4endl << "SqrtS MesonProdThreshold " << SqrtS << " " << MesonProdThreshold
           << G4endl;
    #endif

  } 

  if        ((ProjectilePDGcode == -2212 || ProjectilePDGcode == -2214)&& ( TargetPDGcode == 2212 || TargetPDGcode == 2214 ) ) { 
    X_b *= 5.0; X_c *= 5.0; X_d *= 6.0;  // Pbar P
  } else if ((ProjectilePDGcode == -2212 || ProjectilePDGcode == -2214)&& ( TargetPDGcode == 2112 || TargetPDGcode == 2114 ) ) {
    X_b *= 4.0; X_c *= 4.0; X_d *= 4.0;  // Pbar N
  } else if ((ProjectilePDGcode == -2112 || ProjectilePDGcode == -2114)&& ( TargetPDGcode == 2212 || TargetPDGcode == 2214 ) ) {
    X_b *= 4.0; X_c *= 4.0; X_d *= 4.0;  // NeutrBar P
  } else if ((ProjectilePDGcode == -2112 || ProjectilePDGcode == -2114)&& ( TargetPDGcode == 2112 || TargetPDGcode == 2114 ) ) {
    X_b *= 5.0; X_c *= 5.0; X_d *= 6.0;  // NeutrBar N
  } else if ((ProjectilePDGcode == -3122 || ProjectilePDGcode == -3124)&& ( TargetPDGcode == 2212 || TargetPDGcode == 2214 ) ) {
    X_b *= 3.0; X_c *= 3.0; X_d *= 2.0;  // LambdaBar P
  } else if ((ProjectilePDGcode == -3122 || ProjectilePDGcode == -3124)&& ( TargetPDGcode == 2112 || TargetPDGcode == 2114 ) ) {
    X_b *= 3.0; X_c *= 3.0; X_d *= 2.0;  // LambdaBar N
  } else if ((ProjectilePDGcode == -3112 || ProjectilePDGcode == -3114)&& ( TargetPDGcode == 2212 || TargetPDGcode == 2214 ) ) {
    X_b *= 2.0; X_c *= 2.0; X_d *= 0.0;  // Sigma-Bar P
  } else if ((ProjectilePDGcode == -3112 || ProjectilePDGcode == -3114)&& ( TargetPDGcode == 2112 || TargetPDGcode == 2114 ) ) {
    X_b *= 4.0; X_c *= 4.0; X_d *= 2.0;  // Sigma-Bar N
  } else if ((ProjectilePDGcode == -3212 || ProjectilePDGcode == -3214)&& ( TargetPDGcode == 2212 || TargetPDGcode == 2214 ) ) {
    X_b *= 3.0; X_c *= 3.0; X_d *= 2.0;  // Sigma0Bar P
  } else if ((ProjectilePDGcode == -3212 || ProjectilePDGcode == -3214)&& ( TargetPDGcode == 2112 || TargetPDGcode == 2114 ) ) {
    X_b *= 3.0; X_c *= 3.0; X_d *= 2.0;  // Sigma0Bar N
  } else if ((ProjectilePDGcode == -3222 || ProjectilePDGcode == -3224)&& ( TargetPDGcode == 2212 || TargetPDGcode == 2214 ) ) {
    X_b *= 4.0; X_c *= 4.0; X_d *= 2.0;  // Sigma+Bar P
  } else if ((ProjectilePDGcode == -3222 || ProjectilePDGcode == -3224)&& ( TargetPDGcode == 2112 || TargetPDGcode == 2114 ) ) {
    X_b *= 2.0; X_c *= 2.0; X_d *= 0.0;  // Sigma+Bar N
  } else if ((ProjectilePDGcode == -3312 || ProjectilePDGcode == -3314)&& ( TargetPDGcode == 2212 || TargetPDGcode == 2214 ) ) {
    X_b *= 1.0; X_c *= 1.0; X_d *= 0.0;  // Xi-Bar P
  } else if ((ProjectilePDGcode == -3312 || ProjectilePDGcode == -3314)&& ( TargetPDGcode == 2112 || TargetPDGcode == 2114 ) ) {
    X_b *= 2.0; X_c *= 2.0; X_d *= 0.0;  // Xi-Bar N
  } else if ((ProjectilePDGcode == -3322 || ProjectilePDGcode == -3324)&& ( TargetPDGcode == 2212 || TargetPDGcode == 2214 ) ) {
    X_b *= 2.0; X_c *= 2.0; X_d *= 0.0;  // Xi0Bar P
  } else if ((ProjectilePDGcode == -3322 || ProjectilePDGcode == -3324)&& ( TargetPDGcode == 2112 || TargetPDGcode == 2114 ) ) {
    X_b *= 1.0; X_c *= 1.0; X_d *= 0.0;  // Xi0Bar N
  } else if ( ProjectilePDGcode == -3334 && ( TargetPDGcode == 2212 || TargetPDGcode == 2214 ) ) {
    X_b *= 0.0; X_c *= 0.0; X_d *= 0.0;  // Omega-Bar P
  } else if ( ProjectilePDGcode == -3334 && ( TargetPDGcode == 2112 || TargetPDGcode == 2114 ) ) {
    X_b *= 0.0; X_c *= 0.0; X_d *= 0.0;  // Omega-Bar N
  } else {
    G4cout << "Unknown anti-baryon for FTF annihilation: PDGcodes - "
           << ProjectilePDGcode << " " << TargetPDGcode << G4endl;
  }

  #ifdef debugFTFannih 
  G4cout << "Annih Actual X a b c d " << X_a << " " << X_b << " " << X_c << " " << X_d << G4endl;
  #endif

  G4double Xannihilation = X_a + X_b + X_c + X_d;
//X_a=0.;    // Uzhi
//X_b=0.; 
//X_c=0.;
//X_d=0.;
//Xannihilation = X_a + X_b + X_c + X_d;

  // Projectile unpacking
  G4int AQ[3];
  UnpackBaryon( ProjectilePDGcode, AQ[0], AQ[1], AQ[2] );

  // Target unpacking
  G4int Q[3];
  UnpackBaryon( TargetPDGcode, Q[0], Q[1], Q[2] ); 

  G4double Ksi = G4UniformRand();

  if ( Ksi < X_a / Xannihilation ) {

    // Simulation of 3 anti-quark-quark strings creation
    //    Sampling of anti-quark order in projectile

    #ifdef debugFTFannih 
    G4cout << "Process a, 3 shirt diagram" << G4endl;
    #endif

    G4int SampledCase = G4RandFlat::shootInt( G4long( 6 ) );

    G4int Tmp1( 0 ), Tmp2( 0 );
    if ( SampledCase == 0 ) {                                    
    } else if ( SampledCase == 1 ) { 
      Tmp1 = AQ[1]; AQ[1] = AQ[2]; AQ[2] = Tmp1;
    } else if ( SampledCase == 2 ) { 
      Tmp1 = AQ[0]; AQ[0] = AQ[1]; AQ[1] = Tmp1; 
    } else if ( SampledCase == 3 ) { 
      Tmp1 = AQ[0]; Tmp2 = AQ[1];  AQ[0] = AQ[2]; AQ[1] = Tmp1;  AQ[2] = Tmp2;
    } else if ( SampledCase == 4 ) { 
      Tmp1 = AQ[0]; Tmp2 = AQ[1];  AQ[0] = Tmp2;  AQ[1] = AQ[2]; AQ[2] = Tmp1; 
    } else if ( SampledCase == 5 ) { 
      Tmp1 = AQ[0]; Tmp2 = AQ[1];  AQ[0] = AQ[2]; AQ[1] = Tmp2;  AQ[2] = Tmp1;
    }  

    //  Set the string properties

    //G4cout << "String 1 " << AQ[0] << " " << Q[0] << G4endl;
    projectile->SplitUp();
    projectile->SetFirstParton( AQ[0] );
    projectile->SetSecondParton( Q[0] );
    projectile->SetStatus( 0 );

// Uzhi March 2016 start
G4int aAQ, aQ;
aAQ=std::abs( AQ[0] ); aQ=std::abs( Q[0] );
G4int NewCode;
G4double aKsi = G4UniformRand();

if ( aAQ == aQ )
{
 if ( aAQ != 3 )
 {
  NewCode = 111;                       // Pi0-meson
  if ( aKsi < 0.5 )
  {
   NewCode = 221;                     // Eta -meson
   if ( aKsi < 0.25 ) {NewCode = 331;} // Eta'-meson
  }
 } else
 {
  NewCode = 221;                      // Eta -meson
  if( aKsi < 0.5 ) {NewCode = 331;}    // Eta'-meson
 }
} else
{
 if ( aAQ > aQ ){ NewCode = aAQ*100 + aQ*10 + 1; NewCode *= aAQ/AQ[0]; } 
 else           { NewCode = aQ*100 + aAQ*10 + 1; NewCode *=  aQ/Q[0];  }
}

G4ParticleDefinition* TestParticle = G4ParticleTable::GetParticleTable()->FindParticle( NewCode );
if(!TestParticle) return false;
projectile->SetDefinition( TestParticle );

theParameters->SetProjMinDiffMass( 0.5 );     // Uzhi 2016 M+140 ??
theParameters->SetProjMinNonDiffMass( 0.5 );  // Uzhi 2016 M+140 ??
// Uzhi March 2016 end

    //G4cout << "String 2 " << Q[1] << " " << AQ[1] << G4endl;
    target->SplitUp();
    target->SetFirstParton( Q[1] );
    target->SetSecondParton( AQ[1] );
    target->SetStatus( 0 );

// Uzhi March 2016 Start
aAQ=std::abs( AQ[1] ); aQ=std::abs( Q[1] ); aKsi = G4UniformRand();
if ( aAQ == aQ )
{
 if ( aAQ != 3 )
 {
  NewCode = 111;                       // Pi0-meson
  if ( aKsi < 0.5 )
  {
   NewCode = 221;                     // Eta -meson
   if ( aKsi < 0.25 ) {NewCode = 331;} // Eta'-meson
  }
 } else
 {
  NewCode = 221;                      // Eta -meson
  if( aKsi < 0.5 ) {NewCode = 331;}    // Eta'-meson
 }
} else
{
 if ( aAQ > aQ ){ NewCode = aAQ*100 + aQ*10 + 1; NewCode *= aAQ/AQ[1]; } 
 else           { NewCode = aQ*100 + aAQ*10 + 1; NewCode *=  aQ/Q[1];  }
}
TestParticle = G4ParticleTable::GetParticleTable()->FindParticle( NewCode );
if(!TestParticle) return false;
target->SetDefinition( TestParticle );

theParameters->SetTarMinDiffMass( 0.5 );    // Uzhi 2016 M+140 ??
theParameters->SetTarMinNonDiffMass( 0.5 ); // Uzhi 2016 M+140 ??
// Uzhi March 2016 end


    //G4cout << "String 3 " << AQ[2] << " " << Q[2] << G4endl;
    AdditionalString = new G4DiffractiveSplitableHadron();

// Uzhi March 2016 start
aAQ=std::abs( AQ[2] ); aQ=std::abs( Q[2] ); aKsi = G4UniformRand();

if ( aAQ == aQ )
{
 if ( aAQ != 3 )
 {
  NewCode = 111;                       // Pi0-meson
  if ( aKsi < 0.5 )
  {
   NewCode = 221;                     // Eta -meson
   if ( aKsi < 0.25 ) {NewCode = 331;} // Eta'-meson
  }
 } else
 {
  NewCode = 221;                      // Eta -meson
  if( aKsi < 0.5 ) {NewCode = 331;}    // Eta'-meson
 }
} else
{
 if ( aAQ > aQ ){ NewCode = aAQ*100 + aQ*10 + 1; NewCode *= aAQ/AQ[2]; } 
 else           { NewCode = aQ*100 + aAQ*10 + 1; NewCode *=  aQ/Q[2];  }
}

TestParticle = G4ParticleTable::GetParticleTable()->FindParticle( NewCode );
if(!TestParticle) return false;
AdditionalString->SetDefinition( TestParticle );
// Uzhi March 2016 end

    AdditionalString->SplitUp();
    AdditionalString->SetFirstParton( AQ[2] );
    AdditionalString->SetSecondParton( Q[2] );
    AdditionalString->SetStatus( 0 );
    //G4cout << G4endl << "*AdditionalString in Annih" << AdditionalString << G4endl;

    // Sampling kinematical properties
    // 1 string AQ[0]-Q[0]// 2 string AQ[1]-Q[1]// 3 string AQ[2]-Q[2]

    G4ThreeVector Quark_Mom[6];
    G4double ModMom2[6];  //ModMom[6]

    AveragePt2 = 200.0*200.0; maxPtSquare = S;

    G4double SumMt( 0.0 );
    G4double MassQ2 = 0.0;  // 100.0*100.0*MeV*MeV;
    G4int NumberOfTries( 0 );
    G4double ScaleFactor( 1.0 );

    const G4int maxNumberOfLoops = 1000;
    G4int loopCounter = 0;
    do {
      NumberOfTries++;
      if ( NumberOfTries == 100*(NumberOfTries/100) ) {
        // At large number of tries it would be better to reduce the values of <Pt^2>
        ScaleFactor /= 2.0;
        AveragePt2 *= ScaleFactor;
      }
      G4ThreeVector PtSum( 0.0, 0.0, 0.0 );
      for ( G4int i = 0; i < 6; i++ ) {
        Quark_Mom [i] = GaussianPt( AveragePt2, maxPtSquare );
        PtSum += Quark_Mom[i];
      }
      PtSum /= 6.0;
      SumMt = 0.0;    
      for( G4int i = 0; i < 6; i++ ) {
        Quark_Mom[i] -= PtSum;
        //ModMom[i] = Quark_Mom[i].mag();
        ModMom2[i] = Quark_Mom[i].mag2();
        SumMt += std::sqrt( ModMom2[i] + MassQ2 );
      }
    } while ( ( SumMt > SqrtS ) && 
              ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */
    if ( loopCounter >= maxNumberOfLoops ) {
      return false;
    }

    G4double WminusTarget( 0.0 ), WplusProjectile( 0.0 );

    // Closed is variant with sampling of Xs at minimum
    //G4double SumMod_anti = ModMom[0] + ModMom[1] + ModMom[2];
    //Quark_Mom[0].setZ( ModMom[0]/SumMod_anti );
    //Quark_Mom[1].setZ( ModMom[1]/SumMod_anti );
    //Quark_Mom[2].setZ( ModMom[2]/SumMod_anti );
    //G4double SumMod_bary = ModMom[3] + ModMom[4] + ModMom[5];
    //Quark_Mom[3].setZ( ModMom[3]/SumMod_bary );
    //Quark_Mom[4].setZ( ModMom[4]/SumMod_bary );
    //Quark_Mom[5].setZ( ModMom[5]/SumMod_bary );
    //G4double Alfa = SumMod_anti*SumMod_anti;
    //G4double Beta = SumMod_bary*SumMod_bary;
    //G4double DecayMomentum2 = S*S + Alfa*Alfa + Beta*Beta
    //                          - 2.0*S*Alfa - 2.0*S*Beta - 2.0*Alfa*Beta;  
    //WminusTarget = ( S - Alfa + Beta + std::sqrt( DecayMomentum2 ) )/2.0/SqrtS; 
    //WplusProjectile = SqrtS - Beta/WminusTarget;
    // Closed is variant with sampling of Xs at minimum

    // Sampling X's of anti-baryon
    G4double Alfa_R = 0.5;
    NumberOfTries = 0;
    ScaleFactor = 1.0;
    G4bool Succes( true );

    loopCounter = 0;
    do {

      Succes = true;
      NumberOfTries++;
      if ( NumberOfTries == 100*(NumberOfTries/100) ) { 
        // At large number of tries it would be better to reduce the values of Pt's
        ScaleFactor /= 2.0;
      }

      if ( Alfa_R == 1.0 ) {
        G4double Xaq1 = 1.0 - std::sqrt( G4UniformRand() );
        G4double Xaq2 = (1.0 - Xaq1) * G4UniformRand();
        G4double Xaq3 = 1.0 - Xaq1 - Xaq2;
        Quark_Mom[0].setZ( Xaq1 ); Quark_Mom[1].setZ( Xaq2 ); Quark_Mom[2].setZ( Xaq3 );
      } else {
        G4double Xaq1 = sqr( G4UniformRand() );
        G4double Xaq2 = (1.0 - Xaq1)*sqr( std::sin( pi/2.0*G4UniformRand() ) );
        G4double Xaq3 = 1.0 - Xaq1 - Xaq2;
        Quark_Mom[0].setZ( Xaq1 ); Quark_Mom[1].setZ( Xaq2 ); Quark_Mom[2].setZ( Xaq3 );
      }

      // Sampling X's of baryon
      if ( Alfa_R == 1.0 ) {
        G4double Xq1 = 1.0 - std::sqrt( G4UniformRand() );
        G4double Xq2 = (1.0 - Xq1) * G4UniformRand();
        G4double Xq3 = 1.0 - Xq1 - Xq2;
        Quark_Mom[3].setZ( Xq1 ); Quark_Mom[4].setZ( Xq2 ); Quark_Mom[5].setZ( Xq3 );
      } else {
        G4double Xq1 = sqr( G4UniformRand() );
        G4double Xq2 = (1.0 - Xq1) * sqr( std::sin( pi/2.0*G4UniformRand() ) );
        G4double Xq3 = 1.0 - Xq1 - Xq2;
        Quark_Mom[3].setZ( Xq1 ); Quark_Mom[4].setZ( Xq2 ); Quark_Mom[5].setZ( Xq3 );
      }

      G4double Alfa( 0.0 ), Beta( 0.0 );
      for ( G4int i = 0; i < 3; i++ ) {  // For Anti-baryon
        if ( Quark_Mom[i].getZ() != 0.0 ) { 
          Alfa += ( ScaleFactor * ModMom2[i] + MassQ2 ) / Quark_Mom[i].getZ();
        } else {
          Succes = false;
        }
      } 
      for ( G4int i = 3; i < 6; i++ ) {  // For baryon
        if ( Quark_Mom[i].getZ() != 0.0 ) {
          Beta += ( ScaleFactor * ModMom2[i] + MassQ2 ) / Quark_Mom[i].getZ();
        } else {
          Succes = false;
        }
      } 

      if ( ! Succes ) continue;

      if ( std::sqrt( Alfa ) + std::sqrt( Beta ) > SqrtS ) {
        Succes = false; 
        continue;
      }

      G4double DecayMomentum2 = S*S + Alfa*Alfa + Beta*Beta
                              - 2.0*S*Alfa - 2.0*S*Beta - 2.0*Alfa*Beta;
      
      WminusTarget = ( S - Alfa + Beta + std::sqrt( DecayMomentum2 ) ) / 2.0 / SqrtS; 
      WplusProjectile = SqrtS - Beta/WminusTarget;

    } while ( ( ! Succes ) &&
              ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */
    if ( loopCounter >= maxNumberOfLoops ) {
      return false;
    }

    G4double SqrtScaleF = std::sqrt( ScaleFactor );
    for ( G4int i = 0; i < 3; i++ ) {
      G4double Pz = WplusProjectile * Quark_Mom[i].getZ() / 2.0 -
                    ( ScaleFactor * ModMom2[i] + MassQ2 ) / 
                    ( 2.0 * WplusProjectile * Quark_Mom[i].getZ() ); 
      Quark_Mom[i].setZ( Pz );
      if ( ScaleFactor != 1.0 ) {
        Quark_Mom[i].setX( SqrtScaleF * Quark_Mom[i].getX() ); 
        Quark_Mom[i].setY( SqrtScaleF * Quark_Mom[i].getY() );
      }
    }
    for ( G4int i = 3; i < 6; i++ ) {
      G4double Pz = -WminusTarget * Quark_Mom[i].getZ() / 2.0 +
                     ( ScaleFactor * ModMom2[i] + MassQ2 ) / 
                     ( 2.0 * WminusTarget * Quark_Mom[i].getZ() );
      Quark_Mom[i].setZ( Pz );
      if ( ScaleFactor != 1.0 ) {
        Quark_Mom[i].setX( SqrtScaleF * Quark_Mom[i].getX() ); 
        Quark_Mom[i].setY( SqrtScaleF * Quark_Mom[i].getY() );
      }
    }
    //G4cout << "Sum AQ " << Quark_Mom[0] + Quark_Mom[1] + Quark_Mom[2] << G4endl
    //       << "Sum Q  " << Quark_Mom[3] + Quark_Mom[4] + Quark_Mom[5] << G4endl;

    G4ThreeVector tmp = Quark_Mom[0] + Quark_Mom[3];
    G4LorentzVector Pstring1( tmp, std::sqrt( Quark_Mom[0].mag2() + MassQ2 ) +
                                   std::sqrt( Quark_Mom[3].mag2() + MassQ2 ) );
    G4double Ystring1 = Pstring1.rapidity();

    //G4cout << "Mom 1 string " << G4endl << Quark_Mom[0] << G4endl << Quark_Mom[3] << G4endl
    //       << tmp << " " << tmp.mag() << G4endl;
    //G4cout << "1 str " << Pstring1 << " " << Pstring1.mag() << " " << Ystring1 << G4endl;

    tmp = Quark_Mom[1] + Quark_Mom[4];
    G4LorentzVector Pstring2( tmp, std::sqrt( Quark_Mom[1].mag2() + MassQ2 ) +
                                   std::sqrt( Quark_Mom[4].mag2() + MassQ2 ) );
    G4double Ystring2 = Pstring2.rapidity();

    //G4cout << "Mom 2 string " << G4endl << Quark_Mom[1] << G4endl << Quark_Mom[4] << G4endl
    //       << tmp << " " << tmp.mag() << G4endl;
    //G4cout << "2 str " << Pstring2 << " " << Pstring2.mag() << " " << Ystring2 << G4endl;

    tmp = Quark_Mom[2] + Quark_Mom[5];
    G4LorentzVector Pstring3( tmp, std::sqrt( Quark_Mom[2].mag2() + MassQ2 ) +
                                   std::sqrt( Quark_Mom[5].mag2() + MassQ2 ) );
    G4double Ystring3 = Pstring3.rapidity();

    //G4cout << "Mom 3 string " << G4endl << Quark_Mom[2] << G4endl << Quark_Mom[5] << G4endl
    //       << tmp << " " << tmp.mag() << G4endl;
    //G4cout << "3 str " << Pstring3 << " " << Pstring3.mag() << " " << Ystring3 << G4endl
    //       << "SumE " << Pstring1.e() + Pstring2.e() + Pstring3.e() << G4endl
    //       << Pstring1.mag() << " " <<Pstring2.mag() << " " << Pstring3.mag() << G4endl;
    //G4int Uzhi; G4cin >> Uzhi;

    G4LorentzVector LeftString( 0.0, 0.0, 0.0, 0.0 );
    if ( Ystring1 > Ystring2  &&  Ystring2 > Ystring3 ) {
      Pprojectile = Pstring1;
      LeftString  = Pstring2;
      Ptarget     = Pstring3;
    }
    if ( Ystring1 > Ystring3  &&  Ystring3 > Ystring2 ) {
      Pprojectile = Pstring1;
      LeftString  = Pstring3;
      Ptarget     = Pstring2;
    }

    if ( Ystring2 > Ystring1  &&  Ystring1 > Ystring3 ) {
      Pprojectile = Pstring2;
      LeftString  = Pstring1;
      Ptarget     = Pstring3;
    }  
    if ( Ystring2 > Ystring3  &&  Ystring3 > Ystring1 ) {
      Pprojectile = Pstring2;
      LeftString  = Pstring3;
      Ptarget     = Pstring1;
    }

    if ( Ystring3 > Ystring1  &&  Ystring1 > Ystring2 ) {
      Pprojectile = Pstring3;
      LeftString  = Pstring1;
      Ptarget     = Pstring2;
    }
    if ( Ystring3 > Ystring2  &&  Ystring2 > Ystring1 ) {
      Pprojectile = Pstring3;
      LeftString  = Pstring2;
      Ptarget     = Pstring1;
    }
    //G4cout << "SumP " << Pprojectile + LeftString + Ptarget << " " << SqrtS << G4endl;

    Pprojectile.transform( toLab );
    LeftString.transform( toLab );
    Ptarget.transform( toLab );
    //G4cout << "SumP " << Pprojectile + LeftString + Ptarget << " " << SqrtS << G4endl;

    // Calculation of the creation time
    projectile->SetTimeOfCreation( target->GetTimeOfCreation() );
    projectile->SetPosition( target->GetPosition() );
    AdditionalString->SetTimeOfCreation( target->GetTimeOfCreation() );
    AdditionalString->SetPosition( target->GetPosition() );
    // Creation time and position of target nucleon were determined in
    // ReggeonCascade() of G4FTFModel

    //G4cout << "Mproj " << Pprojectile.mag() << G4endl << "Mtarg " << Ptarget.mag() << G4endl;
    projectile->Set4Momentum( Pprojectile );
    AdditionalString->Set4Momentum( LeftString );
    target->Set4Momentum( Ptarget );
    projectile->IncrementCollisionCount( 1 );
    AdditionalString->IncrementCollisionCount( 1 );
    target->IncrementCollisionCount( 1 );

theParameters->SetProbabilityOfAnnihilation( 0.0 );  // Uzhi March 2016

    return true;

  }  // End of if ( Ksi < X_a / Xannihilation )

  // Simulation of anti-diquark-diquark string creation

  if ( Ksi < (X_a + X_b) / Xannihilation ) {

    #ifdef debugFTFannih 
    G4cout << "Process b, quark - anti-quark annihilation, di-q - anti-di-q string" << G4endl;
    #endif

    G4int CandidatsN( 0 ), CandAQ[9][2], CandQ[9][2];
    G4int LeftAQ1( 0 ), LeftAQ2( 0 ), LeftQ1( 0 ), LeftQ2( 0 );

    for ( G4int iAQ = 0; iAQ < 3; iAQ++ ) {
      for ( G4int iQ = 0; iQ < 3; iQ++ ) {
        if ( -AQ[iAQ] == Q[iQ] ) {
          if ( iAQ == 0 ) { CandAQ[CandidatsN][0] = 1; CandAQ[CandidatsN][1] = 2; }
          if ( iAQ == 1 ) { CandAQ[CandidatsN][0] = 0; CandAQ[CandidatsN][1] = 2; }
          if ( iAQ == 2 ) { CandAQ[CandidatsN][0] = 0; CandAQ[CandidatsN][1] = 1; }
          if ( iQ  == 0 ) { CandQ[CandidatsN][0]  = 1; CandQ[CandidatsN][1]  = 2; }
          if ( iQ  == 1 ) { CandQ[CandidatsN][0]  = 0; CandQ[CandidatsN][1]  = 2; }
          if ( iQ  == 2 ) { CandQ[CandidatsN][0]  = 0; CandQ[CandidatsN][1]  = 1; }
          CandidatsN++;
        }
      }
    }
    //G4cout << "CandidatsN " << CandidatsN << G4endl;

    if ( CandidatsN != 0 ) {
      G4int SampledCase = G4RandFlat::shootInt( G4long( CandidatsN ) );
      LeftAQ1 = AQ[ CandAQ[SampledCase][0] ];
      LeftAQ2 = AQ[ CandAQ[SampledCase][1] ];
      LeftQ1  =  Q[ CandQ[SampledCase][0] ];
      LeftQ2  =  Q[ CandQ[SampledCase][1] ];

      // Build anti-diquark and diquark
      G4int Anti_DQ( 0 ), DQ( 0 );
      if ( std::abs( LeftAQ1 ) > std::abs( LeftAQ2 ) ) { 
        Anti_DQ = 1000*LeftAQ1 + 100*LeftAQ2 - 3;  // 1
      } else {
        Anti_DQ = 1000*LeftAQ2 + 100*LeftAQ1 - 3;  // 1
      }
      //if ( G4UniformRand() > 0.5 ) Anti_DQ -= 2;
      if ( std::abs( LeftQ1 ) > std::abs( LeftQ2 ) ) { 
        DQ = 1000*LeftQ1 + 100*LeftQ2 + 3;  // 1
      } else {
        DQ = 1000*LeftQ2 + 100*LeftQ1 + 3;  // 1
      }
      // if ( G4UniformRand() > 0.5 ) DQ += 2;

      // Set the string properties
      //G4cout << "Left ADiQ DiQ " << Anti_DQ << " " << DQ << G4endl;
      projectile->SplitUp();
      //projectile->SetFirstParton( Anti_DQ );
      //projectile->SetSecondParton( DQ );
      projectile->SetFirstParton( DQ );
      projectile->SetSecondParton( Anti_DQ );
      projectile->SetStatus( 0 );
      target->SetStatus( 4 );  // The target nucleon has annihilated 3->4 Uzhi Oct 2014
      Pprojectile.setPx( 0.0 );
      Pprojectile.setPy( 0.0 );
      Pprojectile.setPz( 0.0 );
      Pprojectile.setE( SqrtS );
      Pprojectile.transform( toLab );
// Uzhi March 2016 if QQ_QQbar will interact Set Mmin, MdifMin

      // Calculation of the creation time
      projectile->SetTimeOfCreation( target->GetTimeOfCreation() );
      projectile->SetPosition( target->GetPosition() );
      // Creation time and position of target nucleon were determined in
      // ReggeonCascade() of G4FTFModel

      //G4cout << "Mproj " << Pprojectile.mag() << G4endl
      //       << "Mtarg " << Ptarget.mag() << G4endl;
      projectile->Set4Momentum( Pprojectile );

      projectile->IncrementCollisionCount( 1 );
      target->IncrementCollisionCount( 1 );

//theParameters->SetProbabilityOfAnnihilation( 0.0 );  // Uzhi March 2016
// In the case baryon and anti-baryon are created. Thus the antibaryon can annihilate latter.

      return true;
    }

  }  // End of if ( Ksi < (X_a + X_b) / Xannihilation )

  if ( Ksi < ( X_a + X_b + X_c ) / Xannihilation ) {

    // Simulation of 2 anti-quark-quark strings creation

    #ifdef debugFTFannih 
    G4cout << "Process c, quark - anti-quark and string junctions annihilation, 2 strings left."
           << G4endl;
    #endif

    G4int CandidatsN( 0 ), CandAQ[9][2], CandQ[9][2];
    G4int LeftAQ1( 0 ), LeftAQ2( 0 ), LeftQ1( 0 ), LeftQ2( 0 );

    for ( G4int iAQ = 0; iAQ < 3; iAQ++ ) {
      for ( G4int iQ = 0; iQ < 3; iQ++ ) {
        if ( -AQ[iAQ] == Q[iQ] ) {
          if ( iAQ == 0 ) { CandAQ[CandidatsN][0] = 1; CandAQ[CandidatsN][1] = 2; }
          if ( iAQ == 1 ) { CandAQ[CandidatsN][0] = 0; CandAQ[CandidatsN][1] = 2; }
          if ( iAQ == 2 ) { CandAQ[CandidatsN][0] = 0; CandAQ[CandidatsN][1] = 1; }
          if ( iQ  == 0 ) { CandQ[CandidatsN][0]  = 1; CandQ[CandidatsN][1]  = 2; }
          if ( iQ  == 1 ) { CandQ[CandidatsN][0]  = 0; CandQ[CandidatsN][1]  = 2; }
          if ( iQ  == 2 ) { CandQ[CandidatsN][0]  = 0; CandQ[CandidatsN][1]  = 1; }
          CandidatsN++;
        }
      }
    }
    //G4cout << "CandidatsN " << CandidatsN << G4endl;

    if ( CandidatsN != 0 ) {
      G4int SampledCase = G4RandFlat::shootInt( G4long( CandidatsN ) );
      LeftAQ1 = AQ[ CandAQ[SampledCase][0] ];
      LeftAQ2 = AQ[ CandAQ[SampledCase][1] ];
      if ( G4UniformRand() < 0.5 ) {
        LeftQ1 = Q[ CandQ[SampledCase][0] ];
        LeftQ2 = Q[ CandQ[SampledCase][1] ];
      } else {
        LeftQ2 = Q[ CandQ[SampledCase][0] ];
        LeftQ1 = Q[ CandQ[SampledCase][1] ];
      }

      // Set the string properties
      //G4cout << "String 1 " << LeftAQ1 << " " << LeftQ1 << G4endl;
      projectile->SplitUp();
      projectile->SetFirstParton( LeftAQ1 );
      projectile->SetSecondParton( LeftQ1 );
      projectile->SetStatus( 0 );

// Uzhi March 2016 start
G4int aAQ, aQ;
aAQ=std::abs( LeftAQ1 ); aQ=std::abs( LeftQ1 );

G4int NewCode;
G4double aKsi = G4UniformRand();

if ( aAQ == aQ )
{
 if ( aAQ != 3 )
 {
  NewCode = 111;                       // Pi0-meson
  if ( aKsi < 0.5 )
  {
   NewCode = 221;                     // Eta -meson
   if ( aKsi < 0.25 ) {NewCode = 331;} // Eta'-meson
  }
 } else
 {
  NewCode = 221;                      // Eta -meson
  if( aKsi < 0.5 ) {NewCode = 331;}    // Eta'-meson
 }
} else
{
 if ( aAQ > aQ ){ NewCode = aAQ*100 + aQ*10 + 1; NewCode *= aAQ/LeftAQ1; } 
 else           { NewCode = aQ*100 + aAQ*10 + 1; NewCode *=  aQ/LeftQ1;  }
}

G4ParticleDefinition* TestParticle = G4ParticleTable::GetParticleTable()->FindParticle( NewCode );
if(!TestParticle) return false;
projectile->SetDefinition( TestParticle );
theParameters->SetProjMinDiffMass( 0.5 );            // (0.5)  // GeV    Uzhi March 2016 ???
theParameters->SetProjMinNonDiffMass( 0.5 );
// Uzhi March 2016 end

      //G4cout << "String 2 " << LeftAQ2 << " " << LeftQ2 << G4endl;
      target->SplitUp();
      target->SetFirstParton( LeftQ2 );
      target->SetSecondParton( LeftAQ2 );
      target->SetStatus( 0 );

// Uzhi March 2016 start
aAQ=std::abs( LeftAQ2 ); aQ=std::abs( LeftQ2 ); aKsi = G4UniformRand();

if ( aAQ == aQ )
{
 if ( aAQ != 3 )
 {
  NewCode = 111;                       // Pi0-meson
  if ( aKsi < 0.5 )
  {
   NewCode = 221;                     // Eta -meson
   if ( aKsi < 0.25 ) {NewCode = 331;} // Eta'-meson
  }
 } else
 {
  NewCode = 221;                      // Eta -meson
  if( aKsi < 0.5 ) {NewCode = 331;}    // Eta'-meson
 }
} else
{
 if ( aAQ > aQ ){ NewCode = aAQ*100 + aQ*10 + 1; NewCode *= aAQ/LeftAQ2; } 
 else           { NewCode = aQ*100 + aAQ*10 + 1; NewCode *=  aQ/LeftQ2;  }
}

TestParticle = G4ParticleTable::GetParticleTable()->FindParticle( NewCode );
if(!TestParticle) return false;
target->SetDefinition( TestParticle );
theParameters->SetTarMinDiffMass( 0.5 );     // Uzhi March 2016 ???
theParameters->SetTarMinNonDiffMass( 0.5 );
// Uzhi March 2016

      // Sampling kinematical properties
      // 1 string LeftAQ1-LeftQ1// 2 string LeftAQ2-LeftQ2
      G4ThreeVector Quark_Mom[4];
      G4double ModMom2[4]; //ModMom[4], 

      AveragePt2 = 200.0*200.0; maxPtSquare = S;

      G4double SumMt( 0.0 );
      G4double MassQ2 = 0.0; //100.0*100.0*MeV*MeV;
      G4int    NumberOfTries( 0 );
      G4double ScaleFactor( 1.0 );

      const G4int maxNumberOfLoops = 1000;
      G4int loopCounter = 0;
      do { 
        NumberOfTries++;
        if ( NumberOfTries == 100*(NumberOfTries/100) ) { 
          // At large number of tries it would be better to reduce the values of <Pt^2>
          ScaleFactor /= 2.0;
          AveragePt2 *= ScaleFactor;
        }
        G4ThreeVector PtSum( 0.0, 0.0, 0.0 );
        for( G4int i = 0; i < 4; i++ ) {
          Quark_Mom[i] = GaussianPt( AveragePt2, maxPtSquare );
          PtSum += Quark_Mom[i];
        }
        PtSum /= 4.0;
        SumMt = 0.0;    
        for ( G4int i = 0; i < 4; i++ ) {
          Quark_Mom[i] -= PtSum;
          //ModMom[i] = Quark_Mom[i].mag();
          ModMom2[i] = Quark_Mom[i].mag2();
          SumMt += std::sqrt( ModMom2[i] + MassQ2 );
        }
      } while ( ( SumMt > SqrtS ) &&  
                ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */
      if ( loopCounter >= maxNumberOfLoops ) {
        return false;
      }

      G4double WminusTarget( 0.0 ), WplusProjectile( 0.0 );

      // Sampling X's of anti-baryon
      G4double Alfa_R = 0.5;
      NumberOfTries = 0;
      ScaleFactor = 1.0;
      G4bool Succes( true );

      loopCounter = 0;
      do { 

        Succes = true;
        NumberOfTries++;
        if ( NumberOfTries == 100*(NumberOfTries/100) ) { 
          // At large number of tries it would be better to reduce the values of Pt's
          ScaleFactor /= 2.0;
        }

        if ( Alfa_R == 1.0 ) {
          G4double Xaq1 = std::sqrt( G4UniformRand() );
          G4double Xaq2 = 1.0 - Xaq1;
          Quark_Mom[0].setZ( Xaq1 ); Quark_Mom[1].setZ( Xaq2 ); 
        } else {
          G4double Xaq1 = sqr( std::sin( pi/2.0*G4UniformRand() ) );
          G4double Xaq2 = 1.0 - Xaq1;
          Quark_Mom[0].setZ( Xaq1 ); Quark_Mom[1].setZ( Xaq2 );
        }

        // Sampling X's of baryon ------------
        if ( Alfa_R == 1.0 ) {
          G4double Xq1 = 1.0 - std::sqrt( G4UniformRand() );
          G4double Xq2 = 1.0 - Xq1;
          Quark_Mom[2].setZ( Xq1 ); Quark_Mom[3].setZ( Xq2 );
        } else {
          G4double Xq1 = sqr( std::sin( pi/2.0*G4UniformRand() ) );
          G4double Xq2 = 1.0 - Xq1;
          Quark_Mom[2].setZ( Xq1 ); Quark_Mom[3].setZ( Xq2 );
        }

        G4double Alfa( 0.0 ), Beta( 0.0 );
        for ( G4int i = 0; i < 2; i++ ) {  // For Anti-baryon
          if ( Quark_Mom[i].getZ() != 0.0 ) {
            Alfa += ( ScaleFactor * ModMom2[i] + MassQ2 ) / Quark_Mom[i].getZ();
          } else {
            Succes = false;
          }
        } 
        for ( G4int i = 2; i < 4; i++ ) {  // For baryon
          if ( Quark_Mom[i].getZ() != 0.0 ) { 
            Beta += ( ScaleFactor * ModMom2[i] + MassQ2 ) / Quark_Mom[i].getZ();
          } else {
            Succes = false;
          }
        } 

        if ( ! Succes ) continue;

        if ( std::sqrt( Alfa ) + std::sqrt( Beta ) > SqrtS ) {
          Succes = false; 
          continue;
        }

        G4double DecayMomentum2 = S*S + Alfa*Alfa + Beta*Beta
                                - 2.0*S*Alfa - 2.0*S*Beta - 2.0*Alfa*Beta;
        WminusTarget = ( S - Alfa + Beta + std::sqrt( DecayMomentum2 ) ) / 2.0 / SqrtS; 
        WplusProjectile = SqrtS - Beta/WminusTarget;

      } while ( ( ! Succes ) &&
                ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */
      if ( loopCounter >= maxNumberOfLoops ) {
        return false;
      }

      G4double SqrtScaleF = std::sqrt( ScaleFactor );

      for ( G4int i = 0; i < 2; i++ ) {
        G4double Pz = WplusProjectile * Quark_Mom[i].getZ() / 2.0 -
                      ( ScaleFactor * ModMom2[i] + MassQ2 ) /
                      ( 2.0 * WplusProjectile * Quark_Mom[i].getZ() ); 
        Quark_Mom[i].setZ( Pz );
        if ( ScaleFactor != 1.0 ) {
          Quark_Mom[i].setX( SqrtScaleF * Quark_Mom[i].getX() ); 
          Quark_Mom[i].setY( SqrtScaleF * Quark_Mom[i].getY() );
        }
        //G4cout << "Anti Q " << i << " " << Quark_Mom[i] << G4endl;
      }
      for ( G4int i = 2; i < 4; i++ ) {
        G4double Pz = -WminusTarget * Quark_Mom[i].getZ() / 2.0 +
                      ( ScaleFactor * ModMom2[i] + MassQ2 ) /
                      ( 2.0 * WminusTarget * Quark_Mom[i].getZ() );
        Quark_Mom[i].setZ( Pz );
        if ( ScaleFactor != 1.0 ) {
          Quark_Mom[i].setX( SqrtScaleF * Quark_Mom[i].getX() ); 
          Quark_Mom[i].setY( SqrtScaleF * Quark_Mom[i].getY() );
        }
        //G4cout << "Bary Q " << i << " " << Quark_Mom[i] << G4endl;
      }
      //G4cout << "Sum AQ " << Quark_Mom[0] + Quark_Mom[1] << G4endl
      //       << "Sum Q  " << Quark_Mom[2] + Quark_Mom[3] << G4endl;

      G4ThreeVector tmp = Quark_Mom[0] + Quark_Mom[2];
      G4LorentzVector Pstring1( tmp, std::sqrt( Quark_Mom[0].mag2() + MassQ2 ) +
                                     std::sqrt( Quark_Mom[2].mag2() + MassQ2 ) );
      G4double Ystring1 = Pstring1.rapidity();

      //G4cout << "Mom 1 string " << G4endl << Quark_Mom[0] << G4endl << Quark_Mom[2] << G4endl
      //       << tmp << " " << tmp.mag() << G4endl;
      //G4cout << "1 str " << Pstring1 << " " << Pstring1.mag() << " " << Ystring1 << G4endl;

      tmp = Quark_Mom[1] + Quark_Mom[3];
      G4LorentzVector Pstring2( tmp, std::sqrt( Quark_Mom[1].mag2() + MassQ2 ) +
                                     std::sqrt( Quark_Mom[3].mag2() + MassQ2 ) );
      G4double Ystring2 = Pstring2.rapidity();

      //G4cout << "Mom 2 string " << G4endl <<Quark_Mom[1] << G4endl << Quark_Mom[3] << G4endl
      //       << tmp << " " << tmp.mag() << G4endl;
      //G4cout << "2 str " << Pstring2 << " " << Pstring2.mag() << " " << Ystring2 << G4endl;

      if ( Ystring1 > Ystring2 ) {
        Pprojectile = Pstring1;
        Ptarget     = Pstring2;
      } else {
        Pprojectile = Pstring2;
        Ptarget     = Pstring1;
      }

      //G4cout << "SumP CMS " << Pprojectile + Ptarget << " " << SqrtS << G4endl;
      Pprojectile.transform( toLab );
      Ptarget.transform( toLab );
      //G4cout << " SumP Lab " << Pprojectile + Ptarget << " " << SqrtS << G4endl;

      // Calculation of the creation time
      projectile->SetTimeOfCreation( target->GetTimeOfCreation() );
      projectile->SetPosition( target->GetPosition() );
      // Creation time and position of target nucleon were determined in
      // ReggeonCascade() of G4FTFModel
      //G4cout << "Mproj " << Pprojectile.mag() << G4endl << "Mtarg " << Ptarget.mag() << G4endl;
      projectile->Set4Momentum( Pprojectile );
      target->Set4Momentum( Ptarget );
      projectile->IncrementCollisionCount( 1 );
      target->IncrementCollisionCount( 1 );

theParameters->SetProbabilityOfAnnihilation( 0.0 );  // Uzhi March 2016

      return true;

    } // End of if ( CandidatsN != 0 )

  } // End of if ( Ksi < ( X_a + X_b + X_c ) / Xannihilation )

  // Simulation of anti-quark-quark string creation

  if ( Ksi < ( X_a + X_b + X_c + X_d ) / Xannihilation ) {

    #ifdef debugFTFannih 
    G4cout << "Process d, only 1 quark - anti-quark string" << G4endl;
    #endif

    G4int CandidatsN( 0 ), CandAQ[36], CandQ[36];
    G4int LeftAQ( 0 ), LeftQ( 0 );

    for ( G4int iAQ1 = 0; iAQ1 < 3; iAQ1++ ) {
      for ( G4int iAQ2 = 0; iAQ2 < 3; iAQ2++ ) {
        if ( iAQ1 != iAQ2 ) {
          for ( G4int iQ1 = 0; iQ1 < 3; iQ1++ ) {
            for ( G4int iQ2 = 0; iQ2 < 3; iQ2++ ) {
              if ( iQ1 != iQ2 ) {
                if ( -AQ[iAQ1] == Q[iQ1]  &&  -AQ[iAQ2] == Q[iQ2] ) {
                  if ( iAQ1 == 0  &&  iAQ2 == 1 ) { CandAQ[CandidatsN] = 2; }
                  if ( iAQ1 == 1  &&  iAQ2 == 0 ) { CandAQ[CandidatsN] = 2; }

                  if ( iAQ1 == 0  &&  iAQ2 == 2 ) { CandAQ[CandidatsN] = 1; }
                  if ( iAQ1 == 2  &&  iAQ2 == 0 ) { CandAQ[CandidatsN] = 1; }

                  if ( iAQ1 == 1  &&  iAQ2 == 2 ) { CandAQ[CandidatsN] = 0; }
                  if ( iAQ1 == 2  &&  iAQ2 == 1 ) { CandAQ[CandidatsN] = 0; }

                  if ( iQ1 == 0   &&   iQ2 == 1 ) { CandQ[CandidatsN]  = 2; }
                  if ( iQ1 == 1   &&   iQ2 == 0 ) { CandQ[CandidatsN]  = 2; }

                  if ( iQ1 == 0   &&   iQ2 == 2 ) { CandQ[CandidatsN]  = 1; }
                  if ( iQ1 == 2   &&   iQ2 == 0 ) { CandQ[CandidatsN]  = 1; }

                  if ( iQ1 == 1   &&   iQ2 == 2 ) { CandQ[CandidatsN]  = 0; }
                  if ( iQ1 == 2   &&   iQ2 == 1 ) { CandQ[CandidatsN]  = 0; }
                  CandidatsN++;
                }
              }
            }
          }
        }
      }
    }

    if ( CandidatsN != 0 ) {
      G4int SampledCase = G4RandFlat::shootInt( G4long( CandidatsN ) );
      LeftAQ = AQ[ CandAQ[SampledCase] ];
      LeftQ  =  Q[ CandQ[SampledCase] ];
      //G4cout << "Left Aq Q " << LeftAQ << " " << LeftQ << G4endl;

      // Set the string properties
      projectile->SplitUp();
      //projectile->SetFirstParton( LeftAQ );
      //projectile->SetSecondParton( LeftQ );
      projectile->SetFirstParton( LeftQ );
      projectile->SetSecondParton( LeftAQ );
      projectile->SetStatus( 0 );

// Uzhi March 2016 start
G4int aAQ, aQ;
aAQ=std::abs( LeftAQ ); aQ=std::abs( LeftQ );

G4int NewCode;
G4double aKsi = G4UniformRand();

if ( aAQ == aQ )
{
 if ( aAQ != 3 )
 {
  NewCode = 111;                       // Pi0-meson
  if ( aKsi < 0.5 )
  {
   NewCode = 221;                     // Eta -meson
   if ( aKsi < 0.25 ) {NewCode = 331;} // Eta'-meson
  }
 } else
 {
  NewCode = 221;                      // Eta -meson
  if( aKsi < 0.5 ) {NewCode = 331;}    // Eta'-meson
 }
} else
{
 if ( aAQ > aQ ){ NewCode = aAQ*100 + aQ*10 + 1; NewCode *= aAQ/LeftAQ; } 
 else           { NewCode = aQ*100 + aAQ*10 + 1; NewCode *=  aQ/LeftQ;  }
}

G4ParticleDefinition* TestParticle = G4ParticleTable::GetParticleTable()->FindParticle( NewCode );
if(!TestParticle) return false;
projectile->SetDefinition( TestParticle );
theParameters->SetProjMinDiffMass( 0.5 );            // (0.5)  // GeV Uzhi March 2016
theParameters->SetProjMinNonDiffMass( 0.5 );
// Uzhi March 2016 end

      target->SetStatus( 4 );  // The target nucleon has annihilated 3->4 Uzhi Oct 2014
      Pprojectile.setPx( 0.0 );
      Pprojectile.setPy( 0.0 );
      Pprojectile.setPz( 0.0 );
      Pprojectile.setE( SqrtS );
      Pprojectile.transform( toLab );

      // Calculation of the creation time
      projectile->SetTimeOfCreation( target->GetTimeOfCreation() );
      projectile->SetPosition( target->GetPosition() );
      // Creation time and position of target nucleon were determined in
      // ReggeonCascade() of G4FTFModel

      //G4cout << "Mproj " << Pprojectile.mag() << G4endl << "Mtarg " << Ptarget.mag() << G4endl;
      projectile->Set4Momentum( Pprojectile );

      projectile->IncrementCollisionCount( 1 );
      target->IncrementCollisionCount( 1 );

theParameters->SetProbabilityOfAnnihilation( 0.0 );  // Uzhi March 2016

      return true;
    }

  }  // End of if ( Ksi < ( X_a + X_b + X_c + X_d ) / Xannihilation )

  //G4cout << "Pr Y " << Pprojectile.rapidity() << " Tr Y " << Ptarget.rapidity() << G4endl;
  return true;
}


//============================================================================

G4double G4FTFAnnihilation::ChooseX( G4double /* Alpha */, G4double /* Beta */ ) const {
  // If for sampling Xs other values of Alfa and Beta instead of 0.5 will be
  // chosen the method will be implemented
  //G4double tmp = Alpha*Beta;
  //tmp *= 1.0;
  return 0.5;
}



//============================================================================

G4ThreeVector G4FTFAnnihilation::GaussianPt( G4double AveragePt2, G4double maxPtSquare ) const {
  //  @@ this method is used in FTFModel as well. Should go somewhere common!
  G4double Pt2( 0.0 );
  if ( AveragePt2 <= 0.0 ) {
    Pt2 = 0.0;
  } else {
    Pt2 = -AveragePt2 * G4Log( 1.0 + G4UniformRand() * 
                                        ( G4Exp( -maxPtSquare/AveragePt2 ) -1.0 ) );
  }
  G4double Pt = std::sqrt( Pt2 );
  G4double phi = G4UniformRand() * twopi;
  return G4ThreeVector ( Pt*std::cos( phi ), Pt*std::sin( phi ), 0.0 );
}


//============================================================================

void G4FTFAnnihilation::UnpackBaryon( G4int IdPDG, G4int& Q1, G4int& Q2, G4int& Q3 ) const {
  G4int AbsId = std::abs( IdPDG );
  Q1 =   AbsId          / 1000;
  Q2 = ( AbsId % 1000 ) / 100;
  Q3 = ( AbsId % 100 )  / 10;     
  if ( IdPDG < 0 ) { Q1 = -Q1; Q2 = -Q2; Q3 = -Q3; }  // Anti-baryon     
  return;
}


//============================================================================

G4FTFAnnihilation::G4FTFAnnihilation( const G4FTFAnnihilation& ) {
  throw G4HadronicException( __FILE__, __LINE__, 
                             "G4FTFAnnihilation copy contructor not meant to be called" );
}


//============================================================================

const G4FTFAnnihilation & G4FTFAnnihilation::operator=( const G4FTFAnnihilation& ) {
  throw G4HadronicException( __FILE__, __LINE__, 
                             "G4FTFAnnihilation = operator not meant to be called" ); 
}


//============================================================================

int G4FTFAnnihilation::operator==( const G4FTFAnnihilation& ) const {
  throw G4HadronicException( __FILE__, __LINE__, 
                             "G4FTFAnnihilation == operator not meant to be called" );
}


//============================================================================

int G4FTFAnnihilation::operator!=( const G4FTFAnnihilation& ) const {
  throw G4HadronicException( __FILE__, __LINE__, 
                             "G4DiffractiveExcitation != operator not meant to be called" );
}
