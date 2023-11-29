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

//#include "UZHI_diffraction.hh"

#include "G4ParticleTable.hh"

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

  #ifdef debugFTFannih 
  G4cout << "---------------------------- Annihilation----------------" << G4endl;
  #endif

  CommonVariables common;

  // Projectile parameters
  common.Pprojectile = projectile->Get4Momentum();
  G4int ProjectilePDGcode = projectile->GetDefinition()->GetPDGEncoding();
  if ( ProjectilePDGcode > 0 ) {
    target->SetStatus( 3 );
    return false;
  } 
  G4double M0projectile2 = common.Pprojectile.mag2();

  // Target parameters
  G4int TargetPDGcode = target->GetDefinition()->GetPDGEncoding();
  common.Ptarget = target->Get4Momentum();
  G4double M0target2 = common.Ptarget.mag2();

  #ifdef debugFTFannih
  G4cout << "PDG codes " << ProjectilePDGcode << " " << TargetPDGcode << G4endl
         << "Pprojec " << common.Pprojectile << " " << common.Pprojectile.mag() << G4endl
         << "Ptarget " << common.Ptarget    << " " << common.Ptarget.mag() << G4endl
         << "M0 proj target " << std::sqrt( M0projectile2 ) 
         << " " << std::sqrt( M0target2 ) << G4endl;
  #endif

  // Kinematical properties of the interactions
  G4LorentzVector Psum = common.Pprojectile + common.Ptarget;  // 4-momentum in CMS
  common.S = Psum.mag2(); 
  common.SqrtS = std::sqrt( common.S );
  #ifdef debugFTFannih
  G4cout << "Psum SqrtS S " << Psum << " " << common.SqrtS << " " << common.S << G4endl;
  #endif

  // Transform momenta to cms and then rotate parallel to z axis
  G4LorentzRotation toCms( -1*Psum.boostVector() );
  G4LorentzVector Ptmp( toCms*common.Pprojectile );
  toCms.rotateZ( -1*Ptmp.phi() );
  toCms.rotateY( -1*Ptmp.theta() );
  common.toLab = toCms.inverse();

  if ( G4UniformRand() <= G4Pow::GetInstance()->powA( 1880.0/common.SqrtS, 4.0 ) ) {
    common.RotateStrings = true;
    common.RandomRotation.rotateZ( 2.0*pi*G4UniformRand() );
    common.RandomRotation.rotateY( std::acos( 2.0*G4UniformRand() - 1.0 ) );
    common.RandomRotation.rotateZ( 2.0*pi*G4UniformRand() );
  }

  G4double MesonProdThreshold = projectile->GetDefinition()->GetPDGMass() +
                                target->GetDefinition()->GetPDGMass() +
                                ( 2.0*140.0 + 16.0 )*MeV;  // 2 Mpi + DeltaE
  G4double Prel2 = sqr(common.S) + sqr(M0projectile2) + sqr(M0target2)
                   - 2.0*( common.S*(M0projectile2 + M0target2) + M0projectile2*M0target2 );
  Prel2 /= common.S;
  G4double X_a = 0.0, X_b = 0.0, X_c = 0.0, X_d = 0.0;
  if ( Prel2 <= 0.0 ) {
    // Annihilation at rest! Values are copied from Parameters
    X_a = 625.1;    // mb  // 3-shirt diagram
    X_b =   0.0;    // mb  // anti-quark-quark annihilation
    X_c =  49.989;  // mb  // 2 Q-Qbar string creation
    X_d =   6.614;  // mb  // One Q-Qbar string
    #ifdef debugFTFannih 
    G4cout << "Annih at Rest X a b c d " << X_a << " " << X_b << " " << X_c << " " << X_d 
           << G4endl;
    #endif
  } else { // Annihilation in flight!
    G4double FlowF = 1.0 / std::sqrt( Prel2 )*GeV;
    // Process cross sections
    X_a = 25.0*FlowF;  // mb 3-shirt diagram
    if ( common.SqrtS < MesonProdThreshold ) {
      X_b = 3.13 + 140.0*G4Pow::GetInstance()->powA( ( MesonProdThreshold - common.SqrtS )/GeV, 2.5 );
    } else {
      X_b = 6.8*GeV / common.SqrtS;  // mb anti-quark-quark annihilation
    }
    if ( projectile->GetDefinition()->GetPDGMass() + target->GetDefinition()->GetPDGMass()
         > common.SqrtS ) {
      X_b = 0.0;
    }
    // This can be in an interaction of low energy anti-baryon with off-shell nuclear nucleon
    X_c = 2.0 * FlowF * sqr( projectile->GetDefinition()->GetPDGMass() +
                             target->GetDefinition()->GetPDGMass() ) / common.S; // mb re-arrangement of
                                                                                 // 2 quarks and 2 anti-quarks
    X_d = 23.3*GeV*GeV / common.S; // mb anti-quark-quark string creation
    #ifdef debugFTFannih 
    G4cout << "Annih in Flight X a b c d " << X_a << " " << X_b << " " << X_c << " " << X_d
           << G4endl << "SqrtS MesonProdThreshold " << common.SqrtS << " " << MesonProdThreshold
           << G4endl;
    #endif
  } 

  G4bool isUnknown = false;
  if ( TargetPDGcode == 2212  ||  TargetPDGcode == 2214 ) {  // Target proton or Delta+
    if        ( ProjectilePDGcode == -2212  ||  ProjectilePDGcode == -2214 ) {  // anti_proton  or anti_Delta+
      X_b *= 5.0; X_c *= 5.0; X_d *= 6.0;  // Pbar P
    } else if ( ProjectilePDGcode == -2112  ||  ProjectilePDGcode == -2114 ) {  // anti_neutron or anti_Delta0
      X_b *= 4.0; X_c *= 4.0; X_d *= 4.0;  // NeutrBar P
    } else if ( ProjectilePDGcode == -3122 ) {                                  // anti_Lambda (no anti_Lambda* in PDG)
      X_b *= 3.0; X_c *= 3.0; X_d *= 2.0;  // LambdaBar P
    } else if ( ProjectilePDGcode == -3112 ) {                                  // anti_Sigma- (no anti_Sigma*- in G4)
      X_b *= 2.0; X_c *= 2.0; X_d *= 0.0;  // Sigma-Bar P
    } else if ( ProjectilePDGcode == -3212 ) {                                  // anti_Sigma0 (no anti_Sigma*0 in G4)
      X_b *= 3.0; X_c *= 3.0; X_d *= 2.0;  // Sigma0Bar P
    } else if ( ProjectilePDGcode == -3222 ) {                                  // anti_Sigma+ (no anti_Sigma*+ in G4)
      X_b *= 4.0; X_c *= 4.0; X_d *= 2.0;  // Sigma+Bar P
    } else if ( ProjectilePDGcode == -3312 ) {                                  // anti_Xi-    (no anti_Xi*-    in G4)
      X_b *= 1.0; X_c *= 1.0; X_d *= 0.0;  // Xi-Bar P
    } else if ( ProjectilePDGcode == -3322 ) {                                  // anti_Xi0    (no anti_Xi*0    in G4)
      X_b *= 2.0; X_c *= 2.0; X_d *= 0.0;  // Xi0Bar P
    } else if ( ProjectilePDGcode == -3334 ) {                                  // anti_Omega- (no anti_Omega*- in PDG)
      X_b *= 0.0; X_c *= 0.0; X_d *= 0.0;  // Omega-Bar P
    } else {
      isUnknown = true;
    }
  } else if ( TargetPDGcode == 2112  ||  TargetPDGcode == 2114 ) {  // Target neutron or Delta0
    if        ( ProjectilePDGcode == -2212  ||  ProjectilePDGcode == -2214 ) {  // anti_proton  or anti_Delta+
      X_b *= 4.0; X_c *= 4.0; X_d *= 4.0;  // Pbar N
    } else if ( ProjectilePDGcode == -2112  ||  ProjectilePDGcode == -2114 ) {  // anti_neutron or anti_Delta0
      X_b *= 5.0; X_c *= 5.0; X_d *= 6.0;  // NeutrBar N
    } else if ( ProjectilePDGcode == -3122 ) {                                  // anti_Lambda (no anti_Lambda* in PDG)
      X_b *= 3.0; X_c *= 3.0; X_d *= 2.0;  // LambdaBar N
    } else if ( ProjectilePDGcode == -3112 ) {                                  // anti_Sigma- (no anti_Sigma*- in G4)
      X_b *= 4.0; X_c *= 4.0; X_d *= 2.0;  // Sigma-Bar N                   
    } else if ( ProjectilePDGcode == -3212 ) {                                  // anti_Sigma0 (no anti_Sigma*0 in G4)
      X_b *= 3.0; X_c *= 3.0; X_d *= 2.0;  // Sigma0Bar N
    } else if ( ProjectilePDGcode == -3222 ) {                                  // anti_Sigma+ (no anti_Sigma*+ in G4)
      X_b *= 2.0; X_c *= 2.0; X_d *= 0.0;  // Sigma+Bar N
    } else if ( ProjectilePDGcode == -3312 ) {                                  // anti_Xi-    (no anti_Xi*-    in G4)
      X_b *= 2.0; X_c *= 2.0; X_d *= 0.0;  // Xi-Bar N
    } else if ( ProjectilePDGcode == -3322 ) {                                  // anti_Xi0    (no anti_Xi*0    in G4)
      X_b *= 1.0; X_c *= 1.0; X_d *= 0.0;  // Xi0Bar N
    } else if ( ProjectilePDGcode == -3334 ) {                                  // anti_Omega- (no anti_Omega*- in PDG)
      X_b *= 0.0; X_c *= 0.0; X_d *= 0.0;  // Omega-Bar N
    } else {
      isUnknown = true;
    }
  } else {
    isUnknown = true;
  }
  if ( isUnknown ) {
    G4cout << "Unknown anti-baryon for FTF annihilation: PDGcodes - "
           << ProjectilePDGcode << " " << TargetPDGcode << G4endl;
  }

  #ifdef debugFTFannih 
  G4cout << "Annih Actual X a b c d " << X_a << " " << X_b << " " << X_c << " " << X_d << G4endl;
  #endif

  G4double Xannihilation = X_a + X_b + X_c + X_d;

  // Projectile unpacking
  UnpackBaryon( ProjectilePDGcode, common.AQ[0], common.AQ[1], common.AQ[2] );

  // Target unpacking
  UnpackBaryon( TargetPDGcode, common.Q[0], common.Q[1], common.Q[2] ); 

  G4double Ksi = G4UniformRand();

  if ( Ksi < X_a / Xannihilation ) {
    return Create3QuarkAntiQuarkStrings( projectile, target, AdditionalString, theParameters, common );
  }

  G4int resultCode = 99;
  if ( Ksi < (X_a + X_b) / Xannihilation ) {
    resultCode = Create1DiquarkAntiDiquarkString( projectile, target, common );
    if ( resultCode == 0 ) {
      return true;
    } else if ( resultCode == 99 ) {
      return false;
    }
  }

  if ( Ksi < ( X_a + X_b + X_c ) / Xannihilation ) {
    resultCode = Create2QuarkAntiQuarkStrings( projectile, target, theParameters, common );
    if ( resultCode == 0 ) {
      return true;
    } else if ( resultCode == 99 ) {
      return false;
    }
  }

  if ( Ksi < ( X_a + X_b + X_c + X_d ) / Xannihilation ) {
    return Create1QuarkAntiQuarkString( projectile, target, theParameters, common );
  }

  return true;
}


//-----------------------------------------------------------------------

G4bool G4FTFAnnihilation::
Create3QuarkAntiQuarkStrings( G4VSplitableHadron* projectile, 
                              G4VSplitableHadron* target,
                              G4VSplitableHadron*& AdditionalString,
                              G4FTFParameters* theParameters,
                              G4FTFAnnihilation::CommonVariables& common ) const {
  // Simulation of 3 anti-quark - quark strings creation

  #ifdef debugFTFannih 
  G4cout << "Process a, 3 shirt diagram" << G4endl;
  #endif

  // Sampling kinematical properties of quark. It can be done before string's creation

  const G4int maxNumberOfLoops = 1000;
  G4double MassQ2 = 0.0;               // Simplest case is considered with Mass_Q = 0.0
                                       // In principle, this must work with Mass_Q != 0.0
  G4double      Quark_Xs[6];
  G4ThreeVector Quark_Mom[6];

  G4double Alfa_R = 0.5;
  G4double AveragePt2 = 200.0*200.0, maxPtSquare = common.S;
  G4double ScaleFactor = 1.0;
  G4double Alfa = 0.0, Beta = 0.0;

  G4int NumberOfTries = 0, loopCounter = 0;

  do {
    // Sampling X's of anti-baryon and baryon
    G4double x1 = 0.0, x2 = 0.0, x3 = 0.0;
    G4double Product = 1.0;
    for ( G4int iCase = 0; iCase < 2; ++iCase ) {  // anti-baryon (1st case), baryon (2nd case)

      G4double r1 = G4UniformRand(), r2 = G4UniformRand();
      if ( Alfa_R == 1.0 ) {
        x1 = 1.0 - std::sqrt( r1 );
        x2 = (1.0 - x1) * r2;
      } else {
        x1 = sqr( r1 );
        x2 = (1.0 - x1) * sqr( std::sin( pi/2.0*r2 ) );
      }
      x3 = 1.0 - x1 - x2;

      G4int index = iCase*3;  // 0 for anti-baryon, 3 for baryon
      Quark_Xs[index] = x1; Quark_Xs[index+1] = x2; Quark_Xs[index+2] = x3;
      Product *= (x1*x2*x3);
    }

    if ( Product == 0.0 ) continue;    

    ++NumberOfTries;
    if ( NumberOfTries == 100*(NumberOfTries/100) ) {
      // After a large number of tries, it is better to reduce the values of <Pt^2>
      ScaleFactor /= 2.0;
      AveragePt2 *= ScaleFactor;
    }

    G4ThreeVector PtSum( 0.0, 0.0, 0.0 );
    for ( G4int i = 0; i < 6; ++i ) {
      Quark_Mom [i] = GaussianPt( AveragePt2, maxPtSquare );
      PtSum += Quark_Mom[i];
    }

    PtSum /= 6.0;
    Alfa = 0.0; Beta = 0.0;

    for ( G4int i = 0; i < 6; ++i ) {  // Loop over the quarks and (anti-)quarks
      Quark_Mom[i] -= PtSum;

      G4double val = ( Quark_Mom[i].mag2() + MassQ2 ) / Quark_Xs[i];
      if ( i < 3 ) {  // anti-baryon
        Alfa += val;
      } else {        // baryon (iCase == 1)
        Beta += val;
      }
    }

  } while ( ( std::sqrt( Alfa ) + std::sqrt( Beta ) > common.SqrtS ) &&
            ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */

  if ( loopCounter >= maxNumberOfLoops ) {
    return false;
  }

  G4double DecayMomentum2 = sqr(common.S) + sqr(Alfa) + sqr(Beta)
                            - 2.0*( common.S*(Alfa + Beta) + Alfa*Beta );

  G4double WminusTarget = 0.0, WplusProjectile = 0.0;
  WminusTarget = ( common.S - Alfa + Beta + std::sqrt( DecayMomentum2 ) ) / 2.0 / common.SqrtS; 
  WplusProjectile = common.SqrtS - Beta/WminusTarget;

  for ( G4int iCase = 0; iCase < 2; ++iCase ) {  // anti-baryon (1st case), baryon (2nd case)
    G4int index = iCase*3;                 // 0 for anti-baryon, 3 for baryon
    G4double w = WplusProjectile;          // for anti-baryon
    if ( iCase == 1 ) w = - WminusTarget;  // for baryon
    for ( G4int i = 0; i < 3; ++i ) {
      G4double Pz = w * Quark_Xs[index+i] / 2.0 -
                    ( Quark_Mom[index+i].mag2() + MassQ2 ) / 
                    ( 2.0 * w * Quark_Xs[index+i] ); 
      Quark_Mom[index+i].setZ( Pz );
    }
  }

  // Sampling of anti-quark order in projectile
  G4int SampledCase = (G4int)G4RandFlat::shootInt( 6 );
  G4int Tmp1 = 0, Tmp2 = 0;
  switch ( SampledCase ) {                                    
    case 1 : Tmp1 = common.AQ[1]; common.AQ[1] = common.AQ[2]; common.AQ[2] = Tmp1; break;
    case 2 : Tmp1 = common.AQ[0]; common.AQ[0] = common.AQ[1]; common.AQ[1] = Tmp1; break; 
    case 3 : Tmp1 = common.AQ[0]; Tmp2         = common.AQ[1]; common.AQ[0] = common.AQ[2]; 
             common.AQ[1] = Tmp1;         common.AQ[2] = Tmp2; break;
    case 4 : Tmp1 = common.AQ[0]; Tmp2         = common.AQ[1]; common.AQ[0] = Tmp2;
             common.AQ[1] = common.AQ[2]; common.AQ[2] = Tmp1; break;
    case 5 : Tmp1 = common.AQ[0]; Tmp2         = common.AQ[1]; common.AQ[0] = common.AQ[2];
             common.AQ[1] = Tmp2;         common.AQ[2] = Tmp1; break;
  }

  // Set the string properties
  // An anti quark - quark pair can have the quantum number of either a scalar meson
  // or a vector meson: the last digit of the PDG code is, respectively, 1 and 3. 
  // For simplicity only scalar is considered here.
  G4int NewCode = 0, antiQuark = 0, quark = 0;
  G4ParticleDefinition* TestParticle = nullptr;
  for ( G4int iString = 0; iString < 3; ++iString ) {  // Loop over the 3 string cases
    if ( iString == 0 ) {
      antiQuark = common.AQ[0];  quark = common.Q[0];
      projectile->SetFirstParton( antiQuark );
      projectile->SetSecondParton( quark );
      projectile->SetStatus( 0 );
    } else if ( iString == 1 ) {
      quark = common.Q[1];  antiQuark = common.AQ[1]; 
      target->SetFirstParton( quark );
      target->SetSecondParton( antiQuark );
      target->SetStatus( 0 );
    } else {  // iString == 2
      antiQuark = common.AQ[2]; quark =  common.Q[2];
    }
    G4int absAntiQuark = std::abs( antiQuark ), absQuark = std::abs( quark );
    G4double aKsi = G4UniformRand();
    if ( absAntiQuark == absQuark ) {
      if ( absAntiQuark != 3 ) {  // Not yet considered the case absAntiQuark 4 (charm) and 5 (bottom)
        NewCode = 111;            // Pi0-meson
        if ( aKsi < 0.5 ) {
          NewCode = 221;          // Eta -meson
          if ( aKsi < 0.25 ) {
            NewCode = 331;        // Eta'-meson
          }
        }
      } else {
        NewCode = 221;            // Eta -meson
        if ( aKsi < 0.5 ) {
          NewCode = 331;          // Eta'-meson
        }
      }
    } else {  // Vector mesons - rho, omega, phi (not yet considered the analogous cases for charm and bottom)
      if ( absAntiQuark > absQuark ) {
        NewCode = absAntiQuark*100 + absQuark*10 + 1; NewCode *= absAntiQuark/antiQuark;
      } else {
        NewCode = absQuark*100 + absAntiQuark*10 + 1; NewCode *= absQuark/quark;
      }
    }
    if ( iString == 2 ) AdditionalString = new G4DiffractiveSplitableHadron();
    TestParticle = G4ParticleTable::GetParticleTable()->FindParticle( NewCode );
    if ( ! TestParticle ) return false;
    if ( iString == 0 ) {
      projectile->SetDefinition( TestParticle );
      theParameters->SetProjMinDiffMass( 0.5 );     // 0.5 GeV : Min diffractive mass of pi-meson
      theParameters->SetProjMinNonDiffMass( 0.5 );  // It must be self-consistent with Parameters
    } else if ( iString == 1 ) {
      target->SetDefinition( TestParticle );
      theParameters->SetTarMinDiffMass( 0.5 );
      theParameters->SetTarMinNonDiffMass( 0.5 );
    } else {  // iString == 2
      AdditionalString->SetDefinition( TestParticle );
      AdditionalString->SetFirstParton( common.AQ[2] );
      AdditionalString->SetSecondParton( common.Q[2] );
      AdditionalString->SetStatus( 0 );
    }
  }  // End of the for loop over the 3 string cases

  // 1st string AQ[0]-Q[0], 2nd string AQ[1]-Q[1], 3rd string AQ[2]-Q[2]

  G4LorentzVector Pstring1, Pstring2, Pstring3;
  G4int QuarkOrder[3] = { 0 };
  G4double YstringMax = 0.0, YstringMin = 0.0;
  for ( G4int i = 0; i < 3; ++i ) {
    G4ThreeVector tmp = Quark_Mom[i] + Quark_Mom[i+3];
    G4LorentzVector Pstring( tmp, std::sqrt( Quark_Mom[i].mag2() + MassQ2 ) +
                                  std::sqrt( Quark_Mom[i+3].mag2() + MassQ2 ) );
    // Add protection for  rapidity = 0.5*ln( (E+Pz)/(E-Pz) )
    G4double Ystring = 0.0;
    if ( Pstring.e() > 1.0e-30 ) {
      if ( Pstring.e() + Pstring.pz() < 1.0e-30 ) {    // Very small numerator   in the logarithm 
        Ystring = -1.0e30;   // A very large negative value (E ~= -Pz)
        if ( Pstring.e() - Pstring.pz() < 1.0e-30 ) {  // Very small denominator in the logarithm
          Ystring = 1.0e30;  // A very large positive value (E ~= Pz)
        } else {  // Normal case
          Ystring = Pstring.rapidity();
        }
      }
    }
    // Keep ordering in rapidity: "1" highest, "2" middle, "3" smallest
    if ( i == 0 ) {
      Pstring1 = Pstring;     YstringMax = Ystring;
      QuarkOrder[0] = 0;
    } else if ( i == 1 ) {
      if ( Ystring > YstringMax ) {
        Pstring2 = Pstring1;  YstringMin = YstringMax;
        Pstring1 = Pstring;   YstringMax = Ystring;
        QuarkOrder[0] = 1; QuarkOrder[1] = 0;
      } else {
        Pstring2 = Pstring;   YstringMin = Ystring;
        QuarkOrder[1] = 1;
      }
    } else {  // i == 2
      if ( Ystring > YstringMax ) {
        Pstring3 = Pstring2;
        Pstring2 = Pstring1;
        Pstring1 = Pstring;
        QuarkOrder[1] = QuarkOrder[0];
        QuarkOrder[2] = QuarkOrder[1];
        QuarkOrder[0] = 2;
      } else if ( Ystring > YstringMin ) {
        Pstring3 = Pstring2;
        Pstring2 = Pstring;
      } else {
        Pstring3 = Pstring;
        QuarkOrder[2] = 2;
      }
    }
  }

  G4LorentzVector Quark_4Mom[6];
  for ( G4int i = 0; i < 6; ++i ) {
    Quark_4Mom[i] = G4LorentzVector( Quark_Mom[i], std::sqrt( Quark_Mom[i].mag2() + MassQ2 ) );
    if ( common.RotateStrings ) Quark_4Mom[i] *= common.RandomRotation;
    Quark_4Mom[i].transform( common.toLab );
  }

  projectile->Splitting();
  projectile->GetNextAntiParton()->Set4Momentum( Quark_4Mom[QuarkOrder[0]] );
  projectile->GetNextParton()->Set4Momentum( Quark_4Mom[QuarkOrder[0]+3] );

  target->Splitting();
  target->GetNextParton()->Set4Momentum( Quark_4Mom[QuarkOrder[2]] );
  target->GetNextAntiParton()->Set4Momentum( Quark_4Mom[QuarkOrder[2]+3] );

  AdditionalString->Splitting();
  AdditionalString->GetNextAntiParton()->Set4Momentum( Quark_4Mom[QuarkOrder[1]] );
  AdditionalString->GetNextParton()->Set4Momentum( Quark_4Mom[QuarkOrder[1]+3] );

  common.Pprojectile = Pstring1;           // Highest rapidity
  common.Ptarget     = Pstring3;           // Lowest rapidity
  G4LorentzVector LeftString( Pstring2 );  // Middle rapidity

  if ( common.RotateStrings ) {
    common.Pprojectile *= common.RandomRotation;
    common.Ptarget     *= common.RandomRotation;
    LeftString         *= common.RandomRotation;
  }

  common.Pprojectile.transform( common.toLab );
  common.Ptarget.transform( common.toLab );
  LeftString.transform( common.toLab );
  
  // Calculation of the creation time
  // Creation time and position of target nucleon were determined in ReggeonCascade() of G4FTFModel
  projectile->SetTimeOfCreation( target->GetTimeOfCreation() );
  projectile->SetPosition( target->GetPosition() );
  AdditionalString->SetTimeOfCreation( target->GetTimeOfCreation() );
  AdditionalString->SetPosition( target->GetPosition() );

  projectile->Set4Momentum( common.Pprojectile );
  AdditionalString->Set4Momentum( LeftString );
  target->Set4Momentum( common.Ptarget );

  projectile->IncrementCollisionCount( 1 );
  AdditionalString->IncrementCollisionCount( 1 );
  target->IncrementCollisionCount( 1 );

  return true;
}


//-----------------------------------------------------------------------

G4int G4FTFAnnihilation::
Create1DiquarkAntiDiquarkString( G4VSplitableHadron* projectile, 
                                 G4VSplitableHadron* target,
                                 G4FTFAnnihilation::CommonVariables& common ) const {
  // Simulation of anti-diquark-diquark string creation.
  // This method returns an integer code - instead of a boolean, with the following meaning:
  //   "0" : successfully ended and nothing else needs to be done;
  //   "1" : successfully completed, but the work needs to be continued;
  //  "99" : unsuccessfully ended, nothing else can be done.

  #ifdef debugFTFannih 
  G4cout << "Process b, quark - anti-quark annihilation, di-q - anti-di-q string" << G4endl;
  #endif

  G4int CandidatsN = 0, CandAQ[9][2] = {}, CandQ[9][2] = {};
  for ( G4int iAQ = 0; iAQ < 3; ++iAQ ) {  // index of the 3 constituent anti-quarks of the antibaryon projectile
    for ( G4int iQ = 0; iQ < 3; ++iQ ) {   // index of the 3 constituent quarks of the target nucleon
      if ( -common.AQ[iAQ] == common.Q[iQ] ) {  // antiquark - quark that can annihilate
        // Here "0", "1", "2" means, respectively, "first", "second" and "third" constituent
        // of the (anti-baryon) projectile or (nucleon) target.
        if ( iAQ == 0 ) { CandAQ[CandidatsN][0] = 1; CandAQ[CandidatsN][1] = 2; }
        if ( iAQ == 1 ) { CandAQ[CandidatsN][0] = 0; CandAQ[CandidatsN][1] = 2; }
        if ( iAQ == 2 ) { CandAQ[CandidatsN][0] = 0; CandAQ[CandidatsN][1] = 1; }
        if ( iQ  == 0 ) { CandQ[CandidatsN][0]  = 1; CandQ[CandidatsN][1]  = 2; }
        if ( iQ  == 1 ) { CandQ[CandidatsN][0]  = 0; CandQ[CandidatsN][1]  = 2; }
        if ( iQ  == 2 ) { CandQ[CandidatsN][0]  = 0; CandQ[CandidatsN][1]  = 1; }
        ++CandidatsN;
      }
    }
  }

  // Remaining two (anti-)quarks that form the (anti-)diquark
  G4int LeftAQ1 = 0, LeftAQ2 = 0, LeftQ1 = 0, LeftQ2 = 0;  
  if ( CandidatsN != 0 ) {
    G4int SampledCase = (G4int)G4RandFlat::shootInt( CandidatsN );
    LeftAQ1 = common.AQ[ CandAQ[SampledCase][0] ];
    LeftAQ2 = common.AQ[ CandAQ[SampledCase][1] ];
    LeftQ1  =  common.Q[ CandQ[SampledCase][0] ];
    LeftQ2  =  common.Q[ CandQ[SampledCase][1] ];

    // Build anti-diquark and diquark : the last digit can be either 3 - for all combinations
    // of anti-quark - anti-quark and quark - quark - or 1 - only when the two anti-quarks
    // or quarks are different. For simplicity, only 3 is considered.
    G4int Anti_DQ = 0, DQ = 0;
    if ( std::abs( LeftAQ1 ) > std::abs( LeftAQ2 ) ) { 
      Anti_DQ = 1000*LeftAQ1 + 100*LeftAQ2 - 3;
    } else {
      Anti_DQ = 1000*LeftAQ2 + 100*LeftAQ1 - 3;
    }
    if ( std::abs( LeftQ1 ) > std::abs( LeftQ2 ) ) { 
      DQ = 1000*LeftQ1 + 100*LeftQ2 + 3;
    } else {
      DQ = 1000*LeftQ2 + 100*LeftQ1 + 3;
    }

    // Set the string properties
    projectile->SetFirstParton( DQ );
    projectile->SetSecondParton( Anti_DQ );

    // It is assumed that quark and di-quark masses are 0.
    G4LorentzVector Pquark  = G4LorentzVector( 0.0, 0.0, -common.SqrtS/2.0, common.SqrtS/2.0 );
    G4LorentzVector Paquark = G4LorentzVector( 0.0, 0.0,  common.SqrtS/2.0, common.SqrtS/2.0 );

    if ( common.RotateStrings ) {
      Pquark  *= common.RandomRotation;
      Paquark *= common.RandomRotation;
    }

    Pquark.transform( common.toLab );
    Paquark.transform( common.toLab );

    projectile->GetNextParton()->Set4Momentum( Pquark );
    projectile->GetNextAntiParton()->Set4Momentum( Paquark );

    projectile->Splitting();

    projectile->SetStatus( 0 );
    target->SetStatus( 4 );  // The target nucleon has annihilated 3->4
    common.Pprojectile.setPx( 0.0 );
    common.Pprojectile.setPy( 0.0 );
    common.Pprojectile.setPz( 0.0 );
    common.Pprojectile.setE( common.SqrtS );
    common.Pprojectile.transform( common.toLab );

    // Calculation of the creation time
    // Creation time and position of target nucleon were determined in ReggeonCascade() of G4FTFModel
    projectile->SetTimeOfCreation( target->GetTimeOfCreation() );
    projectile->SetPosition( target->GetPosition() );
    projectile->Set4Momentum( common.Pprojectile );

    projectile->IncrementCollisionCount( 1 );
    target->IncrementCollisionCount( 1 );

    return 0;  // Completed successfully: nothing else to be done
  }  // End of if ( CandidatsN != 0 )

  // If we allow the string to interact with other nuclear nucleons, we have to
  // set up MinDiffrMass in Parameters, and ascribe a PDGEncoding. To be done yet!
  
  return 1;  // Successfully ended, but the work is not over
}


//-----------------------------------------------------------------------

G4int G4FTFAnnihilation::
Create2QuarkAntiQuarkStrings( G4VSplitableHadron* projectile,
                              G4VSplitableHadron* target,
                              G4FTFParameters* theParameters,
                              G4FTFAnnihilation::CommonVariables& common ) const {
  // Simulation of 2 anti-quark-quark strings creation.
  // This method returns an integer code - instead of a boolean, with the following meaning:
  //   "0" : successfully ended and nothing else needs to be done;
  //   "1" : successfully completed, but the work needs to be continued;
  //  "99" : unsuccessfully ended, nothing else can be done.

  #ifdef debugFTFannih 
  G4cout << "Process c, quark - anti-quark and string junctions annihilation, 2 strings left."
         << G4endl;
  #endif

  // Sampling kinematical properties: 1st string LeftAQ1-LeftQ1, 2nd string LeftAQ2-LeftQ2
  G4ThreeVector Quark_Mom[4];
  G4double Quark_Xs[4];
  G4double AveragePt2 = 200.0*200.0, maxPtSquare = common.S, MassQ2 = 0.0, ScaleFactor = 1.0;
  G4int NumberOfTries = 0, loopCounter = 0;
  const G4int maxNumberOfLoops = 1000;
  G4double Alfa = 0.0, Beta = 0.0;
  G4double WminusTarget = 0.0, WplusProjectile = 0.0, Alfa_R = 0.5;
  do { 
    // Sampling X's of the 2 quarks and 2 anti-quarks

    G4double Product = 1.0;
    for ( G4int iCase = 0; iCase < 2; ++iCase ) {  // Loop over the two strings
      G4double x = 0.0, r = G4UniformRand();
      if ( Alfa_R == 1.0 ) {
        if ( iCase == 0 ) {  // first string
          x = std::sqrt( r );
        } else {             // second string
          x = 1.0 - std::sqrt( r );
        }
      } else {
        x = sqr( std::sin( pi/2.0*r ) );
      }
      G4int index = iCase*2;  // 0 for the first string, 2 for the second string
      Quark_Xs[index] = x ; Quark_Xs[index+1] = 1.0 - x ;  
      Product *= x*(1.0-x);
    }

    if ( Product == 0.0 ) continue;

    ++NumberOfTries;
    if ( NumberOfTries == 100*(NumberOfTries/100) ) { 
      // After a large number of tries, it is better to reduce the values of <Pt^2>
      ScaleFactor /= 2.0;
      AveragePt2 *= ScaleFactor;
    }

    G4ThreeVector PtSum( 0.0, 0.0, 0.0 );
    for( G4int i = 0; i < 4; ++i ) {
      Quark_Mom[i] = GaussianPt( AveragePt2, maxPtSquare );
      PtSum += Quark_Mom[i];
    }

    PtSum /= 4.0;   
    for ( G4int i = 0; i < 4; ++i ) {
      Quark_Mom[i] -= PtSum;
    }

    Alfa = 0.0; Beta = 0.0;
    for ( G4int iCase = 0; iCase < 2; ++iCase ) {
       G4int index = iCase * 2;
       for ( G4int i = 0; i < 2; ++i ) {
          G4double val = ( Quark_Mom[index+i].mag2() + MassQ2 ) / Quark_Xs[index+i];
          if ( iCase == 0 ) {  // first string
            Alfa += val;
          } else {             // second string
            Beta += val;
          }
       }
    }

  } while ( ( std::sqrt( Alfa ) + std::sqrt( Beta ) > common.SqrtS ) &&  
            ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */

  if ( loopCounter >= maxNumberOfLoops ) {
    return 99;  // unsuccessfully ended, nothing else can be done
  }

  G4double DecayMomentum2 = sqr(common.S) + sqr(Alfa) + sqr(Beta)
                            - 2.0*( common.S*(Alfa + Beta) + Alfa*Beta );
  WminusTarget = ( common.S - Alfa + Beta + std::sqrt( DecayMomentum2 ) ) / 2.0 / common.SqrtS; 
  WplusProjectile = common.SqrtS - Beta/WminusTarget;

  for ( G4int iCase = 0; iCase < 2; ++iCase ) {  // Loop over the two strings
    G4int index = iCase*2;  // 0 for the first string, 2 for the second string
    for ( G4int i = 0; i < 2; ++i ) {
      G4double w = WplusProjectile;          // For the first string
      if ( iCase == 1 ) w = - WminusTarget;  // For the second string
      G4double Pz = w * Quark_Xs[index+i] / 2.0
                    - ( Quark_Mom[index+i].mag2() + MassQ2 ) /
                      ( 2.0 * w * Quark_Xs[index+i] ); 
      Quark_Mom[index+i].setZ( Pz );
    }
  }

  G4int CandidatsN = 0, CandAQ[9][2] = {}, CandQ[9][2] = {};
  G4int LeftAQ1 = 0, LeftAQ2 = 0, LeftQ1 = 0, LeftQ2 = 0;
  for ( G4int iAQ = 0; iAQ < 3; ++iAQ ) {  // index of the 3 constituent anti-quarks of the antibaryon projectile
    for ( G4int iQ = 0; iQ < 3; ++iQ ) {   // index of the 3 constituent quarks of the nucleon target
      if ( -common.AQ[iAQ] == common.Q[iQ] ) {  // antiquark - quark that can annihilate
        // Here "0", "1", "2" means, respectively, "first", "second" and "third" constituent
        // of the (anti-baryon) projectile or (nucleon) target.
        if ( iAQ == 0 ) { CandAQ[CandidatsN][0] = 1; CandAQ[CandidatsN][1] = 2; }
        if ( iAQ == 1 ) { CandAQ[CandidatsN][0] = 0; CandAQ[CandidatsN][1] = 2; }
        if ( iAQ == 2 ) { CandAQ[CandidatsN][0] = 0; CandAQ[CandidatsN][1] = 1; }
        if ( iQ  == 0 ) { CandQ[CandidatsN][0]  = 1; CandQ[CandidatsN][1]  = 2; }
        if ( iQ  == 1 ) { CandQ[CandidatsN][0]  = 0; CandQ[CandidatsN][1]  = 2; }
        if ( iQ  == 2 ) { CandQ[CandidatsN][0]  = 0; CandQ[CandidatsN][1]  = 1; }
        ++CandidatsN;
      }
    }
  }

  if ( CandidatsN != 0 ) {
    G4int SampledCase = (G4int)G4RandFlat::shootInt( CandidatsN );
    LeftAQ1 = common.AQ[ CandAQ[SampledCase][0] ];
    LeftAQ2 = common.AQ[ CandAQ[SampledCase][1] ];
    if ( G4UniformRand() < 0.5 ) {
      LeftQ1 = common.Q[ CandQ[SampledCase][0] ];
      LeftQ2 = common.Q[ CandQ[SampledCase][1] ];
    } else {
      LeftQ2 = common.Q[ CandQ[SampledCase][0] ];
      LeftQ1 = common.Q[ CandQ[SampledCase][1] ];
    }

    // Set the string properties
    // An anti quark - quark pair can have the quantum number of either a scalar meson
    // or a vector meson: the last digit of the PDG code is, respectively, 1 and 3. 
    // For simplicity only scalar is considered here.
    G4int NewCode = 0, antiQuark = 0, quark = 0;
    G4ParticleDefinition* TestParticle = nullptr;
    for ( G4int iString = 0; iString < 2; ++iString ) {  // Loop over the 2 string cases
      if ( iString == 0 ) {
        antiQuark = LeftAQ1; quark = LeftQ1;
        projectile->SetFirstParton( antiQuark );
        projectile->SetSecondParton( quark );
        projectile->SetStatus( 0 );
      } else {  // iString == 1
        quark = LeftQ2; antiQuark = LeftAQ2;
        target->SetFirstParton( quark );
        target->SetSecondParton( antiQuark );
        target->SetStatus( 0 );
      }
      G4int absAntiQuark = std::abs( antiQuark ), absQuark = std::abs( quark );
      G4double aKsi = G4UniformRand();
      if ( absAntiQuark == absQuark ) {
        if ( absAntiQuark != 3 ) {
          NewCode = 111;          // Pi0-meson
          if ( aKsi < 0.5 ) {
            NewCode = 221;        // Eta -meson
            if ( aKsi < 0.25 ) {
              NewCode = 331;      // Eta'-meson
            }
          }
        } else {
          NewCode = 221;          // Eta -meson
          if ( aKsi < 0.5 ) {
            NewCode = 331;        // Eta'-meson
          }
        }
      } else {
        if ( absAntiQuark > absQuark ) { 
          NewCode = absAntiQuark*100 + absQuark*10 + 1; NewCode *= absAntiQuark/antiQuark; 
        } else { 
          NewCode = absQuark*100 + absAntiQuark*10 + 1; NewCode *= absQuark/quark;
        }
      }
      TestParticle = G4ParticleTable::GetParticleTable()->FindParticle( NewCode );
      if ( ! TestParticle ) return 99;  // unsuccessfully ended, nothing else can be done
      if ( iString == 0 ) {
        projectile->SetDefinition( TestParticle );
        theParameters->SetProjMinDiffMass( 0.5 );
        theParameters->SetProjMinNonDiffMass( 0.5 );
      } else {  // iString == 1
        target->SetDefinition( TestParticle );
        theParameters->SetTarMinDiffMass( 0.5 );
        theParameters->SetTarMinNonDiffMass( 0.5 );
      }
    }  // End of loop over the 2 string cases

    G4int QuarkOrder[2];
    G4LorentzVector Pstring1, Pstring2;
    G4double Ystring1 = 0.0, Ystring2 = 0.0;

    for ( G4int iCase = 0; iCase < 2; ++iCase ) {  // Loop over the two strings
      G4ThreeVector tmp = Quark_Mom[iCase] + Quark_Mom[iCase+2];
      G4LorentzVector Pstring( tmp, std::sqrt( Quark_Mom[iCase].mag2()   + MassQ2 ) +
                                    std::sqrt( Quark_Mom[iCase+2].mag2() + MassQ2 ) );
      // Add protection for  rapidity = 0.5*ln( (E+Pz)/(E-Pz) )
      G4double Ystring = 0.0;
      if ( Pstring.e() > 1.0e-30 ) {
        if ( Pstring.e() + Pstring.pz() < 1.0e-30 ) {    // Very small numerator   in the logarithm 
          Ystring = -1.0e30;   // A very large negative value (E ~= -Pz)
          if ( Pstring.e() - Pstring.pz() < 1.0e-30 ) {  // Very small denominator in the logarithm
            Ystring = 1.0e30;  // A very large positive value (E ~= Pz)
          } else {  // Normal case
            Ystring = Pstring.rapidity();
          }
        }
      }
      if ( iCase == 0 ) {  // For the first string
        Pstring1 = Pstring; Ystring1 = Ystring;
      } else {             // For the second string
        Pstring2 = Pstring; Ystring2 = Ystring;          
      }
    }       
    if ( Ystring1 > Ystring2 ) {
      common.Pprojectile = Pstring1;  common.Ptarget = Pstring2;
      QuarkOrder[0] = 0; QuarkOrder[1] = 1;
    } else {
      common.Pprojectile = Pstring2;  common.Ptarget = Pstring1;
      QuarkOrder[0] = 1; QuarkOrder[1] = 0;
    }

    if ( common.RotateStrings ) {
      common.Pprojectile *= common.RandomRotation;
      common.Ptarget     *= common.RandomRotation;
    }

    common.Pprojectile.transform( common.toLab );
    common.Ptarget.transform( common.toLab );
    
    G4LorentzVector Quark_4Mom[4];
    for ( G4int i = 0; i < 4; ++i ) {
      Quark_4Mom[i] = G4LorentzVector( Quark_Mom[i], std::sqrt( Quark_Mom[i].mag2() + MassQ2 ) );
      if ( common.RotateStrings ) Quark_4Mom[i] *= common.RandomRotation;
      Quark_4Mom[i].transform( common.toLab );
    }

    projectile->Splitting();
    projectile->GetNextAntiParton()->Set4Momentum( Quark_4Mom[QuarkOrder[0]] );
    projectile->GetNextParton()->Set4Momentum( Quark_4Mom[QuarkOrder[0]+2] );

    target->Splitting();
    target->GetNextParton()->Set4Momentum( Quark_4Mom[QuarkOrder[1]] );
    target->GetNextAntiParton()->Set4Momentum( Quark_4Mom[QuarkOrder[1]+2] );

    // Calculation of the creation time
    // Creation time and position of target nucleon were determined in ReggeonCascade() of G4FTFModel
    projectile->SetTimeOfCreation( target->GetTimeOfCreation() );
    projectile->SetPosition( target->GetPosition() );
    projectile->Set4Momentum( common.Pprojectile );
    target->Set4Momentum( common.Ptarget );

    projectile->IncrementCollisionCount( 1 );
    target->IncrementCollisionCount( 1 );

    return 0;  // Completed successfully: nothing else to be done
  }  // End of if ( CandidatsN != 0 )

  return 1;  // Successfully ended, but the work is not over
}


//-----------------------------------------------------------------------

G4bool G4FTFAnnihilation::
Create1QuarkAntiQuarkString( G4VSplitableHadron* projectile,
                             G4VSplitableHadron* target,
                             G4FTFParameters* theParameters,
                             G4FTFAnnihilation::CommonVariables& common ) const {
  // Simulation of anti-quark - quark string creation

  #ifdef debugFTFannih 
  G4cout << "Process d, only 1 quark - anti-quark string" << G4endl;
  #endif

  // Determine the set of candidates anti-quark - quark pairs that do not annihilate.
  // Here "0", "1", "2" means, respectively, "first", "second" and "third" constituent
  // of the (anti-baryon) projectile or (nucleon) target.
  G4int CandidatsN = 0, CandAQ[36], CandQ[36];
  G4int LeftAQ = 0, LeftQ = 0;
  for ( G4int iAQ1 = 0; iAQ1 < 3; ++iAQ1 ) {
    for ( G4int iAQ2 = 0; iAQ2 < 3; ++iAQ2 ) {
      if ( iAQ1 != iAQ2 ) {
        for ( G4int iQ1 = 0; iQ1 < 3; ++iQ1 ) {
          for ( G4int iQ2 = 0; iQ2 < 3; ++iQ2 ) {
            if ( iQ1 != iQ2 ) {
              if ( -common.AQ[iAQ1] == common.Q[iQ1]  &&  -common.AQ[iAQ2] == common.Q[iQ2] ) {
                if        ( ( iAQ1 == 0  &&  iAQ2 == 1 ) || ( iAQ1 == 1  &&  iAQ2 == 0 ) ) { 
                  CandAQ[CandidatsN] = 2; 
                } else if ( ( iAQ1 == 0  &&  iAQ2 == 2 ) || ( iAQ1 == 2  &&  iAQ2 == 0 ) ) { 
                  CandAQ[CandidatsN] = 1;
                } else if ( ( iAQ1 == 1  &&  iAQ2 == 2 ) || ( iAQ1 == 2  &&  iAQ2 == 1 ) ) {
                  CandAQ[CandidatsN] = 0; 
                }
                if        ( ( iQ1 == 0   &&   iQ2 == 1 ) || ( iQ1 == 1   &&   iQ2 == 0 ) ) {
                  CandQ[CandidatsN]  = 2;
                } else if ( ( iQ1 == 0   &&   iQ2 == 2 ) || ( iQ1 == 2   &&   iQ2 == 0 ) ) {
                  CandQ[CandidatsN]  = 1;
                } else if ( ( iQ1 == 1   &&   iQ2 == 2 ) || ( iQ1 == 2   &&   iQ2 == 1 ) ) {
                  CandQ[CandidatsN]  = 0;
                }
                ++CandidatsN;
              }
            }
          }
        }
      }
    }
  }

  if ( CandidatsN != 0 ) {
    G4int SampledCase = (G4int)G4RandFlat::shootInt( CandidatsN );
    LeftAQ = common.AQ[ CandAQ[SampledCase] ];
    LeftQ  =  common.Q[ CandQ[SampledCase] ];

    // Set the string properties
    projectile->SetFirstParton( LeftQ );
    projectile->SetSecondParton( LeftAQ );
    projectile->SetStatus( 0 );
    G4int aAQ = std::abs( LeftAQ ), aQ = std::abs( LeftQ );
    G4int NewCode = 0;
    G4double aKsi = G4UniformRand();
    // The string can have the quantum number of either a scalar or a vector (whose last digit
    // of the PDG code is, respectively, 1 and 3). For simplicity only scalar is considered here.
    if ( aAQ == aQ ) {
      if ( aAQ != 3 ) {
        NewCode = 111;          // Pi0-meson
        if ( aKsi < 0.5 ) {
          NewCode = 221;        // Eta -meson
          if ( aKsi < 0.25 ) {
            NewCode = 331;      // Eta'-meson
          }
        }
      } else {
        NewCode = 221;          // Eta -meson
        if ( aKsi < 0.5 ) {
          NewCode = 331;        // Eta'-meson
        }
      }
    } else {
      if ( aAQ > aQ ) { 
        NewCode = aAQ*100 + aQ*10 + 1; NewCode *= aAQ/LeftAQ;
      } else { 
        NewCode = aQ*100 + aAQ*10 + 1; NewCode *=  aQ/LeftQ;
      }
    }

    G4ParticleDefinition* TestParticle = G4ParticleTable::GetParticleTable()->FindParticle( NewCode );
    if ( ! TestParticle ) return false;
    projectile->SetDefinition( TestParticle );
    theParameters->SetProjMinDiffMass( 0.5 );
    theParameters->SetProjMinNonDiffMass( 0.5 );

    target->SetStatus( 4 );  // The target nucleon has annihilated 3->4
    common.Pprojectile.setPx( 0.0 );
    common.Pprojectile.setPy( 0.0 );
    common.Pprojectile.setPz( 0.0 );
    common.Pprojectile.setE( common.SqrtS );

    common.Pprojectile.transform( common.toLab );

    G4LorentzVector Pquark  = G4LorentzVector( 0.0, 0.0, -common.SqrtS/2.0, common.SqrtS/2.0 );
    G4LorentzVector Paquark = G4LorentzVector( 0.0, 0.0, +common.SqrtS/2.0, common.SqrtS/2.0 );

    if ( common.RotateStrings ) { 
      Pquark *= common.RandomRotation; Paquark *= common.RandomRotation;
    }  
    Pquark.transform(common.toLab);  projectile->GetNextParton()->Set4Momentum(Pquark);
    Paquark.transform(common.toLab); projectile->GetNextAntiParton()->Set4Momentum(Paquark);

    projectile->Splitting();

    // Calculation of the creation time
    // Creation time and position of target nucleon were determined in ReggeonCascade() of G4FTFModel
    projectile->SetTimeOfCreation( target->GetTimeOfCreation() );
    projectile->SetPosition( target->GetPosition() );
    projectile->Set4Momentum( common.Pprojectile );

    projectile->IncrementCollisionCount( 1 );
    target->IncrementCollisionCount( 1 );

    return true;
  }  // End of if ( CandidatsN != 0 )

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
  G4double Pt2 = 0.0;
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
                             "G4FTFAnnihilation copy constructor not meant to be called" );
}


//============================================================================

const G4FTFAnnihilation & G4FTFAnnihilation::operator=( const G4FTFAnnihilation& ) {
  throw G4HadronicException( __FILE__, __LINE__, 
                             "G4FTFAnnihilation = operator not meant to be called" ); 
}


//============================================================================

G4bool G4FTFAnnihilation::operator==( const G4FTFAnnihilation& ) const {
  throw G4HadronicException( __FILE__, __LINE__, 
                             "G4FTFAnnihilation == operator not meant to be called" );
}


//============================================================================

G4bool G4FTFAnnihilation::operator!=( const G4FTFAnnihilation& ) const {
  throw G4HadronicException( __FILE__, __LINE__, 
                             "G4DiffractiveExcitation != operator not meant to be called" );
}
