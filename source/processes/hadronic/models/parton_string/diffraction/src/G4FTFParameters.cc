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
// $Id: G4FTFParameters.cc 102029 2016-12-16 14:53:08Z gcosmo $
// GEANT4 tag $Name:  $
//

#include <utility>                                        

#include "G4FTFParameters.hh"

#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleDefinition.hh"

#include "G4Proton.hh"
#include "G4Neutron.hh"

#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"

#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

//============================================================================

//#define debugFTFparams


//============================================================================

G4FTFParameters::G4FTFParameters() :
  FTFhNcmsEnergy( 0.0 ), 
  FTFxsManager( 0 ),
  FTFXtotal( 0.0 ), FTFXelastic( 0.0 ), FTFXinelastic( 0.0 ), FTFXannihilation( 0.0 ),
  ProbabilityOfAnnihilation( 0.0 ), ProbabilityOfElasticScatt( 0.0 ),
  RadiusOfHNinteractions2( 0.0 ), FTFSlope( 0.0 ), 
  AvaragePt2ofElasticScattering( 0.0 ), FTFGamma0( 0.0 ),
  DeltaProbAtQuarkExchange( 0.0 ), ProbOfSameQuarkExchange( 0.0 ), 
  ProjMinDiffMass( 0.0 ), ProjMinNonDiffMass( 0.0 ), ProbLogDistrPrD(0.0),
  TarMinDiffMass( 0.0 ), TarMinNonDiffMass( 0.0 ),
  AveragePt2( 0.0 ), ProbLogDistr( 0.0 ),
  Pt2kink( 0.0 ),
  MaxNumberOfCollisions( 0.0 ), ProbOfInelInteraction( 0.0 ), 
  CofNuclearDestructionPr( 0.0 ), CofNuclearDestruction( 0.0 ),
  R2ofNuclearDestruction( 0.0 ), ExcitationEnergyPerWoundedNucleon( 0.0 ),
  DofNuclearDestruction( 0.0 ), Pt2ofNuclearDestruction( 0.0 ), MaxPt2ofNuclearDestruction( 0.0 ) 
{
  for ( G4int i = 0; i < 4; i++ ) {
    for ( G4int j = 0; j < 7; j++ ) {
      ProcParams[i][j] = 0.0;
    }
  }
}


//============================================================================

G4FTFParameters::~G4FTFParameters() {}


//============================================================================

G4ThreadLocal bool G4FTFParameters::chipsComponentXSisInitialized = false;
G4ThreadLocal G4ChipsComponentXS* G4FTFParameters::chipsComponentXSinstance = 0;


//============================================================================

G4FTFParameters::G4FTFParameters( const G4ParticleDefinition* particle, 
                                  G4int theA, G4int theZ, G4double PlabPerParticle ) :
  FTFhNcmsEnergy( 0.0 ), 
  FTFxsManager( 0 ),
  FTFXtotal( 0.0 ), FTFXelastic( 0.0 ), FTFXinelastic( 0.0 ), FTFXannihilation( 0.0 ),
  ProbabilityOfAnnihilation( 0.0 ), ProbabilityOfElasticScatt( 0.0 ),
  RadiusOfHNinteractions2( 0.0 ), FTFSlope( 0.0 ), 
  AvaragePt2ofElasticScattering( 0.0 ), FTFGamma0( 0.0 ),
  DeltaProbAtQuarkExchange( 0.0 ), ProbOfSameQuarkExchange( 0.0 ), 
  ProjMinDiffMass( 0.0 ), ProjMinNonDiffMass( 0.0 ), ProbLogDistrPrD(0.0),
  TarMinDiffMass( 0.0 ), TarMinNonDiffMass( 0.0 ),
  AveragePt2( 0.0 ), ProbLogDistr( 0.0 ),
  Pt2kink( 0.0 ),
  MaxNumberOfCollisions( 0.0 ), ProbOfInelInteraction( 0.0 ), 
  CofNuclearDestructionPr( 0.0 ), CofNuclearDestruction( 0.0 ),
  R2ofNuclearDestruction( 0.0 ), ExcitationEnergyPerWoundedNucleon( 0.0 ),
  DofNuclearDestruction( 0.0 ), Pt2ofNuclearDestruction( 0.0 ), MaxPt2ofNuclearDestruction( 0.0 ) 
{
  for ( G4int i = 0; i < 4; i++ ) {
    for ( G4int j = 0; j < 7; j++ ) {
      ProcParams[i][j] = 0.0;
    }
  }

  G4int    ProjectilePDGcode    = particle->GetPDGEncoding();
  G4int    ProjectileabsPDGcode = std::abs( ProjectilePDGcode );
  G4double ProjectileMass       = particle->GetPDGMass();
  G4double ProjectileMass2      = ProjectileMass * ProjectileMass;

  G4int ProjectileBaryonNumber( 0 ), AbsProjectileBaryonNumber( 0 ), AbsProjectileCharge( 0 );
  G4bool ProjectileIsNucleus = false;

  if ( std::abs( particle->GetBaryonNumber() ) > 1 ) {  // The projectile is a nucleus
    ProjectileIsNucleus       = true;
    ProjectileBaryonNumber    = particle->GetBaryonNumber();
    AbsProjectileBaryonNumber = std::abs( ProjectileBaryonNumber );
    AbsProjectileCharge       = G4int( particle->GetPDGCharge() );
    if ( ProjectileBaryonNumber > 1 ) {
      ProjectilePDGcode = 2212; ProjectileabsPDGcode = 2212;   // Proton
    } else { 
      ProjectilePDGcode = -2212; ProjectileabsPDGcode = 2212;  // Anti-Proton
    }
    ProjectileMass  = G4Proton::Proton()->GetPDGMass();
    ProjectileMass2 = sqr( ProjectileMass );
  } 

  G4double TargetMass  = G4Proton::Proton()->GetPDGMass();
  G4double TargetMass2 = TargetMass * TargetMass;

  G4double Plab = PlabPerParticle;
  G4double Elab = std::sqrt( Plab*Plab + ProjectileMass2 );
  G4double KineticEnergy = Elab - ProjectileMass;

  G4double S = ProjectileMass2 + TargetMass2 + 2.0*TargetMass*Elab;

  #ifdef debugFTFparams
  G4cout << "--------- FTF Parameters --------------" << G4endl << "Proj Plab " 
         << ProjectilePDGcode << " " << Plab << G4endl << "Mass KinE " << ProjectileMass
         << " " << KineticEnergy << G4endl << " A Z " << theA << " " << theZ << G4endl;
  #endif

  G4double Ylab, Xtotal, Xelastic, Xannihilation;
  G4int NumberOfTargetNucleons;

  Ylab = 0.5 * G4Log( (Elab + Plab)/(Elab - Plab) );

  G4double ECMSsqr = S/GeV/GeV;
  G4double SqrtS   = std::sqrt( S )/GeV;

  #ifdef debugFTFparams
  G4cout << "Sqrt(s) " << SqrtS << G4endl;
  #endif

  TargetMass     /= GeV; TargetMass2     /= (GeV*GeV);
  ProjectileMass /= GeV; ProjectileMass2 /= (GeV*GeV);

  // Andrea Dotti (13Jan2013):
  // The following lines are changed for G4MT. Originally the code was:
  //   static G4ChipsComponentXS* _instance = new G4ChipsComponentXS();  // Witek Pokorski
  // Note the code could go back at original if _instance could be shared among threads
  if ( ! chipsComponentXSisInitialized ) {
    chipsComponentXSisInitialized = true;
    chipsComponentXSinstance = new G4ChipsComponentXS();
  }
  G4ChipsComponentXS* _instance = chipsComponentXSinstance;
  FTFxsManager = _instance;

  Plab /= GeV;
  G4double Xftf = 0.0;

  G4int NumberOfTargetProtons  = theZ; 
  G4int NumberOfTargetNeutrons = theA - theZ;
  NumberOfTargetNucleons = NumberOfTargetProtons + NumberOfTargetNeutrons;

  if ( ProjectilePDGcode == 2212  ||  ProjectilePDGcode == 2112 ) {  // Projectile is nucleon        
    G4ParticleDefinition* Proton = G4Proton::Proton();                                          //ALB 
    G4double XtotPP = FTFxsManager->GetTotalElementCrossSection( Proton, KineticEnergy, 1, 0 ); //ALB

    G4ParticleDefinition* Neutron = G4Neutron::Neutron();
    G4double XtotPN = FTFxsManager->GetTotalElementCrossSection( Neutron, KineticEnergy, 1, 0 ); //ALB
    G4double XelPP  = FTFxsManager->GetElasticElementCrossSection( Proton, KineticEnergy, 1, 0 );
    G4double XelPN  = FTFxsManager->GetElasticElementCrossSection( Neutron, KineticEnergy, 1, 0 );

    #ifdef debugFTFparams
    G4cout << "XsPP " << XtotPP/millibarn << " " << XelPP/millibarn << G4endl
           << "XsPN " << XtotPN/millibarn << " " << XelPN/millibarn << G4endl;
    #endif

    if ( ! ProjectileIsNucleus ) {  // Projectile is hadron
      Xtotal   = ( NumberOfTargetProtons * XtotPP + NumberOfTargetNeutrons * XtotPN ) /
                 NumberOfTargetNucleons;
      Xelastic = ( NumberOfTargetProtons * XelPP  + NumberOfTargetNeutrons * XelPN  ) / 
                 NumberOfTargetNucleons;
    } else {  // Projectile is a nucleus
      Xtotal = ( 
                  AbsProjectileCharge * NumberOfTargetProtons * XtotPP + 
                  ( AbsProjectileBaryonNumber - AbsProjectileCharge ) *
                      NumberOfTargetNeutrons * XtotPP 
                + 
                  ( AbsProjectileCharge * NumberOfTargetNeutrons +
                    ( AbsProjectileBaryonNumber - AbsProjectileCharge ) *
                        NumberOfTargetProtons ) * XtotPN
                ) / ( AbsProjectileBaryonNumber * NumberOfTargetNucleons );
      Xelastic= (
                  AbsProjectileCharge * NumberOfTargetProtons * XelPP + 
                  ( AbsProjectileBaryonNumber - AbsProjectileCharge ) *
                      NumberOfTargetNeutrons * XelPP 
                 + 
                  ( AbsProjectileCharge * NumberOfTargetNeutrons +
                    ( AbsProjectileBaryonNumber - AbsProjectileCharge ) *
                        NumberOfTargetProtons ) * XelPN
                ) / ( AbsProjectileBaryonNumber * NumberOfTargetNucleons );
    }

    Xannihilation = 0.0;
    Xtotal /= millibarn;
    Xelastic /= millibarn;

  } else if ( ProjectilePDGcode < -1000 ) { // Projectile is anti_baryon

    G4double X_a( 0.0 ), X_b( 0.0 ), X_c( 0.0 ), X_d( 0.0 );
    G4double MesonProdThreshold = ProjectileMass + TargetMass + 
                                  ( 2.0 * 0.14 + 0.016 ); // 2 Mpi + DeltaE;

    if ( PlabPerParticle < 40.0*MeV ) { // Low energy limits. Projectile at rest.
      Xtotal =   1512.9;    // mb
      Xelastic =  473.2;    // mb
      X_a =       625.1;    // mb
      X_b =         9.780;  // mb
      X_c =        49.989;  // mb
      X_d =         6.614;  // mb
    } else { // Total and elastic cross section of PbarP interactions a'la Arkhipov
      G4double LogS = G4Log( ECMSsqr / 33.0625 );
      G4double Xasmpt = 36.04 + 0.304*LogS*LogS;              // mb
      LogS = G4Log( SqrtS / 20.74 );
      G4double Basmpt = 11.92 + 0.3036*LogS*LogS;             // GeV^(-2)
      G4double R0 = std::sqrt( 0.40874044*Xasmpt - Basmpt );  // GeV^(-1)

      G4double FlowF = SqrtS / std::sqrt( ECMSsqr*ECMSsqr + ProjectileMass2*ProjectileMass2 +
                                          TargetMass2*TargetMass2 - 2.0*ECMSsqr*ProjectileMass2
                                          - 2.0*ECMSsqr*TargetMass2 
                                          - 2.0*ProjectileMass2*TargetMass2 );

      Xtotal = Xasmpt * ( 1.0 + 13.55*FlowF/R0/R0/R0*
                                (1.0 - 4.47/SqrtS + 12.38/ECMSsqr - 12.43/SqrtS/ECMSsqr) );    // mb

      Xasmpt = 4.4 + 0.101*LogS*LogS;  // mb
      Xelastic = Xasmpt * ( 1.0 + 59.27*FlowF/R0/R0/R0*
                                  (1.0 - 6.95/SqrtS + 23.54/ECMSsqr - 25.34/SqrtS/ECMSsqr ) ); // mb

      //G4cout << "Param Xtotal Xelastic " << Xtotal << " " << Xelastic << G4endl
      //       << "FlowF " << FlowF << " SqrtS " << SqrtS << G4endl
      //       << "Param Xelastic-NaN " << Xelastic << " " 
      //       << 1.5*16.654/pow(ECMSsqr/2.176/2.176,2.2) << " " << ECMSsqr << G4endl;

      X_a = 25.0*FlowF;  // mb, 3-shirts diagram

      if ( SqrtS < MesonProdThreshold ) {
        X_b = 3.13 + 140.0*G4Pow::GetInstance()->powA( MesonProdThreshold - SqrtS, 2.5 );  // mb anti-quark-quark annihilation
        Xelastic -= 3.0*X_b;  // Xel-X(PbarP->NNbar)
      } else {
        X_b = 6.8/SqrtS;  // mb anti-quark-quark annihilation
        Xelastic -= 3.0*X_b;  // Xel-X(PbarP->NNbar)
      }

      X_c = 2.0*FlowF*sqr( ProjectileMass + TargetMass )/ECMSsqr;  // mb rearrangement

      X_d = 23.3/ECMSsqr;  // mb anti-quark-quark string creation
    }

    //G4cout << "Param Xtotal Xelastic " << Xtotal << " " << Xelastic << G4endl
    //       << "Para a b c d " << X_a << " " << X_b << " " << X_c << " " << X_d << G4endl;
    //       << "Para a b c d " << X_a << " " << 5.*X_b << " " << 5.*X_c << " " << 6.*X_d 
    //       << G4endl;

    G4double Xann_on_P( 0.0), Xann_on_N( 0.0 );

    if ( ProjectilePDGcode == -2212 ) {              // Pbar+P/N
      Xann_on_P = X_a + X_b*5.0 + X_c*5.0 + X_d*6.0; 
      Xann_on_N = X_a + X_b*4.0 + X_c*4.0 + X_d*4.0;
    } else if ( ProjectilePDGcode == -2112 ) {       // NeutrBar+P/N
      Xann_on_P = X_a + X_b*4.0 + X_c*4.0 + X_d*4.0;
      Xann_on_N = X_a + X_b*5.0 + X_c*5.0 + X_d*6.0;
    } else if ( ProjectilePDGcode == -3122 ) {       // LambdaBar+P/N
      Xann_on_P = X_a + X_b*3.0 + X_c*3.0 + X_d*2.0;
      Xann_on_N = X_a + X_b*3.0 + X_c*3.0 + X_d*2.0;
    } else if ( ProjectilePDGcode == -3112 ) {       // Sigma-Bar+P/N
      Xann_on_P = X_a + X_b*2.0 + X_c*2.0 + X_d*0.0;
      Xann_on_N = X_a + X_b*4.0 + X_c*4.0 + X_d*2.0;
    } else if ( ProjectilePDGcode == -3212 ) {       // Sigma0Bar+P/N
      Xann_on_P = X_a + X_b*3.0 + X_c*3.0 + X_d*2.0;
      Xann_on_N = X_a + X_b*3.0 + X_c*3.0 + X_d*2.0;
    } else if ( ProjectilePDGcode == -3222 ) {       // Sigma+Bar+P/N
      Xann_on_P = X_a + X_b*4.0 + X_c*4.0 + X_d*2.0;
      Xann_on_N = X_a + X_b*2.0 + X_c*2.0 + X_d*0.0;
    } else if ( ProjectilePDGcode == -3312 ) {       // Xi-Bar+P/N
      Xann_on_P = X_a + X_b*1.0 + X_c*1.0 + X_d*0.0;
      Xann_on_N = X_a + X_b*2.0 + X_c*2.0 + X_d*0.0;
    } else if ( ProjectilePDGcode == -3322 ) {       // Xi0Bar+P/N
      Xann_on_P = X_a + X_b*2.0 + X_c*2.0 + X_d*0.0;
      Xann_on_N = X_a + X_b*1.0 + X_c*1.0 + X_d*0.0;
    } else if ( ProjectilePDGcode == -3334 ) {       // Omega-Bar+P/N
      Xann_on_P = X_a + X_b*0.0 + X_c*0.0 + X_d*0.0;
      Xann_on_N = X_a + X_b*0.0 + X_c*0.0 + X_d*0.0;
    } else {
      G4cout << "Unknown anti-baryon for FTF annihilation" << G4endl;
    }

    //G4cout << "Sum          " << Xann_on_P << G4endl;

    if ( ! ProjectileIsNucleus ) {  // Projectile is anti-baryon
      Xannihilation = ( NumberOfTargetProtons * Xann_on_P  + NumberOfTargetNeutrons * Xann_on_N  )
                      / NumberOfTargetNucleons;
    } else {  // Projectile is a nucleus
      Xannihilation = (
                        ( AbsProjectileCharge * NumberOfTargetProtons + 
                          ( AbsProjectileBaryonNumber - AbsProjectileCharge ) *
                          NumberOfTargetNeutrons ) * Xann_on_P 
                       + 
                        ( AbsProjectileCharge * NumberOfTargetNeutrons +
                          ( AbsProjectileBaryonNumber - AbsProjectileCharge ) *
                          NumberOfTargetProtons ) * Xann_on_N
                      ) / ( AbsProjectileBaryonNumber * NumberOfTargetNucleons );
    }

    //G4double Xftf = 0.0;  
    MesonProdThreshold = ProjectileMass + TargetMass + (0.14 + 0.08); // Mpi + DeltaE
    if ( SqrtS > MesonProdThreshold ) {
      Xftf = 36.0 * ( 1.0 - MesonProdThreshold/SqrtS );
    }

    Xtotal = Xelastic + Xannihilation + Xftf;

    #ifdef debugFTFparams
    G4cout << "Plab Xtotal, Xelastic  Xinel Xftf " << Plab << " " << Xtotal << " " << Xelastic
           << " " << Xtotal - Xelastic << " " << Xtotal - Xelastic - Xannihilation << G4endl
           << "Plab Xelastic/Xtotal,  Xann/Xin " << Plab << " " << Xelastic/Xtotal << " " 
           << Xannihilation/(Xtotal - Xelastic) << G4endl;
    #endif

  } else if ( ProjectilePDGcode == 211 ) {  // Projectile is PionPlus

    G4double XtotPiP = FTFxsManager->GetTotalElementCrossSection( particle, KineticEnergy, 1, 0 ); 
    G4ParticleDefinition* PionMinus = G4PionMinus::PionMinus();
    G4double XtotPiN = FTFxsManager->GetTotalElementCrossSection( PionMinus, KineticEnergy, 1, 0 ); 
    G4double XelPiP  = FTFxsManager->GetElasticElementCrossSection( particle, KineticEnergy, 1, 0 ); 
    G4double XelPiN  = FTFxsManager->GetElasticElementCrossSection( PionMinus, KineticEnergy, 1, 0 );
    Xtotal   = ( NumberOfTargetProtons * XtotPiP + NumberOfTargetNeutrons * XtotPiN ) 
               / NumberOfTargetNucleons;
    Xelastic = ( NumberOfTargetProtons * XelPiP  + NumberOfTargetNeutrons * XelPiN ) 
               / NumberOfTargetNucleons; 
    Xannihilation = 0.0;
    Xtotal /= millibarn;
    Xelastic /= millibarn;
  
  } else if ( ProjectilePDGcode == -211 ) {  // Projectile is PionMinus
 
    G4double XtotPiP = FTFxsManager->GetTotalElementCrossSection( particle, KineticEnergy, 1, 0 );
    G4ParticleDefinition* PionPlus = G4PionPlus::PionPlus();
    G4double XtotPiN = FTFxsManager->GetTotalElementCrossSection( PionPlus, KineticEnergy, 1, 0 );   
    G4double XelPiP  = FTFxsManager->GetElasticElementCrossSection( particle, KineticEnergy, 1, 0 );
    G4double XelPiN  = FTFxsManager->GetElasticElementCrossSection( PionPlus, KineticEnergy, 1, 0 );
    Xtotal   = ( NumberOfTargetProtons * XtotPiP + NumberOfTargetNeutrons * XtotPiN )
               / NumberOfTargetNucleons;
    Xelastic = ( NumberOfTargetProtons * XelPiP  + NumberOfTargetNeutrons * XelPiN )
               / NumberOfTargetNucleons;
    Xannihilation = 0.0;
    Xtotal /= millibarn;
    Xelastic /= millibarn;
      
  } else if ( ProjectilePDGcode == 111 ) {  // Projectile is PionZero
      
    G4ParticleDefinition* PionPlus = G4PionPlus::PionPlus();
    G4double XtotPipP = FTFxsManager->GetTotalElementCrossSection( PionPlus, KineticEnergy, 1, 0 );
    G4ParticleDefinition* PionMinus = G4PionMinus::PionMinus();
    G4double XtotPimP = FTFxsManager->GetTotalElementCrossSection( PionMinus, KineticEnergy, 1, 0 ); 
    G4double XelPipP = FTFxsManager->GetElasticElementCrossSection( PionPlus, KineticEnergy, 1, 0 );
    G4double XelPimP = FTFxsManager->GetElasticElementCrossSection( PionMinus, KineticEnergy, 1, 0 );
    G4double XtotPiP = ( XtotPipP + XtotPimP ) / 2.0;
    G4double XtotPiN = XtotPiP;
    G4double XelPiP = ( XelPipP  + XelPimP ) / 2.0;
    G4double XelPiN = XelPiP;
    Xtotal   = ( NumberOfTargetProtons * XtotPiP + NumberOfTargetNeutrons * XtotPiN )
               / NumberOfTargetNucleons;
    Xelastic = ( NumberOfTargetProtons * XelPiP  + NumberOfTargetNeutrons * XelPiN )
               / NumberOfTargetNucleons; 
    Xannihilation = 0.0;
    Xtotal /= millibarn;
    Xelastic /= millibarn;
      
  } else if ( ProjectilePDGcode == 321 ) {  // Projectile is KaonPlus

    G4double XtotKP = FTFxsManager->GetTotalElementCrossSection( particle, KineticEnergy, 1, 0 );
    G4ParticleDefinition* KaonMinus = G4KaonMinus::KaonMinus();
    G4double XtotKN = FTFxsManager->GetTotalElementCrossSection( KaonMinus, KineticEnergy, 1, 0 );
    G4double XelKP  = FTFxsManager->GetElasticElementCrossSection( particle, KineticEnergy, 1, 0 );
    G4double XelKN  = FTFxsManager->GetElasticElementCrossSection( KaonMinus, KineticEnergy, 1, 0 );
    Xtotal   = ( NumberOfTargetProtons * XtotKP + NumberOfTargetNeutrons * XtotKN )
               / NumberOfTargetNucleons;
    Xelastic = ( NumberOfTargetProtons * XelKP  + NumberOfTargetNeutrons * XelKN )
               / NumberOfTargetNucleons;
    Xannihilation = 0.0;
    Xtotal /= millibarn;
    Xelastic /= millibarn;

  } else if ( ProjectilePDGcode == -321 ) {  // Projectile is KaonMinus

    G4double XtotKP = FTFxsManager->GetTotalElementCrossSection( particle, KineticEnergy, 1, 0 );
    G4ParticleDefinition* KaonPlus = G4KaonPlus::KaonPlus();
    G4double XtotKN = FTFxsManager->GetTotalElementCrossSection( KaonPlus, KineticEnergy, 1, 0 );
    G4double XelKP  = FTFxsManager->GetElasticElementCrossSection( particle, KineticEnergy, 1, 0 );
    G4double XelKN  = FTFxsManager->GetElasticElementCrossSection( KaonPlus, KineticEnergy, 1, 0 ); 
    Xtotal   = ( NumberOfTargetProtons * XtotKP + NumberOfTargetNeutrons * XtotKN )
               / NumberOfTargetNucleons;
    Xelastic = ( NumberOfTargetProtons * XelKP  + NumberOfTargetNeutrons * XelKN )
               / NumberOfTargetNucleons;
    Xannihilation = 0.0;
    Xtotal /= millibarn;
    Xelastic /= millibarn;

  } else if ( ProjectilePDGcode == 311  ||  ProjectilePDGcode == 130  ||  
              ProjectilePDGcode == 310 ) {  // Projectile is KaonZero

    G4ParticleDefinition* KaonPlus = G4KaonPlus::KaonPlus();
    G4double XtotKpP = FTFxsManager->GetTotalElementCrossSection( KaonPlus, KineticEnergy, 1, 0 );
    G4ParticleDefinition* KaonMinus = G4KaonMinus::KaonMinus();
    G4double XtotKmP = FTFxsManager->GetTotalElementCrossSection( KaonMinus, KineticEnergy, 1, 0 );
    G4double XelKpP = FTFxsManager->GetElasticElementCrossSection( KaonPlus, KineticEnergy, 1, 0 );
    G4double XelKmP = FTFxsManager->GetElasticElementCrossSection( KaonMinus, KineticEnergy, 1, 0 );
    G4double XtotKP = ( XtotKpP + XtotKmP ) / 2.0;
    G4double XtotKN = XtotKP;
    G4double XelKP = ( XelKpP + XelKmP ) / 2.0; 
    G4double XelKN = XelKP;
    Xtotal   = ( NumberOfTargetProtons * XtotKP + NumberOfTargetNeutrons * XtotKN )
               / NumberOfTargetNucleons;
    Xelastic = ( NumberOfTargetProtons * XelKP  + NumberOfTargetNeutrons * XelKN )
               / NumberOfTargetNucleons;
    Xannihilation = 0.0;
    Xtotal /= millibarn;
    Xelastic /= millibarn;

  } else {  // Projectile is undefined, Nucleon assumed

    G4ParticleDefinition* Proton = G4Proton::Proton();
    G4double XtotPP = FTFxsManager->GetTotalElementCrossSection( Proton, KineticEnergy, 1, 0 );
    G4ParticleDefinition* Neutron = G4Neutron::Neutron();
    G4double XtotPN = FTFxsManager->GetTotalElementCrossSection( Neutron, KineticEnergy, 1, 0 );
    G4double XelPP  = FTFxsManager->GetElasticElementCrossSection( Proton, KineticEnergy, 1, 0 );
    G4double XelPN  = FTFxsManager->GetElasticElementCrossSection( Neutron, KineticEnergy, 1, 0 );
    Xtotal   = ( NumberOfTargetProtons  * XtotPP + NumberOfTargetNeutrons * XtotPN )
               / NumberOfTargetNucleons;
    Xelastic = ( NumberOfTargetProtons  * XelPP  + NumberOfTargetNeutrons * XelPN )
               / NumberOfTargetNucleons;
    Xannihilation = 0.0;
    Xtotal /= millibarn;
    Xelastic /= millibarn;

  };

  // Geometrical parameters
  SetTotalCrossSection( Xtotal );
  SetElastisCrossSection( Xelastic );
  SetInelasticCrossSection( Xtotal - Xelastic );

  //G4cout << "Plab Xtotal, Xelastic Xinel Xftf " << Plab << " " << Xtotal << " " << Xelastic 
  //       << " " << Xtotal - Xelastic << " " << Xtotal - Xelastic - Xannihilation << G4endl;
  //if (Xtotal - Xelastic != 0.0 ) {
  //  G4cout << "Plab Xelastic/Xtotal,  Xann/Xin " << Plab << " " << Xelastic/Xtotal 
  //         << " " << Xannihilation / (Xtotal - Xelastic) << G4endl;
  //} else {
  //  G4cout << "Plab Xelastic/Xtotal,  Xann     " << Plab << " " << Xelastic/Xtotal
  //         << " " << Xannihilation << G4endl;
  //}
  //G4int Uzhi; G4cin >> Uzhi;

  // Interactions with elastic and inelastic collisions
  SetProbabilityOfElasticScatt( Xtotal, Xelastic );
  SetRadiusOfHNinteractions2( Xtotal/pi/10.0 );
  if ( Xtotal - Xelastic == 0.0 ) {
    SetProbabilityOfAnnihilation( 0.0 );
  } else {
    SetProbabilityOfAnnihilation( Xannihilation / (Xtotal - Xelastic) );
  }

  // No elastic scattering 
  //SetProbabilityOfElasticScatt( Xtotal, 0.0 );
  //SetRadiusOfHNinteractions2( (Xtotal - Xelastic)/pi/10.0 );
  //SetProbabilityOfAnnihilation( 1.0 );
  //SetProbabilityOfAnnihilation( 0.0 );

  SetSlope( Xtotal*Xtotal/16.0/pi/Xelastic/0.3894 ); // Slope parameter of elastic scattering
                                                     // (GeV/c)^(-2))
  //G4cout << "Slope " << GetSlope() << G4endl;
  SetGamma0( GetSlope()*Xtotal/10.0/2.0/pi );

  // Parameters of elastic scattering
  // Gaussian parametrization of elastic scattering amplitude assumed
  SetAvaragePt2ofElasticScattering( 1.0/( Xtotal*Xtotal/16.0/pi/Xelastic/0.3894 )*GeV*GeV );
//  G4cout << "AvaragePt2ofElasticScattering " << GetAvaragePt2ofElasticScattering() << G4endl;

  // Parameters of excitations

  G4double Xinel = Xtotal - Xelastic;  // Uzhi 25.04.2012
  //G4cout << "Param ProjectilePDGcode " << ProjectilePDGcode << G4endl;

  if ( ProjectilePDGcode > 1000 ) {  // Projectile is baryon
    //        Proc#   A1      B1            A2       B2   A3   Atop       Ymin
//    SetParams( 0,     13.71, 1.75,          -214.5, 4.25, 0.0, 0.5  ,     1.1 );  // Qexchange without Exc.
    SetParams( 0,     13.71, 1.75,          -30.69, 3.0 , 0.0, 1.0  ,     0.93 );   // Qexchange without Exc.
    SetParams( 1,     25.0 , 1.0 ,          -50.34, 1.5 , 0.0, 0.0  ,     1.4 );    // Qexchange with    Exc.
    if( Xinel > 0.) {
      SetParams( 2, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,     0.93);// Projectile diffraction
      SetParams( 3, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,     0.93);// Target diffraction
      SetParams( 4,      0.6 , 0.0 ,           -1.20, 0.5 , 0.0, 0.0  ,     1.4 );// Qexchange with Exc. Additional multiply
    } else {
      SetParams( 2, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0  ,     0.0);
      SetParams( 3, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0  ,     0.0);
      SetParams( 4, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0  ,     0.0);
    }

    if ( AbsProjectileBaryonNumber > 1  ||  NumberOfTargetNucleons > 1 ) {
// It is not decided what to do with diffraction dissociation in Had-Nucl and Nucl-Nucl interactions
      SetParams( 2,       0.0, 0.0 ,           0.0  , 0.0 , 0.0, 0.0   , -100.0  );  // Projectile diffraction
//      SetParams( 3,       0.0, 0.0 ,           0.0  , 0.0 , 0.0, 0.0   , -100.0  );  // Target diffraction
    }

    SetDeltaProbAtQuarkExchange( 0.0 );
    if ( NumberOfTargetNucleons > 26 ) {
      SetProbOfSameQuarkExchange( 1.0);
    } else {
      SetProbOfSameQuarkExchange( 0.0 );
    }
    SetProjMinDiffMass( 1.16 );     // GeV 
    SetProjMinNonDiffMass( 1.16 );  // GeV 
    SetTarMinDiffMass( 1.16 );      // GeV
    SetTarMinNonDiffMass( 1.16 );   // GeV 
    SetAveragePt2( 0.3 );           // GeV^2                    // Uzhi Oct  2014
    SetProbLogDistrPrD( 0.6 );                                  // Uzhi June 2016, 0.5 0.3
    SetProbLogDistr(0.6);                                       // Uzhi June 2016  0.4->0.6

  } else if( ProjectilePDGcode < -1000 ) {  // Projectile is anti_baryon

    //        Proc#   A1      B1            A2       B2   A3   Atop       Ymin
    SetParams( 0,      0.0 , 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,  1000.0  );  // Qexchange without Exc. 
    SetParams( 1,      0.0 , 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,  1000.0  );  // Qexchange with    Exc.
    if( Xinel > 0.) {
      SetParams( 2, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,     0.93);  // Projectile diffraction
      SetParams( 3, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,     0.93);  // Target diffraction
      SetParams( 4,       1.0, 0.0 ,             0.0, 0.0 , 0.0, 0.0  ,    0.93 );  // Qexchange with    Exc. Additional multiply
    } else {
      SetParams( 2, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0  ,     0.0);
      SetParams( 3, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0  ,     0.0);
      SetParams( 4, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0  ,     0.0);
    }

    if ( AbsProjectileBaryonNumber > 1  ||  NumberOfTargetNucleons > 1 ) {
      SetParams( 2,      0.0 , 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,  -100.0  );  // Projectile diffraction
//      SetParams( 3,      0.0 , 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,  -100.0  );  // Target diffraction

    }
    SetDeltaProbAtQuarkExchange( 0.0 );
    SetProbOfSameQuarkExchange( 0.0 );
    SetProjMinDiffMass( ProjectileMass + 0.22 );     // GeV 
    SetProjMinNonDiffMass( ProjectileMass + 0.22 );  // GeV
    SetTarMinDiffMass( TargetMass + 0.22 );          // GeV
    SetTarMinNonDiffMass( TargetMass + 0.22 );       // GeV
    SetAveragePt2( 0.3 );                           // GeV^2    // Uzhi Oct 2014
    SetProbLogDistrPrD( 0.6 );                                  // Uzhi June 2016  0.4->0.6 
    SetProbLogDistr( 0.6 );					// Uzhi June 2016  0.4->0.6
            
  } else if ( ProjectileabsPDGcode == 211  ||  ProjectilePDGcode ==  111 ) {  // Projectile is Pion 

    //        Proc#   A1      B1            A2       B2   A3   Atop       Ymin
    SetParams( 0,  720.0,    2.5 ,         2.3 ,     1.0,    0.,   1. ,       2.7 ); 
    SetParams( 1,  12.87,    0.5 ,       -44.91,     1.0,    0.,   0. ,       2.5 );
    SetParams( 2,  0.086,    0.  ,        -0.3 ,     0.5,    0.,   0. ,       2.5 ); 
    SetParams( 3,   32.8,    1.0 ,      -114.5 ,     1.5, 0.084,   0. ,       2.5 );
    SetParams( 4,    1.0,    0.0 ,        -3.49,     0.5,   0.0,   0. ,       2.5 );  // Qexchange with    Exc. Additional multiply

    if ( AbsProjectileBaryonNumber > 1  ||  NumberOfTargetNucleons > 1 ) {
      SetParams( 2,      0.0 , 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,  -100.0  );  // Projectile diffraction
//      SetParams( 3,      0.0 , 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,  -100.0  );  // Target diffraction
    }

    SetDeltaProbAtQuarkExchange( 0.56 );  // (0.35)
    SetProjMinDiffMass( 0.5 );            // (0.5)  // GeV
    SetProjMinNonDiffMass( 0.5 );         // (0.5)  // GeV 
    SetTarMinDiffMass( 1.16 );                      // GeV
    SetTarMinNonDiffMass( 1.16 );                   // GeV
    SetAveragePt2( 0.3 );                          // GeV^2     // Uzhi Oct 2014
    SetProbLogDistrPrD( 0.6 );                                  // Uzhi June 2016  0.4->0.6
    SetProbLogDistr( 0.6 );					// Uzhi June 2016  0.4->0.6

  } else if ( ProjectileabsPDGcode == 321  ||  ProjectileabsPDGcode == 311  || 
              ProjectilePDGcode == 130     ||  ProjectilePDGcode == 310 ) {  // Projectile is Kaon

    //        Proc#   A1      B1            A2       B2   A3   Atop       Ymin
    SetParams( 0,     60.0 , 2.5 ,           0.0  , 0.0 , 0.0, 0.0  ,  -100.0  );  // Qexchange without Exc. 
    SetParams( 1,      6.0 , 1.0 ,         -24.33 , 2.0 , 0.0, 0.0  ,     1.40 );  // Qexchange with    Exc.
    SetParams( 2,      2.76, 1.2 ,         -22.5  , 2.7 ,0.04, 0.0  ,     1.40 );  // Projectile diffraction
    SetParams( 3,      1.09, 0.5 ,          -8.88 , 2.  ,0.05, 0.0  ,     1.40 );  // Target diffraction
    SetParams( 4,       1.0, 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,     0.93 );  // Qexchange with    Exc. Additional multiply

    if ( AbsProjectileBaryonNumber > 1  ||  NumberOfTargetNucleons > 1 ) {
      SetParams( 2,      0.0 , 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,  -100.0  );  // Projectile diffraction
//      SetParams( 3,      0.0 , 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,  -100.0  );  // Target diffraction
    }
    SetDeltaProbAtQuarkExchange( 0.6 );
    SetProjMinDiffMass( 0.7 );     // (1.4) // (0.7) // GeV 
    SetProjMinNonDiffMass( 0.7 );  // (1.4) // (0.7) // GeV 
    SetTarMinDiffMass( 1.16 );                       // GeV
    SetTarMinNonDiffMass( 1.16 );                    // GeV
    SetAveragePt2( 0.3 );                           // GeV^2    // Uzhi Oct  2014
    SetProbLogDistrPrD( 0.6 );                                  // Uzhi June 2016  0.4->0.6
    SetProbLogDistr( 0.6 );                                     // Uzhi June 2016  0.4->0.6

   } else {  // Projectile is undefined, Nucleon assumed

    //        Proc#   A1      B1            A2       B2   A3   Atop       Ymin
//    SetParams( 0,     13.71, 1.75,          -214.5, 4.25, 0.0, 0.5  ,     1.1 );// Qexchange without Exc. May 2016
    SetParams( 0,     13.71, 1.75,          -30.69, 3.0 , 0.0, 1.0  ,     0.93 ); // Qexchange without Exc.
    SetParams( 1,      25.0, 1.0,          -50.34,  1.5 , 0.0, 0.0  ,     1.4 );  // Qexchange with    Exc.
    if( Xinel > 0.) {
      SetParams( 2, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,   0.93);  // Projectile diffraction
      SetParams( 3, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,   0.93);  // Target diffraction
      SetParams( 4,       1.0, 0.0 ,          -2.01 , 0.5 , 0.0, 0.0  ,   1.4 );  // Qexchange with    Exc. Additional multiply
    } else {
      SetParams( 2, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0  ,     0.0);
      SetParams( 3, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0  ,     0.0);
      SetParams( 4, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0  ,     0.0);
    }
    if ( AbsProjectileBaryonNumber > 1  ||  NumberOfTargetNucleons > 1 ) {
      SetParams( 2,      0.0 , 0.0 ,            0.0 , 0.0 , 0.0, 0.0  ,  -100.0  );  // Projectile diffraction
//      SetParams( 3,      0.0 , 0.0 ,            0.0 , 0.0 , 0.0, 0.0  ,  -100.0  );  // Target diffraction
    }
    SetDeltaProbAtQuarkExchange( 0.0 );              // 7 June 2011
    SetProbOfSameQuarkExchange( 0.0 );
    SetProjMinDiffMass( ProjectileMass + 0.22 );     // GeV 
    SetProjMinNonDiffMass( ProjectileMass + 0.22 );  // GeV
    SetTarMinDiffMass( TargetMass + 0.22 );          // GeV
    SetTarMinNonDiffMass( TargetMass + 0.22 );       // GeV
    SetAveragePt2( 0.3 );                           // GeV^2    // Uzhi Oct 2014
    SetProbLogDistrPrD( 0.6 );                                  // Uzhi June 2016  0.4->0.6
    SetProbLogDistr( 0.6 );					// Uzhi June 2016  0.4->0.6

  }

 // Set parameters of a string kink
 //  SetPt2Kink( 6.0*GeV*GeV );                                   // Uzhi Oct 2014
  SetPt2Kink( 0.0*GeV*GeV );                                      // Uzhi Oct 2014   // Uzhi to switch off kinky strings
  G4double Puubar( 1.0/3.0 ), Pddbar( 1.0/3.0 ), Pssbar( 1.0/3.0 );  // SU(3) symmetry
//  G4double Puubar( 0.41 ), Pddbar( 0.41 ), Pssbar( 0.18 );  // Broken SU(3) symmetry
  SetQuarkProbabilitiesAtGluonSplitUp( Puubar, Pddbar, Pssbar );

 // Set parameters of nuclear destruction
 if ( ProjectileabsPDGcode < 1000 ) {  // Meson projectile
   SetMaxNumberOfCollisions( Plab, 2.0 );  //  3.0 )
   SetCofNuclearDestruction( 0.00481*G4double(NumberOfTargetNucleons)*             // Uzhi 3.05.2015
           G4Exp( 4.0*(Ylab - 2.1) )/( 1.0 + G4Exp( 4.0*(Ylab - 2.1) ) ) );
   SetR2ofNuclearDestruction( 1.5*fermi*fermi );
   SetDofNuclearDestruction( 0.3 );
   SetPt2ofNuclearDestruction( ( 0.035 + 0.04*G4Exp( 4.0*(Ylab - 2.5) )/
                                         ( 1.0 + G4Exp( 4.0*(Ylab - 2.5) ) ) )*GeV*GeV );
   SetMaxPt2ofNuclearDestruction( 1.0*GeV*GeV );
   SetExcitationEnergyPerWoundedNucleon( 40.0*MeV );                               // Uzhi March 2015: 100 -> 40
 } else if ( ProjectilePDGcode < -1000 ) {  // for anti-baryon projectile
   SetMaxNumberOfCollisions( Plab, 2.0 );  // 3.0 )
   SetCofNuclearDestruction( 0.00481*G4double(NumberOfTargetNucleons)*             // Uzhi 3.05.2015
           G4Exp( 4.0*(Ylab - 2.1) )/( 1.0 + G4Exp( 4.0*(Ylab - 2.1) ) ) );
   SetR2ofNuclearDestruction( 1.5*fermi*fermi );
   SetDofNuclearDestruction( 0.3 );
   SetPt2ofNuclearDestruction( ( 0.035 + 0.04*G4Exp( 4.0*(Ylab - 2.5) )/
                                         ( 1.0 + G4Exp( 4.0*(Ylab - 2.5) ) ) )*GeV*GeV );
   SetMaxPt2ofNuclearDestruction( 1.0*GeV*GeV );
   SetExcitationEnergyPerWoundedNucleon( 40.0*MeV );                              // Uzhi March 2015: 100 -> 20
   if ( Plab < 2.0 ) {  // 2 GeV/c
     // For slow anti-baryon we have to garanty putting on mass-shell
     SetCofNuclearDestruction( 0.0 );
     SetR2ofNuclearDestruction( 1.5*fermi*fermi );
     SetDofNuclearDestruction( 0.01 );
     SetPt2ofNuclearDestruction( 0.035*GeV*GeV );
     SetMaxPt2ofNuclearDestruction( 0.04*GeV*GeV );
   }
 } else {  // Projectile baryon assumed
   SetMaxNumberOfCollisions( Plab, 2.0 ); // 3.0 )

   SetCofNuclearDestructionPr( 0.00481*G4double(AbsProjectileBaryonNumber)*
           G4Exp( 4.0*(Ylab - 2.1) )/( 1.0 + G4Exp( 4.0*(Ylab - 2.1) ) ) );

   SetCofNuclearDestruction(   0.00481*G4double(NumberOfTargetNucleons)*
           G4Exp( 4.0*(Ylab - 2.1) )/( 1.0 + G4Exp( 4.0*(Ylab - 2.1) ) ) );
   SetR2ofNuclearDestruction( 1.5*fermi*fermi );
   SetDofNuclearDestruction( 0.3 );
   SetPt2ofNuclearDestruction( ( 0.035 + 0.04*G4Exp( 4.0*(Ylab - 2.5) )/
                                         ( 1.0 + G4Exp( 4.0*(Ylab - 2.5) ) ) )*GeV*GeV);  
//G4cout<<"Pt2 "<<std::sqrt(GetPt2ofNuclearDestruction())<<" "<<std::sqrt(GetPt2ofNuclearDestruction()/2.)<<G4endl;
//{G4int Uzhi; G4cin>>Uzhi;}
   SetMaxPt2ofNuclearDestruction( 9.0*GeV*GeV );
   SetExcitationEnergyPerWoundedNucleon( 40.0*MeV );                              // Uzhi March 2015: 100 -> 40
 }

 //SetCofNuclearDestruction( 0.47*G4Exp( 2.0*(Ylab - 2.5) )/( 1.0 + G4Exp( 2.0*(Ylab - 2.5) ) ) ); 
 //SetPt2ofNuclearDestruction( ( 0.035 + 0.1*G4Exp( 4.0*(Ylab - 3.0) )/( 1.0 + G4Exp( 4.0*(Ylab - 3.0) ) ) )*GeV*GeV );

 //SetMagQuarkExchange( 120.0 );       // 210.0 PipP
 //SetSlopeQuarkExchange( 2.0 );
 //SetDeltaProbAtQuarkExchange( 0.6 );
 //SetProjMinDiffMass( 0.7 );          // GeV 1.1
 //SetProjMinNonDiffMass( 0.7 );       // GeV
 //SetProbabilityOfProjDiff( 0.0);     // 0.85*G4Pow::GetInstance()->powA( s/GeV/GeV, -0.5 ) ); // 40/32 X-dif/X-inel
 //SetTarMinDiffMass( 1.1 );           // GeV
 //SetTarMinNonDiffMass( 1.1 );        // GeV
 //SetProbabilityOfTarDiff( 0.0 );     // 0.85*G4Pow::GetInstance()->powA( s/GeV/GeV, -0.5 ) ); // 40/32 X-dif/X-inel

//SetAveragePt2( 0.0 );               // GeV^2   0.3
 //------------------------------------
//SetProbabilityOfElasticScatt( 1.0, 0.0); //1.0);                            //(Xtotal, Xelastic);
//SetProbabilityOfProjDiff( 1.0*0.62*G4Pow::GetInstance()->powA( s/GeV/GeV, -0.51 ) );  // 0->1
//SetProbabilityOfTarDiff( 4.0*0.62*G4Pow::GetInstance()->powA( s/GeV/GeV, -0.51 ) );   // 2->4
//SetAveragePt2( 0.3 );                                               // (0.15)
//SetAvaragePt2ofElasticScattering( 0.0 );

//SetMaxNumberOfCollisions( Plab, 6.0 ); //(4.0*(Plab + 0.01), Plab); // 6.0 );
//SetAveragePt2( 0.15 );
//SetCofNuclearDestruction(1.);//( 0.5 );                           // (0.25)             
//SetExcitationEnergyPerWoundedNucleon(0.);//( 30.0*MeV );                    // (75.0*MeV) 
//SetDofNuclearDestruction(0.6);//( 0.2 ); //0.4                                     // 0.3 0.5
//SetPt2ofNuclearDestruction(0.);//(2.*0.075*GeV*GeV); //( 0.3*GeV*GeV ); // (0.168*GeV*GeV) 

//SetMaxNumberOfCollisions( Plab, 78.0 ); // 3.0 )

//G4cout << "Cnd " << GetCofNuclearDestruction() << G4endl;
//G4cout << "Dnd " << GetDofNuclearDestruction() << G4endl;
//G4cout << "Pt2 " << GetPt2ofNuclearDestruction()/GeV/GeV << G4endl;
//G4int Uzhi; G4cin >> Uzhi;

} 


//============================================================================

G4double G4FTFParameters::GetProcProb( const G4int ProcN, const G4double y ) {
  G4double Prob( 0.0 );
  if ( y < ProcParams[ProcN][6] ) {
    Prob = ProcParams[ProcN][5]; 
    if(Prob < 0.) Prob=0.;
    return Prob;
  }
  Prob = ProcParams[ProcN][0] * G4Exp( -ProcParams[ProcN][1]*y ) +
         ProcParams[ProcN][2] * G4Exp( -ProcParams[ProcN][3]*y ) +
         ProcParams[ProcN][4];
  if(Prob < 0.) Prob=0.;
  return Prob;
}
