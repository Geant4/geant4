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

#include "G4CrossSectionDataSetRegistry.hh"
#include "G4VComponentCrossSection.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4LundStringFragmentation.hh"

#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

#include "G4HadronicDeveloperParameters.hh"

//============================================================================

//#define debugFTFparams

//============================================================================

G4FTFParameters::G4FTFParameters() 
{
  StringMass = new G4LundStringFragmentation;  // for estimation of min. mass of diffr. states
  Reset();
  csGGinstance = 
    G4CrossSectionDataSetRegistry::Instance()->GetComponentCrossSection("Glauber-Gribov");
  if (!csGGinstance) {
    csGGinstance = new G4ComponentGGHadronNucleusXsc();
  }

  // Set parameters of a string kink
  SetPt2Kink( 0.0*GeV*GeV );  // To switch off kinky strings (bad results obtained with 6.0*GeV*GeV)
  G4double Puubar( 1.0/3.0 ), Pddbar( 1.0/3.0 ), Pssbar( 1.0/3.0 );  // SU(3) symmetry
  //G4double Puubar( 0.41 ), Pddbar( 0.41 ), Pssbar( 0.18 );         // Broken SU(3) symmetry
  SetQuarkProbabilitiesAtGluonSplitUp( Puubar, Pddbar, Pssbar );
}

//============================================================================

void G4FTFParameters::InitForInteraction( const G4ParticleDefinition* particle, 
                                          G4int theA, G4int theZ, G4double PlabPerParticle ) 
{
  Reset();

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
    AbsProjectileCharge       = std::abs( G4int( particle->GetPDGCharge() ) );
    if ( ProjectileBaryonNumber > 1 ) {
      ProjectilePDGcode = 2212; ProjectileabsPDGcode = 2212;  // Proton
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

  G4double Ylab, Xtotal( 0.0 ), Xelastic( 0.0 ), Xannihilation( 0.0 );
  G4int NumberOfTargetNucleons;

  Ylab = 0.5 * G4Log( (Elab + Plab)/(Elab - Plab) );

  G4double ECMSsqr = S/GeV/GeV;
  G4double SqrtS   = std::sqrt( S )/GeV;

  #ifdef debugFTFparams
  G4cout << "Sqrt(s) " << SqrtS << G4endl;
  #endif

  TargetMass     /= GeV; TargetMass2     /= (GeV*GeV);
  ProjectileMass /= GeV; ProjectileMass2 /= (GeV*GeV);

  Plab /= GeV;
  G4double Xftf = 0.0;

  G4int NumberOfTargetProtons  = theZ; 
  G4int NumberOfTargetNeutrons = theA - theZ;
  NumberOfTargetNucleons = NumberOfTargetProtons + NumberOfTargetNeutrons;

  // ---------- hadron projectile ----------------
  if ( AbsProjectileBaryonNumber <= 1 ) {  // Projectile is hadron or baryon 

    // Interaction on P
    G4double xTtP = csGGinstance->GetTotalIsotopeCrossSection( particle, KineticEnergy, 1, 1);
    G4double xElP = csGGinstance->GetElasticIsotopeCrossSection(particle, KineticEnergy, 1, 1);

    // Interaction on N
    G4double xTtN = csGGinstance->GetTotalIsotopeCrossSection( particle, KineticEnergy, 0, 1);
    G4double xElN = csGGinstance->GetElasticIsotopeCrossSection(particle, KineticEnergy, 0, 1);

    // Average properties of h+N interactions
    Xtotal   = ( NumberOfTargetProtons * xTtP + NumberOfTargetNeutrons * xTtN ) / NumberOfTargetNucleons;
    Xelastic = ( NumberOfTargetProtons * xElP + NumberOfTargetNeutrons * xElN ) / NumberOfTargetNucleons; 
    Xannihilation = 0.0;

    Xtotal /= millibarn;
    Xelastic /= millibarn;

    #ifdef debugFTFparams
    G4cout<<"Estimated cross sections (total and elastic) of h+N interactions "<<Xtotal<<" "<<Xelastic<<" (mb)"<<G4endl;
    #endif
  }

  // ---------- nucleus projectile ----------------
  if ( ProjectileIsNucleus  &&  ProjectileBaryonNumber > 1 ) {

    #ifdef debugFTFparams
    G4cout<<"Projectile is a nucleus: A and Z - "<<ProjectileBaryonNumber<<" "<<ProjectileCharge<<G4endl;
    #endif

    const G4ParticleDefinition* Proton = G4Proton::Proton(); 
    // Interaction on P
    G4double XtotPP = csGGinstance->GetTotalIsotopeCrossSection(Proton, KineticEnergy, 1, 1);
    G4double XelPP  = csGGinstance->GetElasticIsotopeCrossSection(Proton, KineticEnergy, 1, 1);

    const G4ParticleDefinition* Neutron = G4Neutron::Neutron();
    // Interaction on N
    G4double XtotPN = csGGinstance->GetTotalIsotopeCrossSection(Neutron, KineticEnergy, 0, 1);
    G4double XelPN  = csGGinstance->GetElasticIsotopeCrossSection(Neutron, KineticEnergy, 0, 1);

    #ifdef debugFTFparams
    G4cout << "XsPP (total and elastic) " << XtotPP/millibarn << " " << XelPP/millibarn <<" (mb)"<< G4endl
           << "XsPN (total and elastic) " << XtotPN/millibarn << " " << XelPN/millibarn <<" (mb)"<< G4endl;
    #endif

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

    Xannihilation = 0.0;
    Xtotal /= millibarn;
    Xelastic /= millibarn;
  }

  // ---------- The projectile is anti-baryon or anti-nucleus ----------------
  //                     anti  Sigma^0_c                  anti Delta^-
  if ( ProjectilePDGcode >= -4112  &&  ProjectilePDGcode <= -1114 ) {
    // Only non-strange and strange baryons are considered

    #ifdef debugFTFparams
    G4cout<<"Projectile is a anti-baryon or anti-nucleus  - "<<ProjectileBaryonNumber<<" "<<ProjectileCharge<<G4endl;
    G4cout<<"(Only non-strange and strange baryons are considered)"<<G4endl;
    #endif

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
      G4double Xasmpt = 36.04 + 0.304*LogS*LogS;  // mb
      LogS = G4Log( SqrtS / 20.74 );
      G4double Basmpt = 11.92 + 0.3036*LogS*LogS;  // GeV^(-2)
      G4double R0 = std::sqrt( 0.40874044*Xasmpt - Basmpt );  // GeV^(-1)

      G4double FlowF = SqrtS / std::sqrt( ECMSsqr*ECMSsqr + ProjectileMass2*ProjectileMass2 +
                                          TargetMass2*TargetMass2 - 2.0*ECMSsqr*ProjectileMass2
                                          - 2.0*ECMSsqr*TargetMass2 
                                          - 2.0*ProjectileMass2*TargetMass2 );

      Xtotal = Xasmpt * ( 1.0 + 13.55*FlowF/R0/R0/R0*
                                (1.0 - 4.47/SqrtS + 12.38/ECMSsqr - 12.43/SqrtS/ECMSsqr) );  // mb

      Xasmpt = 4.4 + 0.101*LogS*LogS;  // mb
      Xelastic = Xasmpt * ( 1.0 + 59.27*FlowF/R0/R0/R0*
                                  (1.0 - 6.95/SqrtS + 23.54/ECMSsqr - 25.34/SqrtS/ECMSsqr ) ); // mb

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
           << " " << Xtotal - Xelastic << " " << Xtotal - Xelastic - Xannihilation << " (mb)"<< G4endl
           << "Plab Xelastic/Xtotal,  Xann/Xin " << Plab << " " << Xelastic/Xtotal << " " 
           << Xannihilation/(Xtotal - Xelastic) << G4endl;
    #endif

  }

  if ( Xtotal == 0.0 ) {  // Projectile is undefined, Nucleon assumed

    const G4ParticleDefinition* Proton = G4Proton::Proton(); 
    // Interaction on P
    G4double XtotPP = csGGinstance->GetTotalIsotopeCrossSection(Proton, KineticEnergy, 1, 1);
    G4double XelPP  = csGGinstance->GetElasticIsotopeCrossSection(Proton, KineticEnergy, 1, 1);

    // Interaction on N
    G4double XtotPN = csGGinstance->GetTotalIsotopeCrossSection(Proton, KineticEnergy, 0, 1);
    G4double XelPN  = csGGinstance->GetElasticIsotopeCrossSection(Proton, KineticEnergy, 0, 1);

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

  // Interactions with elastic and inelastic collisions
  SetProbabilityOfElasticScatt( Xtotal, Xelastic );

  SetRadiusOfHNinteractions2( Xtotal/pi/10.0 );

  if ( ( Xtotal - Xelastic ) == 0.0 ) {
    SetProbabilityOfAnnihilation( 0.0 );
  } else {
    SetProbabilityOfAnnihilation( Xannihilation / (Xtotal - Xelastic) );
  }

  if(Xelastic > 0.0) {
    SetSlope( Xtotal*Xtotal/16.0/pi/Xelastic/0.3894 );// Slope parameter of elastic scattering
                                                      // (GeV/c)^(-2))
    // Parameters of elastic scattering
    // Gaussian parametrization of elastic scattering amplitude assumed
    SetAvaragePt2ofElasticScattering( 1.0/( Xtotal*Xtotal/16.0/pi/Xelastic/0.3894 )*GeV*GeV );
  } else {
    SetSlope(1.0);
    SetAvaragePt2ofElasticScattering( 0.0);
  }
  SetGamma0( GetSlope()*Xtotal/10.0/2.0/pi );

  G4double Xinel = Xtotal - Xelastic;

  #ifdef debugFTFparams
  G4cout<< "Slope of hN elastic scattering" << GetSlope() << G4endl;
  G4cout << "AvaragePt2ofElasticScattering " << GetAvaragePt2ofElasticScattering() << G4endl;
  G4cout<<"Parameters of excitation for projectile "<<ProjectilePDGcode<< G4endl;
  #endif

  if ( (ProjectilePDGcode == 2212) || (ProjectilePDGcode == 2112) ) {  // Projectile is proton or neutron

    // A process probability is parameterized as Prob = A_1*exp(-A_2*y) + A_3*exp(-A_4*y) + A_top
    // y is a rapidity of a partcle in the target nucleus. Ymin is a minimal rapidity below it X=0

    //        Proc#   A1      B1            A2       B2   A3   Atop       Ymin
    /* original hadr-string-diff-V10-03-07 (similar to 10.3.x) 
    SetParams( 0,     13.71, 1.75,          -214.5, 4.25, 0.0, 0.5  ,     1.1 );  // Qexchange without Exc.
    SetParams( 1,      25.0, 1.0,           -50.34, 1.5 , 0.0, 0.0  ,     1.4 );  // Qexchange with    Exc.
    SetParams( 2, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,     0.93);  // Projectile diffraction
    SetParams( 3, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,     0.93);  // Target diffraction
    SetParams( 4,       1.0, 0.0 ,          -2.01 , 0.5 , 0.0, 0.0  ,     1.4 );  // Qexchange with Exc. Additional multiplier
    */

    //         Proc#
    SetParams( 0,     fParCollBaryonProj.GetProc0A1(), fParCollBaryonProj.GetProc0B1(),          
		      fParCollBaryonProj.GetProc0A2(), fParCollBaryonProj.GetProc0B2(), 
		      fParCollBaryonProj.GetProc0A3(), fParCollBaryonProj.GetProc0Atop(),     
		      fParCollBaryonProj.GetProc0Ymin() );   // Qexchange without Exc.
    SetParams( 1,     fParCollBaryonProj.GetProc1A1(), fParCollBaryonProj.GetProc1B1(),          
		      fParCollBaryonProj.GetProc1A2(), fParCollBaryonProj.GetProc1B2(), 
		      fParCollBaryonProj.GetProc1A3(), fParCollBaryonProj.GetProc1Atop(),     
		      fParCollBaryonProj.GetProc1Ymin() );   // Qexchange with Exc.
    if ( Xinel > 0.0 ) {
      SetParams( 2, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,     0.93 );  // Projectile diffraction
      SetParams( 3, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,     0.93 );  // Target diffraction

      SetParams( 4,     fParCollBaryonProj.GetProc4A1(), fParCollBaryonProj.GetProc4B1(),          
		        fParCollBaryonProj.GetProc4A2(), fParCollBaryonProj.GetProc4B2(), 
		        fParCollBaryonProj.GetProc4A3(), fParCollBaryonProj.GetProc4Atop(),     
		        fParCollBaryonProj.GetProc4Ymin() );  // Qexchange with Exc. Additional multiplier
    } else {  // if Xinel=0., zero everything out (obviously)
      SetParams( 2, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0  ,     0.0);
      SetParams( 3, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0  ,     0.0);
      SetParams( 4, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0  ,     0.0);
    }

    if ( AbsProjectileBaryonNumber > 10  ||  NumberOfTargetNucleons > 10 ) {
      //
      // It is not decided what to do with diffraction dissociation in Had-Nucl and Nucl-Nucl interactions
      // For the moment both ProjDiffDisso & TgtDiffDisso for A > 10 are set to false,
      // so both projectile and target diffraction are turned OFF
      //
      if ( ! fParCollBaryonProj.IsProjDiffDissociation() )
         SetParams( 2,       0.0, 0.0 ,           0.0  , 0.0 , 0.0, 0.0   , -100.0  );  // Projectile diffraction
      if ( ! fParCollBaryonProj.IsTgtDiffDissociation() )
         SetParams( 3,       0.0, 0.0 ,           0.0  , 0.0 , 0.0, 0.0   , -100.0  );  // Target diffraction
    }

    SetDeltaProbAtQuarkExchange( fParCollBaryonProj.GetDeltaProbAtQuarkExchange() );

    if ( NumberOfTargetNucleons > 26 ) {
      SetProbOfSameQuarkExchange( 1.0 );
    } else {
      SetProbOfSameQuarkExchange( fParCollBaryonProj.GetProbOfSameQuarkExchange() );
    }

    SetProjMinDiffMass(    fParCollBaryonProj.GetProjMinDiffMass() );     // GeV 
    SetProjMinNonDiffMass( fParCollBaryonProj.GetProjMinNonDiffMass() );  // GeV 

    SetTarMinDiffMass(     fParCollBaryonProj.GetTgtMinDiffMass() );      // GeV
    SetTarMinNonDiffMass(  fParCollBaryonProj.GetTgtMinNonDiffMass() );   // GeV 

    SetAveragePt2(         fParCollBaryonProj.GetAveragePt2() );          // GeV^2       
    SetProbLogDistrPrD(    fParCollBaryonProj.GetProbLogDistrPrD() ); 
    SetProbLogDistr(       fParCollBaryonProj.GetProbLogDistr() ); 

  } else if ( ProjectilePDGcode == -2212  ||  ProjectilePDGcode == -2112 ) {  // Projectile is anti_proton or anti_neutron

    //        Proc#   A1      B1            A2       B2   A3   Atop       Ymin
    SetParams( 0,      0.0 , 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,  1000.0  );  // Qexchange without Exc. 
    SetParams( 1,      0.0 , 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,  1000.0  );  // Qexchange with    Exc.
    if ( Xinel > 0.) {
      SetParams( 2, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,     0.93 );  // Projectile diffraction
      SetParams( 3, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,     0.93 );  // Target diffraction
      SetParams( 4,       1.0, 0.0 ,             0.0, 0.0 , 0.0, 0.0  ,     0.93 );  // Qexchange with    Exc. Additional multiply
    } else {
      SetParams( 2, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0  ,     0.0 );
      SetParams( 3, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0  ,     0.0 );
      SetParams( 4, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0  ,     0.0 );
    }

    if ( AbsProjectileBaryonNumber > 10  ||  NumberOfTargetNucleons > 10 ) {
      //
      // It is not decided what to do with diffraction dissociation in Had-Nucl and Nucl-Nucl interactions
      // For the moment both ProjDiffDisso & TgtDiffDisso are set to false,
      // so both projectile and target diffraction are turned OFF
      //
      if ( ! fParCollBaryonProj.IsProjDiffDissociation() )
         SetParams( 2,       0.0, 0.0 ,           0.0  , 0.0 , 0.0, 0.0   , -100.0  );  // Projectile diffraction
      if ( ! fParCollBaryonProj.IsTgtDiffDissociation() )
         SetParams( 3,       0.0, 0.0 ,           0.0  , 0.0 , 0.0, 0.0   , -100.0  );  // Target diffraction
    }

    SetDeltaProbAtQuarkExchange( 0.0 );
    SetProbOfSameQuarkExchange( 0.0 );
    SetProjMinDiffMass( ProjectileMass + 0.22 );     // GeV 
    SetProjMinNonDiffMass( ProjectileMass + 0.22 );  // GeV
    SetTarMinDiffMass( TargetMass + 0.22 );          // GeV
    SetTarMinNonDiffMass( TargetMass + 0.22 );       // GeV
    SetAveragePt2( 0.3 );                            // GeV^2
    SetProbLogDistrPrD( 0.55 );
    SetProbLogDistr( 0.55 );
            
  } else if ( ProjectileabsPDGcode == 211  ||  ProjectilePDGcode ==  111 ) {  // Projectile is Pion 

    //        Proc#   A1      B1            A2       B2      A3   Atop        Ymin
    /* --> original code
    SetParams( 0,  150.0,    1.8 ,       -247.3,     2.3,    0.,   1. ,       2.3 ); 
    SetParams( 1,   5.77,    0.6 ,        -5.77,     0.8,    0.,   0. ,       0.0 );
    SetParams( 2,   2.27,    0.5 ,     -98052.0,     4.0,    0.,   0. ,       3.0 ); 
    SetParams( 3,    7.0,    0.9,        -85.28,     1.9,  0.08,   0. ,       2.2 );
    SetParams( 4,    1.0,    0.0 ,       -11.02,     1.0,   0.0,   0. ,       2.4 );  // Qexchange with    Exc. Additional multiply
    */
    //         Proc#
    SetParams( 0,     fParCollPionProj.GetProc0A1(), fParCollPionProj.GetProc0B1(),          
		      fParCollPionProj.GetProc0A2(), fParCollPionProj.GetProc0B2(), 
		      fParCollPionProj.GetProc0A3(), fParCollPionProj.GetProc0Atop(),     
		      fParCollPionProj.GetProc0Ymin() );   // Qexchange without Exc.
    SetParams( 1,     fParCollPionProj.GetProc1A1(), fParCollPionProj.GetProc1B1(),          
		      fParCollPionProj.GetProc1A2(), fParCollPionProj.GetProc1B2(), 
		      fParCollPionProj.GetProc1A3(), fParCollPionProj.GetProc1Atop(),     
		      fParCollPionProj.GetProc1Ymin() );   // Qexchange with Exc.
    SetParams( 2,     fParCollPionProj.GetProc2A1(), fParCollPionProj.GetProc2B1(),          
		      fParCollPionProj.GetProc2A2(), fParCollPionProj.GetProc2B2(), 
		      fParCollPionProj.GetProc2A3(), fParCollPionProj.GetProc2Atop(),     
		      fParCollPionProj.GetProc2Ymin() );   // Projectile diffraction
    SetParams( 3,     fParCollPionProj.GetProc3A1(), fParCollPionProj.GetProc3B1(),          
		      fParCollPionProj.GetProc3A2(), fParCollPionProj.GetProc3B2(), 
		      fParCollPionProj.GetProc3A3(), fParCollPionProj.GetProc3Atop(),     
		      fParCollPionProj.GetProc3Ymin() );   // Target diffraction
    SetParams( 4,     fParCollPionProj.GetProc4A1(), fParCollPionProj.GetProc4B1(),          
		      fParCollPionProj.GetProc4A2(), fParCollPionProj.GetProc4B2(), 
		      fParCollPionProj.GetProc4A3(), fParCollPionProj.GetProc4Atop(),     
		      fParCollPionProj.GetProc4Ymin() );   // Qexchange with Exc. Additional multiply

    // NOTE: how can it be |ProjectileBaryonNumber| > 10 if projectile is a pion ??? 
    //
    if ( AbsProjectileBaryonNumber > 10  ||  NumberOfTargetNucleons > 10 ) {
       if ( ! fParCollPionProj.IsProjDiffDissociation() )
          SetParams( 2,      0.0 , 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,  -100.0 );  // Projectile diffraction
       if ( ! fParCollPionProj.IsTgtDiffDissociation() )
          SetParams( 3,      0.0 , 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,  -100.0 );  // Target diffraction
    }

    /* original code -->
    SetDeltaProbAtQuarkExchange( 0.56 );
    SetProjMinDiffMass( 1.0 );            // GeV
    SetProjMinNonDiffMass( 1.0 );         // GeV 
    SetTarMinDiffMass( 1.16 );            // GeV
    SetTarMinNonDiffMass( 1.16 );         // GeV
    SetAveragePt2( 0.3 );                 // GeV^2
    SetProbLogDistrPrD( 0.55 );
    SetProbLogDistr( 0.55 );
    */

    // JVY update, Aug.8, 2018 --> Feb.14, 2019
    //
    SetDeltaProbAtQuarkExchange( fParCollPionProj.GetDeltaProbAtQuarkExchange() );
    SetProjMinDiffMass( fParCollPionProj.GetProjMinDiffMass() );                  // GeV 
    SetProjMinNonDiffMass( fParCollPionProj.GetProjMinNonDiffMass() );            // GeV 
    SetTarMinDiffMass( fParCollPionProj.GetTgtMinDiffMass() );                    // GeV
    SetTarMinNonDiffMass( fParCollPionProj.GetTgtMinNonDiffMass() );              // GeV 
    SetAveragePt2( fParCollPionProj.GetAveragePt2() );                            // GeV^2       
    SetProbLogDistrPrD( fParCollPionProj.GetProbLogDistrPrD() ); 
    SetProbLogDistr( fParCollPionProj.GetProbLogDistr() ); 
    
    // ---> end update

  } else if ( ProjectileabsPDGcode == 321  ||  ProjectileabsPDGcode == 311  || 
              ProjectilePDGcode == 130     ||  ProjectilePDGcode == 310 ) {  // Projectile is Kaon

    //        Proc#   A1      B1            A2       B2   A3   Atop       Ymin
    SetParams( 0,     60.0 , 2.5 ,           0.0  , 0.0 , 0.0, 0.0  ,  -100.0  );  // Qexchange without Exc. 
    SetParams( 1,      6.0 , 1.0 ,         -24.33 , 2.0 , 0.0, 0.0  ,     1.40 );  // Qexchange with    Exc.
    SetParams( 2,      2.76, 1.2 ,         -22.5  , 2.7 ,0.04, 0.0  ,     1.40 );  // Projectile diffraction
    SetParams( 3,      1.09, 0.5 ,          -8.88 , 2.  ,0.05, 0.0  ,     1.40 );  // Target diffraction
    SetParams( 4,       1.0, 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,     0.93 );  // Qexchange with    Exc. Additional multiply
    if ( AbsProjectileBaryonNumber > 10  ||  NumberOfTargetNucleons > 10 ) {
      SetParams( 2,      0.0 , 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,  -100.0  );  // Projectile diffraction
      SetParams( 3,      0.0 , 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,  -100.0  );  // Target diffraction
    }

    SetDeltaProbAtQuarkExchange( 0.6 );
    SetProjMinDiffMass( 0.7 );     // GeV 
    SetProjMinNonDiffMass( 0.7 );  // GeV 
    SetTarMinDiffMass( 1.16 );     // GeV
    SetTarMinNonDiffMass( 1.16 );  // GeV
    SetAveragePt2( 0.3 );          // GeV^2
    SetProbLogDistrPrD( 0.55 );
    SetProbLogDistr( 0.55 );

  } else {  // Projectile is not p, n, Pi0, Pi+, Pi-, K+, K-, K0 or their anti-particles

    if ( ProjectileabsPDGcode > 1000 ) {  // The projectile is a baryon as P or N
      //        Proc#   A1      B1            A2       B2   A3   Atop       Ymin
      SetParams( 0,     13.71, 1.75,          -30.69, 3.0 , 0.0, 1.0  ,     0.93 ); // Qexchange without Exc.
      SetParams( 1,      25.0, 1.0,          -50.34,  1.5 , 0.0, 0.0  ,     1.4 );  // Qexchange with    Exc.
      if ( Xinel > 0.) {
        SetParams( 2, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,   0.93);  // Projectile diffraction
        SetParams( 3, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,   0.93);  // Target diffraction
        SetParams( 4,       1.0, 0.0 ,          -2.01 , 0.5 , 0.0, 0.0  ,   1.4 );  // Qexchange with    Exc. Additional multiply
      } else {
        SetParams( 2, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0  ,     0.0);
        SetParams( 3, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0  ,     0.0);
        SetParams( 4, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0  ,     0.0);
      }

    } else {  // The projectile is a meson as K+-0
      //        Proc#   A1      B1            A2       B2   A3   Atop       Ymin
      SetParams( 0,     60.0 , 2.5 ,           0.0  , 0.0 , 0.0, 0.0  ,  -100.0  );  // Qexchange without Exc. 
      SetParams( 1,      6.0 , 1.0 ,         -24.33 , 2.0 , 0.0, 0.0  ,     1.40 );  // Qexchange with    Exc.
      SetParams( 2,      2.76, 1.2 ,         -22.5  , 2.7 ,0.04, 0.0  ,     1.40 );  // Projectile diffraction
      SetParams( 3,      1.09, 0.5 ,          -8.88 , 2.  ,0.05, 0.0  ,     1.40 );  // Target diffraction
      SetParams( 4,       1.0, 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,     0.93 );  // Qexchange with    Exc. Additional multiply
    }

    if ( AbsProjectileBaryonNumber > 10  ||  NumberOfTargetNucleons > 10 ) {
      SetParams( 2,      0.0 , 0.0 ,            0.0 , 0.0 , 0.0, 0.0  ,  -100.0  );  // Projectile diffraction
      SetParams( 3,      0.0 , 0.0 ,            0.0 , 0.0 , 0.0, 0.0  ,  -100.0  );  // Target diffraction
    }

    SetDeltaProbAtQuarkExchange( 0.0 );
    SetProbOfSameQuarkExchange( 0.0 );

    SetProjMinDiffMass(    GetMinMass(particle)/GeV );
    SetProjMinNonDiffMass( GetMinMass(particle)/GeV );

    const G4ParticleDefinition* Neutron = G4Neutron::Neutron();
    SetTarMinDiffMass(    GetMinMass(Neutron)/GeV );
    SetTarMinNonDiffMass( GetMinMass(Neutron)/GeV );

    SetAveragePt2( 0.3 );  // GeV^2
    SetProbLogDistrPrD( 0.55 );
    SetProbLogDistr( 0.55 );

  }

  #ifdef debugFTFparams
  G4cout<<"DeltaProbAtQuarkExchange "<< GetDeltaProbAtQuarkExchange() << G4endl;
  G4cout<<"ProbOfSameQuarkExchange  "<< GetProbOfSameQuarkExchange()  << G4endl;
  G4cout<<"ProjMinDiffMass          "<< GetProjMinDiffMass()/GeV <<" GeV"<< G4endl;
  G4cout<<"ProjMinNonDiffMass       "<< GetProjMinNonDiffMass()  <<" GeV"<< G4endl;
  G4cout<<"TarMinDiffMass           "<< GetTarMinDiffMass()      <<" GeV"<< G4endl;
  G4cout<<"TarMinNonDiffMass        "<< GetTarMinNonDiffMass()   <<" GeV"<< G4endl;
  G4cout<<"AveragePt2               "<< GetAveragePt2()          <<" GeV^2"<< G4endl;
  G4cout<<"ProbLogDistrPrD          "<< GetProbLogDistrPrD() << G4endl;
  G4cout<<"ProbLogDistrTrD          "<< GetProbLogDistr()    << G4endl;
  #endif

  // Set parameters of nuclear destruction

  if ( ProjectileabsPDGcode < 1000 ) {  // Meson projectile

    SetMaxNumberOfCollisions( Plab, 2.0 );  //  3.0 )
    //
    // target destruction
    //
    /* original code --->
    SetCofNuclearDestruction( 0.00481*G4double(NumberOfTargetNucleons)*
            G4Exp( 4.0*(Ylab - 2.1) )/( 1.0 + G4Exp( 4.0*(Ylab - 2.1) ) ) );
    
    SetR2ofNuclearDestruction( 1.5*fermi*fermi );
    SetDofNuclearDestruction( 0.3 );
    SetPt2ofNuclearDestruction( ( 0.035 + 0.04*G4Exp( 4.0*(Ylab - 2.5) )/
                                         ( 1.0 + G4Exp( 4.0*(Ylab - 2.5) ) ) )*GeV*GeV );
    SetMaxPt2ofNuclearDestruction( 1.0*GeV*GeV );
    SetExcitationEnergyPerWoundedNucleon( 40.0*MeV );
    */
    double coeff = fParCollMesonProj.GetNuclearTgtDestructP1();
    // 
    // NOTE (JVY): Set this switch to false/true on line 138
    //
    if ( fParCollMesonProj.IsNuclearTgtDestructP1_ADEP() )       
    {                                                             
      coeff *= G4double(NumberOfTargetNucleons);                  
    }                                                             
    double exfactor = G4Exp( fParCollMesonProj.GetNuclearTgtDestructP2()
                           * (Ylab-fParCollMesonProj.GetNuclearTgtDestructP3()) );
    coeff *= exfactor;
    coeff /= ( 1.+ exfactor );

    SetCofNuclearDestruction( coeff );

    SetR2ofNuclearDestruction( fParCollMesonProj.GetR2ofNuclearDestruct() );
    SetDofNuclearDestruction( fParCollMesonProj.GetDofNuclearDestruct() );
    coeff = fParCollMesonProj.GetPt2NuclearDestructP2();
    exfactor = G4Exp( fParCollMesonProj.GetPt2NuclearDestructP3()
                    * (Ylab-fParCollMesonProj.GetPt2NuclearDestructP4()) );
    coeff *= exfactor;
    coeff /= ( 1. + exfactor );
    SetPt2ofNuclearDestruction( (fParCollMesonProj.GetPt2NuclearDestructP1()+coeff)*CLHEP::GeV*CLHEP::GeV ); 

    SetMaxPt2ofNuclearDestruction( fParCollMesonProj.GetMaxPt2ofNuclearDestruct() );
    SetExcitationEnergyPerWoundedNucleon( fParCollMesonProj.GetExciEnergyPerWoundedNucleon() ); 

  } else if ( ProjectilePDGcode == -2212  ||  ProjectilePDGcode == -2112 ) {  // for anti-baryon projectile

    SetMaxNumberOfCollisions( Plab, 2.0 );

    SetCofNuclearDestruction( 0.00481*G4double(NumberOfTargetNucleons)*
           G4Exp( 4.0*(Ylab - 2.1) )/( 1.0 + G4Exp( 4.0*(Ylab - 2.1) ) ) );
    SetR2ofNuclearDestruction( 1.5*fermi*fermi );
    SetDofNuclearDestruction( 0.3 );
    SetPt2ofNuclearDestruction( ( 0.035 + 0.04*G4Exp( 4.0*(Ylab - 2.5) )/
                                         ( 1.0 + G4Exp( 4.0*(Ylab - 2.5) ) ) )*GeV*GeV );
    SetMaxPt2ofNuclearDestruction( 1.0*GeV*GeV );
    SetExcitationEnergyPerWoundedNucleon( 40.0*MeV );
    if ( Plab < 2.0 ) {  // 2 GeV/c
      // For slow anti-baryon we have to garanty putting on mass-shell
      SetCofNuclearDestruction( 0.0 );
      SetR2ofNuclearDestruction( 1.5*fermi*fermi ); // this is equivalent to setting a few line above
                                                    // is it even necessary ?
      SetDofNuclearDestruction( 0.01 );
      SetPt2ofNuclearDestruction( 0.035*GeV*GeV );
      SetMaxPt2ofNuclearDestruction( 0.04*GeV*GeV );
    }

  } else {  // Projectile baryon assumed

    // NOTE (JVY) FIXME !!! Will decide later how/if to make this one configurable...
    //
    SetMaxNumberOfCollisions( Plab, 2.0 );

    // projectile destruction - does NOT really matter for particle projectile, only for a nucleus projectile
    //
    double coeff = 0.;
    coeff = fParCollBaryonProj.GetNuclearProjDestructP1();
    // 
    // NOTE (JVY): Set this switch to false/true on line 136
    //
    if ( fParCollBaryonProj.IsNuclearProjDestructP1_NBRNDEP() )   
    {                                                             
      coeff *= G4double(AbsProjectileBaryonNumber);               
    }                                                             
    double exfactor = G4Exp( fParCollBaryonProj.GetNuclearProjDestructP2()*(Ylab-fParCollBaryonProj.GetNuclearProjDestructP3()) ); 
    coeff *= exfactor;
    coeff /= ( 1.+ exfactor );
    SetCofNuclearDestructionPr( coeff );

    // target desctruction
    //
    coeff = fParCollBaryonProj.GetNuclearTgtDestructP1();
    // 
    // NOTE (JVY): Set this switch to false/true on line 138
    //
    if ( fParCollBaryonProj.IsNuclearTgtDestructP1_ADEP() )       
    {                                                             
      coeff *= G4double(NumberOfTargetNucleons);                  
    }                                                             
    exfactor = G4Exp( fParCollBaryonProj.GetNuclearTgtDestructP2()*(Ylab-fParCollBaryonProj.GetNuclearTgtDestructP3()) );
    coeff *= exfactor;
    coeff /= ( 1.+ exfactor );
    SetCofNuclearDestruction(   coeff );

    SetR2ofNuclearDestruction( fParCollBaryonProj.GetR2ofNuclearDestruct() );
    SetDofNuclearDestruction( fParCollBaryonProj.GetDofNuclearDestruct() );

    coeff = fParCollBaryonProj.GetPt2NuclearDestructP2();
    exfactor = G4Exp( fParCollBaryonProj.GetPt2NuclearDestructP3()*(Ylab-fParCollBaryonProj.GetPt2NuclearDestructP4())  );
    coeff *= exfactor;
    coeff /= ( 1. + exfactor );
    SetPt2ofNuclearDestruction( (fParCollBaryonProj.GetPt2NuclearDestructP1()+coeff)*CLHEP::GeV*CLHEP::GeV ); 

    SetMaxPt2ofNuclearDestruction( fParCollBaryonProj.GetMaxPt2ofNuclearDestruct() );
    SetExcitationEnergyPerWoundedNucleon( fParCollBaryonProj.GetExciEnergyPerWoundedNucleon() ); 

  }

  #ifdef debugFTFparams
  G4cout<<"CofNuclearDestructionPr           "<< GetCofNuclearDestructionPr() << G4endl;
  G4cout<<"CofNuclearDestructionTr           "<< GetCofNuclearDestruction()   << G4endl;
  G4cout<<"R2ofNuclearDestruction            "<< GetR2ofNuclearDestruction()/fermi/fermi  <<" fermi^2"<< G4endl;
  G4cout<<"DofNuclearDestruction             "<< GetDofNuclearDestruction()  << G4endl;
  G4cout<<"Pt2ofNuclearDestruction           "<< GetPt2ofNuclearDestruction()/GeV/GeV <<" GeV^2"<< G4endl;
  G4cout<<"ExcitationEnergyPerWoundedNucleon "<< GetExcitationEnergyPerWoundedNucleon()   <<" MeV"<< G4endl;
  #endif

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
  //SetProbabilityOfElasticScatt( 1.0, 1.0);                            //(Xtotal, Xelastic);
  //SetProbabilityOfProjDiff( 1.0*0.62*G4Pow::GetInstance()->powA( s/GeV/GeV, -0.51 ) );  // 0->1
  //SetProbabilityOfTarDiff( 4.0*0.62*G4Pow::GetInstance()->powA( s/GeV/GeV, -0.51 ) );   // 2->4
  //SetAveragePt2( 0.3 );                                               // (0.15)
  //SetAvaragePt2ofElasticScattering( 0.0 );

  //SetMaxNumberOfCollisions( Plab, 6.0 ); //(4.0*(Plab + 0.01), Plab); // 6.0 );
  //SetAveragePt2( 0.15 );
  //SetCofNuclearDestruction(-1.);//( 0.75 );                           // (0.25)             
  //SetExcitationEnergyPerWoundedNucleon(0.);//( 30.0*MeV );            // (75.0*MeV) 
  //SetDofNuclearDestruction(0.);//( 0.2 ); //0.4                       // 0.3 0.5

  /*
  SetAveragePt2(0.3);
  SetCofNuclearDestructionPr(0.);
  SetCofNuclearDestruction(0.);             //( 0.5 ); (0.25)             
  SetExcitationEnergyPerWoundedNucleon(0.); // 30.0*MeV; (75.0*MeV) 
  SetDofNuclearDestruction(0.);             // 0.2; 0.4; 0.3; 0.5
  SetPt2ofNuclearDestruction(0.);           //(2.*0.075*GeV*GeV); ( 0.3*GeV*GeV ); (0.168*GeV*GeV) 
  */

  //SetExcitationEnergyPerWoundedNucleon(0.001);
  //SetPt2Kink( 0.0*GeV*GeV );

  //SetRadiusOfHNinteractions2( Xtotal/pi/10.0 /2.);
  //SetRadiusOfHNinteractions2( (Xtotal - Xelastic)/pi/10.0 );
  //SetProbabilityOfElasticScatt( 1.0, 0.0);
  /*
  G4cout << "Pt2 " << GetAveragePt2()<<" "<<GetAveragePt2()/GeV/GeV<<G4endl;
  G4cout << "Cnd " << GetCofNuclearDestruction() << G4endl;
  G4cout << "Dnd " << GetDofNuclearDestruction() << G4endl;
  G4cout << "Pt2 " << GetPt2ofNuclearDestruction()/GeV/GeV << G4endl;
  */

} 

//============================================================================

G4double G4FTFParameters::GetMinMass( const G4ParticleDefinition* aParticle ) {
   // The code is used for estimating the minimal string mass produced in diffraction dissociation.
   // The indices used for minMassQDiQStr must be between 1 and 5, corresponding to the 5 considered
   // quarks: d, u, s, c and b; enforcing this explicitly avoids compilation errors.
   G4double EstimatedMass = 0.0;
   G4int partID = std::abs(aParticle->GetPDGEncoding());
   G4int Qleft  = std::max(  partID/100,     1 );
   G4int Qright = std::max( (partID/ 10)%10, 1 );
   if        ( Qleft < 6  &&  Qright < 6 ) {  // Q-Qbar string
     EstimatedMass = StringMass->minMassQQbarStr[Qleft-1][Qright-1];
   } else if ( Qleft < 6  &&  Qright > 6 ) {  // Q - DiQ string
     G4int q1 = std::max( std::min( Qright/10, 5 ), 1 );
     G4int q2 = std::max( std::min( Qright%10, 5 ), 1 );
     EstimatedMass = StringMass->minMassQDiQStr[Qleft-1][q1-1][q2-1];
   } else if ( Qleft > 6  &&  Qright < 6 ) {  // DiQ - Q string
     G4int q1 = std::max( std::min( Qleft/10, 5 ), 1 );
     G4int q2 = std::max( std::min( Qleft%10, 5 ), 1 );
     EstimatedMass = StringMass->minMassQDiQStr[Qright-1][q1-1][q2-1];
   }
   return EstimatedMass;
}

//============================================================================

G4double G4FTFParameters::GetProcProb( const G4int ProcN, const G4double y ) {
  G4double Prob( 0.0 );
  if ( y < ProcParams[ProcN][6] ) {
    Prob = ProcParams[ProcN][5]; 
    if (Prob < 0.) Prob=0.;
    return Prob;
  }
  Prob = ProcParams[ProcN][0] * G4Exp( -ProcParams[ProcN][1]*y ) +
         ProcParams[ProcN][2] * G4Exp( -ProcParams[ProcN][3]*y ) +
         ProcParams[ProcN][4];
  if (Prob < 0.) Prob=0.;
  return Prob;
}

//============================================================================

G4HadronicDeveloperParameters& HDP = G4HadronicDeveloperParameters::GetInstance();


class G4FTFSettingDefaultHDP 
{
  public:   

    // ctor
    G4FTFSettingDefaultHDP() {

    //==================================================================
    //
    // Cross sections for elementary processes
    //
    // these are for Inelastic interactions, i.e. Xinelastic=(Xtotal-Xelastix)>0.
    // for elastic, all the A's & B's, Atop & Ymin are zeros
    // general formula: Pp = A1*exp(B1*Y) + A2*exp(B2*Y) + A3
    // but if Y<Ymin, then Pp=max(0.,Atop)
    // for details, see also G4FTFParameters::GetProcProb( ProcN, y )
    //
    // baryons
    //
    /* JVY, Oct. 31, 2017: Per Alberto R. & Vladimir U., keep this group of parameters FIXED
    // Process=0 --> Qexchg w/o excitation
    //
    HDP.SetDefault( "FTF_BARYON_PROC0_A1",  13.71 );
    HDP.SetDefault( "FTF_BARYON_PROC0_B1",   1.75 );
    HDP.SetDefault( "FTF_BARYON_PROC0_A2", -30.69 ); 
    HDP.SetDefault( "FTF_BARYON_PROC0_B2",   3.0  ); 
    HDP.SetDefault( "FTF_BARYON_PROC0_A3",   0.0  );
    HDP.SetDefault( "FTF_BARYON_PROC0_ATOP", 1.0  ); 
    HDP.SetDefault( "FTF_BARYON_PROC0_YMIN", 0.93 ); 
    //
    // Process=1 --> Qexchg w/excitation
    //
    HDP.SetDefault( "FTF_BARYON_PROC1_A1",  25.   );
    HDP.SetDefault( "FTF_BARYON_PROC1_B1",   1.   );
    HDP.SetDefault( "FTF_BARYON_PROC1_A2", -50.34 );
    HDP.SetDefault( "FTF_BARYON_PROC1_B2",   1.5  );
    HDP.SetDefault( "FTF_BARYON_PROC1_A3",   0.   );
    HDP.SetDefault( "FTF_BARYON_PROC1_ATOP", 0.   );
    HDP.SetDefault( "FTF_BARYON_PROC1_YMIN", 1.4  );
    */
    //
    // NOTE: Process #2 & 3 are projectile & target diffraction
    //       they have more complex definition of A1 & A2 
    //      (see example below)
    // SetParams( 2, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,     0.93);// Projectile diffraction
    // SetParams( 3, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,     0.93);// Target diffraction
    //
    // Also, for ( AbsProjectileBaryonNumber > 10 ||  NumberOfTargetNucleons > 10 )
    // projectile and/or target diffraction (dissociation) may be switched ON/OFF 
    //
    HDP.SetDefault( "FTF_BARYON_DIFF_DISSO_PROJ", false );
    HDP.SetDefault( "FTF_BARYON_DIFF_DISSO_TGT",  false ); // as in hadr-string-diff-V10-03-07
    //
    /* JVY, Oct. 31, 2017: Per Alberto R. & Vladimir U., keep this group of parameters FIXED
    // Process=4 --> Qexchg w/additional multiplier in excitation 
    //
    HDP.SetDefault( "FTF_BARYON_PROC4_A1",  0.6 ); 
    HDP.SetDefault( "FTF_BARYON_PROC4_B1",  0.  );
    HDP.SetDefault( "FTF_BARYON_PROC4_A2", -1.2 ); 
    HDP.SetDefault( "FTF_BARYON_PROC4_B2",  0.5 );
    HDP.SetDefault( "FTF_BARYON_PROC4_A3",  0.  );
    HDP.SetDefault( "FTF_BARYON_PROC4_ATOP",0.  );
    HDP.SetDefault( "FTF_BARYON_PROC4_YMIN",1.4 );
    */
    //
    // Parameters of participating hadron (baryon) excitation
    //
    HDP.SetDefault( "FTF_BARYON_DELTA_PROB_QEXCHG", 0. );
    HDP.SetDefault( "FTF_BARYON_PROB_SAME_QEXCHG", 0. );
    HDP.SetDefault( "FTF_BARYON_DIFF_M_PROJ", 1.16, 1.16, 3. );     // it's supposed to be in GeV but do NOT do (*CLHEP::GeV) 
                                                                    // because it'll be done in the G4FTFParameters::SetProjMinDiffMass
    HDP.SetDefault( "FTF_BARYON_NONDIFF_M_PROJ", 1.16, 1.16, 3. );  // do NOT (*CLHEP::GeV) - same as above
    HDP.SetDefault( "FTF_BARYON_DIFF_M_TGT", 1.16, 1.16, 3. );      // do NOT (*CLHEP::GeV) - same as above
    HDP.SetDefault( "FTF_BARYON_NONDIFF_M_TGT", 1.16, 1.16, 3. );   // do NOT (*CLHEP::GeV) - same as above
    HDP.SetDefault( "FTF_BARYON_AVRG_PT2", 0.3, 0.08, 1. );         //  do NOT (*CLHEP::GeV*CLHEP::GeV) 
    //
    // JVY, Oct. 6, 2017: Per Alberto R., keep these two settings fixed (for now)
    //
    // HDP.SetDefault( "FTF_BARYON_PROB_DISTR_PROJ", 0.3 ); 
    // HDP.SetDefault( "FTF_BARYON_PROB_DISTR_TGT", 0.3 ); 

    
    // pions
    //
    // JVY, Aug.8, 2018 --> Feb.14, 2019 --> June 25, 2019: 
    // Parameters of participating hadron (pions) excitation
    //
    /* JVY, June 25, 2019: For now, keep this group of parameters FIXED
    // Process=0 --> Qexchg w/o excitation
    //
    HDP.SetDefault( "FTF_PION_PROC0_A1", 150.0  );
    HDP.SetDefault( "FTF_PION_PROC0_B1",   1.8  );
    HDP.SetDefault( "FTF_PION_PROC0_A2",-247.3  ); 
    HDP.SetDefault( "FTF_PION_PROC0_B2",   2.3 ); 
    HDP.SetDefault( "FTF_PION_PROC0_A3",   0.0  );
    HDP.SetDefault( "FTF_PION_PROC0_ATOP", 1.0  ); 
    HDP.SetDefault( "FTF_PION_PROC0_YMIN", 2.3  ); 
    //
    // Process=1 --> Qexchg w/excitation
    //
    HDP.SetDefault( "FTF_PION_PROC1_A1",   5.77 );
    HDP.SetDefault( "FTF_PION_PROC1_B1",   0.6  );
    HDP.SetDefault( "FTF_PION_PROC1_A2",  -5.77 );
    HDP.SetDefault( "FTF_PION_PROC1_B2",   0.8  );
    HDP.SetDefault( "FTF_PION_PROC1_A3",   0.   );
    HDP.SetDefault( "FTF_PION_PROC1_ATOP", 0.   );
    HDP.SetDefault( "FTF_PION_PROC1_YMIN", 0.0  );
    //
    // NOTE: Process #2 & 3 are projectile & target diffraction
    //
    // Process=2 --> Projectile diffraction
    //
    HDP.SetDefault( "FTF_PION_PROC2_A1",    2.27 );
    HDP.SetDefault( "FTF_PION_PROC2_B1",    0.5  );
    HDP.SetDefault( "FTF_PION_PROC2_A2", -98052.0);
    HDP.SetDefault( "FTF_PION_PROC2_B2",    4.0  );
    HDP.SetDefault( "FTF_PION_PROC2_A3",    0.   );
    HDP.SetDefault( "FTF_PION_PROC2_ATOP",  0.   );
    HDP.SetDefault( "FTF_PION_PROC2_YMIN",  3.0  );
    //
    // Process=3 --> Target diffraction
    //
    HDP.SetDefault( "FTF_PION_PROC3_A1",    7.0 );
    HDP.SetDefault( "FTF_PION_PROC3_B1",    0.9 );
    HDP.SetDefault( "FTF_PION_PROC3_A2",  -85.28);
    HDP.SetDefault( "FTF_PION_PROC3_B2",    1.9 );
    HDP.SetDefault( "FTF_PION_PROC3_A3",    0.08);
    HDP.SetDefault( "FTF_PION_PROC3_ATOP",  0.  );
    HDP.SetDefault( "FTF_PION_PROC3_YMIN",  2.2 );
    */
    //
    // projectile and/or target diffraction (dissociation) may be switched ON/OFF 
    //
    // NOTE: Both projectile and target diffraction are turned OFF if
    // a) Number of Target Nucleons > 10 (NumberOfTargetNucleons > 10)
    //    OR
    // b) Projectile Baryon Number > 10 (AbsProjectileBaryonNumber > 10)
    //    ... which is "strange" because projectile is a pion !!!... so it's always OFF
    //
    HDP.SetDefault( "FTF_PION_DIFF_DISSO_PROJ", false );
    HDP.SetDefault( "FTF_PION_DIFF_DISSO_TGT",  false ); 
    //
    /* JVY, June 25, 2019: For now keep this group of parameters FIXED
    // Process=4 --> Qexchg w/additional multiplier in excitation 
    //
    HDP.SetDefault( "FTF_PION_PROC4_A1",  1.0 ); 
    HDP.SetDefault( "FTF_PION_PROC4_B1",  0.  );
    HDP.SetDefault( "FTF_PION_PROC4_A2",-11.02); 
    HDP.SetDefault( "FTF_PION_PROC4_B2",  1.0 );
    HDP.SetDefault( "FTF_PION_PROC4_A3",  0.  );
    HDP.SetDefault( "FTF_PION_PROC4_ATOP",0.  );
    HDP.SetDefault( "FTF_PION_PROC4_YMIN",2.4 );
    */
    //    
    // NOTE; As of geant4-10-05, all these settings beloe are correct 
    //       (and are the same as they were in 10.4.ref06) 
    
    HDP.SetDefault( "FTF_PION_DELTA_PROB_QEXCHG", 0.56 );    // in the past used to be 0.35
    HDP.SetDefault( "FTF_PION_DIFF_M_PROJ", 1.0, 0.5, 3. );  // in the past used to be 0.5...
    HDP.SetDefault( "FTF_PION_NONDIFF_M_PROJ", 1.0, 0.5, 3. ); // so why not set 0.5 as a low limit ?... 
                                                               // ... or perhaps even lower ?
    HDP.SetDefault( "FTF_PION_DIFF_M_TGT", 1.16, 1.16, 3. );   // All (NON)DIFF_M's are supposed to be in GeV but do NOT do (*CLHEP::GeV) 
                                                               // because it'll be done in the G4FTFParameters::SetProjMinDiffMass
    HDP.SetDefault( "FTF_PION_NONDIFF_M_TGT", 1.16, 1.16, 3. );
    HDP.SetDefault( "FTF_PION_AVRG_PT2", 0.3, 0.08, 1. );      //  do NOT (*CLHEP::GeV*CLHEP::GeV)

    //==================================================================
    //
    // nuclear destruction 
    //
    // NOTE: Settings of most of these parameters are the same
    //       for different types of projectile hadron
    //       However, we decided to introduce separate variables
    //       and configuration cards for each type of projectile
    //
    // baryons
    //
    // projectile destruction
    //
    HDP.SetDefault( "FTF_BARYON_NUCDESTR_P1_PROJ", 1., 0., 1. ); // in principle, it should be 1./NBRN - FIXME later !
    HDP.SetDefault( "FTF_BARYON_NUCDESTR_P1_NBRN_PROJ", false );
    //
    // for now, keep fixed p2 & p3 for the proj destruction
    // they're defined explicitly in G4FTFParamCollection ctor
    //
    // target destruction
    //
    HDP.SetDefault( "FTF_BARYON_NUCDESTR_P1_TGT", 1., 0., 1. );   
    HDP.SetDefault( "FTF_BARYON_NUCDESTR_P1_ADEP_TGT", false );          
    HDP.SetDefault( "FTF_BARYON_NUCDESTR_P2_TGT", 4.0, 2., 16. );
    HDP.SetDefault( "FTF_BARYON_NUCDESTR_P3_TGT", 2.1, 0., 4. );
    //
    HDP.SetDefault( "FTF_BARYON_PT2_NUCDESTR_P1", 0.035, 0., 0.25 ); 
    HDP.SetDefault( "FTF_BARYON_PT2_NUCDESTR_P2", 0.04, 0., 0.25 ); 
    HDP.SetDefault( "FTF_BARYON_PT2_NUCDESTR_P3", 4.0, 2., 16. ); 
    HDP.SetDefault( "FTF_BARYON_PT2_NUCDESTR_P4", 2.5, 0., 4. ); 
    //
    HDP.SetDefault( "FTF_BARYON_NUCDESTR_R2", 1.5*CLHEP::fermi*CLHEP::fermi, 0.5*CLHEP::fermi*CLHEP::fermi, 2.*CLHEP::fermi*CLHEP::fermi  );
    HDP.SetDefault( "FTF_BARYON_EXCI_E_PER_WNDNUCLN", 40.*CLHEP::MeV, 0., 100.*CLHEP::MeV );
    HDP.SetDefault( "FTF_BARYON_NUCDESTR_DISP", 0.3, 0.1, 0.4 );
    //
    // JVY, Oct. 6, 2017: Per Alberto R., this is just a technical parameter,
    //                    and it should NOT be changed
    //
    // HDP.SetDefault( "FTF_BARYON_NUCDESTR_MAXPT2", 1. * CLHEP::GeV*CLHEP::GeV  ); 	 


    // mesons - these parameters are common for pions, kaons, etc. (per original code)
    //
    // NOTE: *NO* projectile destruction for mesons !!!
    //
    // target destruction
    //
    HDP.SetDefault( "FTF_MESON_NUCDESTR_P1_TGT", 0.00481, 0., 1. );    
    HDP.SetDefault( "FTF_MESON_NUCDESTR_P1_ADEP_TGT", true );           
    HDP.SetDefault( "FTF_MESON_NUCDESTR_P2_TGT", 4.0, 2., 16. );
    HDP.SetDefault( "FTF_MESON_NUCDESTR_P3_TGT", 2.1, 0., 4. );
    //
    HDP.SetDefault( "FTF_MESON_PT2_NUCDESTR_P1", 0.035, 0., 0.25 ); 
    HDP.SetDefault( "FTF_MESON_PT2_NUCDESTR_P2", 0.04, 0., 0.25 ); 
    HDP.SetDefault( "FTF_MESON_PT2_NUCDESTR_P3", 4.0, 2., 16. ); 
    HDP.SetDefault( "FTF_MESON_PT2_NUCDESTR_P4", 2.5, 0., 4. ); 
    //
    HDP.SetDefault( "FTF_MESON_NUCDESTR_R2", 1.5*CLHEP::fermi*CLHEP::fermi, 
                                              0.5*CLHEP::fermi*CLHEP::fermi, 
					      2.*CLHEP::fermi*CLHEP::fermi  );
    HDP.SetDefault( "FTF_MESON_EXCI_E_PER_WNDNUCLN", 40.*CLHEP::MeV, 0., 100.*CLHEP::MeV );
    HDP.SetDefault( "FTF_MESON_NUCDESTR_DISP", 0.3, 0.1, 0.4 );

  }
};


G4FTFSettingDefaultHDP FTFDefaultsHDP;  


G4FTFParamCollection::G4FTFParamCollection()
{

  Reset(); // zero out everything
   
  //
  // keep the 2 parameters below fixed for now (i.e. do not take them from HDP)
  //
  fNuclearProjDestructP2 = 4.0;
  fNuclearProjDestructP3 = 2.1;
  
}


void G4FTFParamCollection::Reset()
{
  // parameters of excitation

  // Proc=0 --> Qexchg w/o excitation
  fProc0A1 = 0.; 
  fProc0B1 = 0.;
  fProc0A2 = 0.;
  fProc0B2 = 0.;
  fProc0A3 = 0.;
  fProc0Atop = 0.;
  fProc0Ymin = 0.;

  // Proc=1 --> Qexchg w/excitation
  fProc1A1 = 0.;
  fProc1B1 = 0.;
  fProc1A2 = 0.;
  fProc1B2 = 0.;
  fProc1A3 = 0.;
  fProc1Atop = 0.;
  fProc1Ymin = 0.;

  // Proc=2 & Proc=3 for ( AbsProjectileBaryonNumber > 1  ||  NumberOfTargetNucleons > 1 )
  // Do NOT do anything as it's set once and for all !!!

  fProjDiffDissociation = false;
  fTgtDiffDissociation = false;

  // Proc=4 --> Qexchg w/additional multiplier in excitation  
  fProc4A1 = 0.;
  fProc4B1 = 0.;
  fProc4A2 = 0.;
  fProc4B2 = 0.; 
  fProc4A3 = 0.;
  fProc4Atop = 0.;
  fProc4Ymin = 0.;

  // parameters of participating baryon excitation

  fDeltaProbAtQuarkExchange = 0.; 
  fProbOfSameQuarkExchange = 0.;
  fProjMinDiffMass = 0.;
  fProjMinNonDiffMass = 0.;
  fTgtMinDiffMass = 0.;
  fTgtMinNonDiffMass = 0.;
  fAveragePt2 = 0.;
  fProbLogDistrPrD = 0.;
  fProbLogDistr = 0.;
   
  // parameters of nuclear distruction

  // COMMONs
  fNuclearProjDestructP1 = 0.;
  fNuclearTgtDestructP1 = 0.;
  fNuclearProjDestructP2 = 0.;
  fNuclearProjDestructP3 = 0.;
  fNuclearTgtDestructP2 = 0.;
  fNuclearTgtDestructP3 = 0.;
  fPt2NuclearDestructP1 = 0.;
  fPt2NuclearDestructP2 = 0.;
  fPt2NuclearDestructP3 = 0.;
  fPt2NuclearDestructP4 = 0.; 

  // baryons
  fR2ofNuclearDestruct = 0.;
  fExciEnergyPerWoundedNucleon = 0.;
  fDofNuclearDestruct = 0.;
  fMaxPt2ofNuclearDestruct = 0.;

  return;
}

//============================================================================

G4FTFParamCollBaryonProj::G4FTFParamCollBaryonProj()
   : G4FTFParamCollection()
{

  // parameters of participating hadron (baryon) excitation
  //
  // baryons projectile
  //
  // Proc=0 --> Qexchg w/o excitation
  //
  /* As of Oct. 31, 2017 keep these fixed
  HDP.DeveloperGet( "FTF_BARYON_PROC0_A1",   fProc0A1 );
  HDP.DeveloperGet( "FTF_BARYON_PROC0_B1",   fProc0B1 );
  HDP.DeveloperGet( "FTF_BARYON_PROC0_A2",   fProc0A2 ); 
  HDP.DeveloperGet( "FTF_BARYON_PROC0_B2",   fProc0B2 ); 
  HDP.DeveloperGet( "FTF_BARYON_PROC0_A3",   fProc0A3 );
  HDP.DeveloperGet( "FTF_BARYON_PROC0_ATOP", fProc0Atop ); 
  HDP.DeveloperGet( "FTF_BARYON_PROC0_YMIN", fProc0Ymin );
  */ 
  //
  fProc0A1 =  13.71; 
  fProc0B1 =   1.75;
  fProc0A2 = -30.69; // (or -214.5 as in Doc ?)
  fProc0B2 =   3.;   // ( or 4. as in Doc ?)
  fProc0A3 =   0.;
  fProc0Atop = 1.;   // ( or 0.5 as in Doc ?)
  fProc0Ymin = 0.93; // (or 1.1 as in Doc ?)
  //
  // Proc=1 --> Qexchg w/excitation
  //
  /* As of Oct. 31, 2017 keep these fixed
  HDP.DeveloperGet( "FTF_BARYON_PROC1_A1",   fProc1A1 );
  HDP.DeveloperGet( "FTF_BARYON_PROC1_B1",   fProc1B1 );
  HDP.DeveloperGet( "FTF_BARYON_PROC1_A2",   fProc1A2 ); 
  HDP.DeveloperGet( "FTF_BARYON_PROC1_B2",   fProc1B2 ); 
  HDP.DeveloperGet( "FTF_BARYON_PROC1_A3",   fProc1A3 );
  HDP.DeveloperGet( "FTF_BARYON_PROC1_ATOP", fProc1Atop ); 
  HDP.DeveloperGet( "FTF_BARYON_PROC1_YMIN", fProc1Ymin ); 
  */
  //
  fProc1A1 =  25.;
  fProc1B1 =   1.;
  fProc1A2 = -50.34;
  fProc1B2 =   1.5;
  fProc1A3 =   0.;
  fProc1Atop = 0.;
  fProc1Ymin = 1.4;
  //
  // Proc=2 & Proc=3 for the case ( AbsProjectileBaryonNumber > 10 ||  NumberOfTargetNucleons > 10 )
  // (diffraction dissociation)
  // NOTE-1: used to be ( AbsProjectileBaryonNumber > 1 ||  NumberOfTargetNucleons > 1 )...
  // NOTE-2: As of 10.5, both are set to false (via HDP)
  //
  HDP.DeveloperGet( "FTF_BARYON_DIFF_DISSO_PROJ", fProjDiffDissociation );
  HDP.DeveloperGet( "FTF_BARYON_DIFF_DISSO_TGT",  fTgtDiffDissociation );
  //
  //
  // Proc=4 --> Qexchg "w/additional multiplier" in excitation 
  //
  /* As of Oct. 31, 2017 keep these fixed
  HDP.DeveloperGet( "FTF_BARYON_PROC4_A1",   fProc4A1 );
  HDP.DeveloperGet( "FTF_BARYON_PROC4_B1",   fProc4B1 );
  HDP.DeveloperGet( "FTF_BARYON_PROC4_A2",   fProc4A2 ); 
  HDP.DeveloperGet( "FTF_BARYON_PROC4_B2",   fProc4B2 ); 
  HDP.DeveloperGet( "FTF_BARYON_PROC4_A3",   fProc4A3 );
  HDP.DeveloperGet( "FTF_BARYON_PROC4_ATOP", fProc4Atop ); 
  HDP.DeveloperGet( "FTF_BARYON_PROC4_YMIN", fProc4Ymin );
  */ 
  // 
  fProc4A1 =   0.6; // (or 1. as in Doc ?)
  fProc4B1 =   0.;
  fProc4A2 =  -1.2; // (or -2.01 as in Doc ?)
  fProc4B2 =   0.5; 
  fProc4A3 =   0.;
  fProc4Atop = 0.;
  fProc4Ymin = 1.4;
  //
  //
  HDP.DeveloperGet( "FTF_BARYON_DELTA_PROB_QEXCHG", fDeltaProbAtQuarkExchange );
  HDP.DeveloperGet( "FTF_BARYON_PROB_SAME_QEXCHG", fProbOfSameQuarkExchange );
  HDP.DeveloperGet( "FTF_BARYON_DIFF_M_PROJ", fProjMinDiffMass );
  HDP.DeveloperGet( "FTF_BARYON_NONDIFF_M_PROJ", fProjMinNonDiffMass );
  HDP.DeveloperGet( "FTF_BARYON_DIFF_M_TGT", fTgtMinDiffMass );
  HDP.DeveloperGet( "FTF_BARYON_NONDIFF_M_TGT", fTgtMinNonDiffMass );
  HDP.DeveloperGet( "FTF_BARYON_AVRG_PT2", fAveragePt2 );
  //
  // JVY - Per Alberto R., we're curretly keeping these two settings fixed,
  // thus they're defined here explicitly, rather than via HDP
  //
  // HDP.DeveloperGet( "FTF_BARYON_PROB_DISTR_PROJ", fProbLogDistrPrD );
  // HDP.DeveloperGet( "FTF_BARYON_PROB_DISTR_TGT", fProbLogDistr );
  fProbLogDistrPrD          = 0.55; 
  fProbLogDistr             = 0.55; 
   
  // nuclear destruction
  //
  // ---> LATER !!! ---> fBaryonMaxNumberOfCollisions = 2.;
  //

  HDP.DeveloperGet( "FTF_BARYON_NUCDESTR_P1_PROJ", fNuclearProjDestructP1 );
  HDP.DeveloperGet( "FTF_BARYON_NUCDESTR_P1_NBRN_PROJ",fNuclearProjDestructP1_NBRNDEP );
  //
  // keep the 2 parameters below fixed for now (i.e. do not take them from HDP)
  //
  fNuclearProjDestructP2 = 4.0;
  fNuclearProjDestructP3 = 2.1;

  HDP.DeveloperGet( "FTF_BARYON_NUCDESTR_P1_TGT", fNuclearTgtDestructP1 );
  HDP.DeveloperGet( "FTF_BARYON_NUCDESTR_P1_ADEP_TGT", fNuclearTgtDestructP1_ADEP );
  HDP.DeveloperGet( "FTF_BARYON_NUCDESTR_P2_TGT", fNuclearTgtDestructP2 );
  HDP.DeveloperGet( "FTF_BARYON_NUCDESTR_P3_TGT", fNuclearTgtDestructP3 );
  //
  HDP.DeveloperGet( "FTF_BARYON_PT2_NUCDESTR_P1", fPt2NuclearDestructP1 ); 
  HDP.DeveloperGet( "FTF_BARYON_PT2_NUCDESTR_P2", fPt2NuclearDestructP2 ); 
  HDP.DeveloperGet( "FTF_BARYON_PT2_NUCDESTR_P3", fPt2NuclearDestructP3 ); 
  HDP.DeveloperGet( "FTF_BARYON_PT2_NUCDESTR_P4", fPt2NuclearDestructP4 ); 
  //
  HDP.DeveloperGet( "FTF_BARYON_NUCDESTR_R2", fR2ofNuclearDestruct );
  HDP.DeveloperGet( "FTF_BARYON_EXCI_E_PER_WNDNUCLN", fExciEnergyPerWoundedNucleon );
  //
  HDP.DeveloperGet( "FTF_BARYON_NUCDESTR_DISP", fDofNuclearDestruct ); // NOTE: "Dof" means "Dispersion of..."
  // 
  // NOTE-1: this parameter has changed from 1. to 9. between 10.2 and 10.3.ref07 !!!
  //         ... then it went back to 1. for the 10.4-candidate... 
  // NOTE-2: this is a "technical" parameter, it should not be changed; this is why
  //         it is defined explicitly rather than via HDP
  // --> HDP.DeveloperGet( "FTF_BARYON_NUCDESTR_MAXPT2", fMaxPt2ofNuclearDestruct );
  fMaxPt2ofNuclearDestruct     = 9. * CLHEP::GeV*CLHEP::GeV; 
}


G4FTFParamCollMesonProj::G4FTFParamCollMesonProj()
   : G4FTFParamCollection()
{
  // nuclear destruction
  //
  // JVY, Aug.8, 2018 --> Feb.14, 2018 --> June 25, 2019: 
  // These parameters are common for all mesons
  //

  HDP.DeveloperGet( "FTF_MESON_NUCDESTR_P1_TGT", fNuclearTgtDestructP1 );
  HDP.DeveloperGet( "FTF_MESON_NUCDESTR_P1_ADEP_TGT", fNuclearTgtDestructP1_ADEP );
  HDP.DeveloperGet( "FTF_MESON_NUCDESTR_P2_TGT", fNuclearTgtDestructP2 );
  HDP.DeveloperGet( "FTF_MESON_NUCDESTR_P3_TGT", fNuclearTgtDestructP3 );
  //
  HDP.DeveloperGet( "FTF_MESON_PT2_NUCDESTR_P1", fPt2NuclearDestructP1 ); 
  HDP.DeveloperGet( "FTF_MESON_PT2_NUCDESTR_P2", fPt2NuclearDestructP2 ); 
  HDP.DeveloperGet( "FTF_MESON_PT2_NUCDESTR_P3", fPt2NuclearDestructP3 ); 
  HDP.DeveloperGet( "FTF_MESON_PT2_NUCDESTR_P4", fPt2NuclearDestructP4 ); 
  //
  HDP.DeveloperGet( "FTF_MESON_NUCDESTR_R2", fR2ofNuclearDestruct );
  HDP.DeveloperGet( "FTF_MESON_EXCI_E_PER_WNDNUCLN", fExciEnergyPerWoundedNucleon );
  HDP.DeveloperGet( "FTF_MESON_NUCDESTR_DISP", fDofNuclearDestruct ); // NOTE: "Dof" means "Dispersion of..." 
  //
  // NOTE: it is a "technical" parameter, it should not be changed; 
  //       this is why it is defined explicitly rather than via HDP
  //
  fMaxPt2ofNuclearDestruct     = 1. * CLHEP::GeV*CLHEP::GeV; 
}


G4FTFParamCollPionProj::G4FTFParamCollPionProj()
   : G4FTFParamCollMesonProj()
{
  // parameters of participating pion excitation (pi+/- or pi0)
  //
  // Proc=0 --> Qexchg w/o excitation
  //
  /* As of June 25, 2019 keep these fixed
  HDP.DeveloperGet( "FTF_PION_PROC0_A1",   fProc0A1 );
  HDP.DeveloperGet( "FTF_PION_PROC0_B1",   fProc0B1 );
  HDP.DeveloperGet( "FTF_PION_PROC0_A2",   fProc0A2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC0_B2",   fProc0B2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC0_A3",   fProc0A3 );
  HDP.DeveloperGet( "FTF_PION_PROC0_ATOP", fProc0Atop ); 
  HDP.DeveloperGet( "FTF_PION_PROC0_YMIN", fProc0Ymin );
  */ 
  //
  fProc0A1 = 150.0; 
  fProc0B1 =   1.8;
  fProc0A2 =-247.3; 
  fProc0B2 =   2.3;   // ( or 4. as in Doc ?)
  fProc0A3 =   0.;
  fProc0Atop = 1.;   // ( or 0.5 as in Doc ?)
  fProc0Ymin = 2.3; // (or 1.1 as in Doc ?)
  //
  // Proc=1 --> Qexchg w/excitation
  //
  /* As of Oct. 31, 2017 keep these fixed
  HDP.DeveloperGet( "FTF_PION_PROC1_A1",   fProc1A1 );
  HDP.DeveloperGet( "FTF_PION_PROC1_B1",   fProc1B1 );
  HDP.DeveloperGet( "FTF_PION_PROC1_A2",   fProc1A2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC1_B2",   fProc1B2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC1_A3",   fProc1A3 );
  HDP.DeveloperGet( "FTF_PION_PROC1_ATOP", fProc1Atop ); 
  HDP.DeveloperGet( "FTF_PION_PROC1_YMIN", fProc1Ymin ); 
  */
  //
  fProc1A1 =   5.77;
  fProc1B1 =   0.6;
  fProc1A2 =  -5.77;
  fProc1B2 =   0.8;
  fProc1A3 =   0.;
  fProc1Atop = 0.;
  fProc1Ymin = 0.0;
  //
  // Proc=2 --> Projectile diffraction
  //
  /* As of Oct. 31, 2017 keep these fixed
  HDP.DeveloperGet( "FTF_PION_PROC2_A1",   fProc2A1 );
  HDP.DeveloperGet( "FTF_PION_PROC2_B1",   fProc2B1 );
  HDP.DeveloperGet( "FTF_PION_PROC2_A2",   fProc2A2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC2_B2",   fProc2B2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC2_A3",   fProc2A3 );
  HDP.DeveloperGet( "FTF_PION_PROC2_ATOP", fProc2Atop ); 
  HDP.DeveloperGet( "FTF_PION_PROC2_YMIN", fProc2Ymin ); 
  */
  //
  fProc2A1 =   2.27;
  fProc2B1 =   0.5;
  fProc2A2 =-98052.0;
  fProc2B2 =   4.0;
  fProc2A3 =   0.;
  fProc2Atop = 0.;
  fProc2Ymin = 3.0;
  //
  // Proc=3 --> Target diffraction
  //
  /* As of Oct. 31, 2017 keep these fixed
  HDP.DeveloperGet( "FTF_PION_PROC3_A1",   fProc3A1 );
  HDP.DeveloperGet( "FTF_PION_PROC3_B1",   fProc3B1 );
  HDP.DeveloperGet( "FTF_PION_PROC3_A2",   fProc3A2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC3_B2",   fProc3B2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC3_A3",   fProc3A3 );
  HDP.DeveloperGet( "FTF_PION_PROC3_ATOP", fProc3Atop ); 
  HDP.DeveloperGet( "FTF_PION_PROC3_YMIN", fProc3Ymin ); 
  */
  //
  fProc3A1 =   7.0;
  fProc3B1 =   0.9;
  fProc3A2 = -85.28;
  fProc3B2 =   1.9;
  fProc3A3 =   0.08;
  fProc3Atop = 0.;
  fProc3Ymin = 2.2;
  //
  // for Proc2 & Proc3, pprojectile or target diffraction can be turned ON/OFF
  // if num.baryons >10 (which is strange for projectile which is pion !!!)
  //
  HDP.DeveloperGet( "FTF_PION_DIFF_DISSO_PROJ", fProjDiffDissociation );
  HDP.DeveloperGet( "FTF_PION_DIFF_DISSO_TGT",  fTgtDiffDissociation );
  //
  // Proc=4 --> Qexchg "w/additional multiplier" in excitation 
  //
  /* As of Oct. 31, 2017 keep these fixed
  HDP.DeveloperGet( "FTF_PION_PROC4_A1",   fProc4A1 );
  HDP.DeveloperGet( "FTF_PION_PROC4_B1",   fProc4B1 );
  HDP.DeveloperGet( "FTF_PION_PROC4_A2",   fProc4A2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC4_B2",   fProc4B2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC4_A3",   fProc4A3 );
  HDP.DeveloperGet( "FTF_PION_PROC4_ATOP", fProc4Atop ); 
  HDP.DeveloperGet( "FTF_PION_PROC4_YMIN", fProc4Ymin );
  */ 
  // 
  fProc4A1 =   1.0; 
  fProc4B1 =   0.;
  fProc4A2 = -11.02; 
  fProc4B2 =   1.0; 
  fProc4A3 =   0.;
  fProc4Atop = 0.;
  fProc4Ymin = 2.4;
  //
  //
  HDP.DeveloperGet( "FTF_PION_DELTA_PROB_QEXCHG", fDeltaProbAtQuarkExchange );
  HDP.DeveloperGet( "FTF_PION_DIFF_M_PROJ", fProjMinDiffMass );
  HDP.DeveloperGet( "FTF_PION_NONDIFF_M_PROJ", fProjMinNonDiffMass );
  HDP.DeveloperGet( "FTF_PION_DIFF_M_TGT", fTgtMinDiffMass );
  HDP.DeveloperGet( "FTF_PION_NONDIFF_M_TGT", fTgtMinNonDiffMass );
  HDP.DeveloperGet( "FTF_PION_AVRG_PT2", fAveragePt2 );
  //
  fProbOfSameQuarkExchange  = 0.; // This does NOT seem to apply to the pion case 
  //
  // currently keep these two parameters fixed
  // thus they're defined here explicitly, rather than via HDP
  //
  fProbLogDistrPrD          = 0.55; 
  fProbLogDistr             = 0.55; 
}


//============================================================================

G4FTFParameters::~G4FTFParameters() {
  if ( StringMass ) delete StringMass;
}

//============================================================================

void G4FTFParameters::Reset()
{
  FTFhNcmsEnergy  = 0.0; 
  FTFXtotal = 0.0;
  FTFXelastic = 0.0;
  FTFXinelastic = 0.0; 
  FTFXannihilation = 0.0;
  ProbabilityOfAnnihilation = 0.0; 
  ProbabilityOfElasticScatt  = 0.0;
  RadiusOfHNinteractions2 = 0.0; 
  FTFSlope = 0.0;  
  AvaragePt2ofElasticScattering = 0.0; 
  FTFGamma0 = 0.0;
  DeltaProbAtQuarkExchange = 0.0; 
  ProbOfSameQuarkExchange = 0.0; 
  ProjMinDiffMass = 0.0; 
  ProjMinNonDiffMass = 0.0; 
  ProbLogDistrPrD = 0.0;
  TarMinDiffMass = 0.0; 
  TarMinNonDiffMass = 0.0;
  AveragePt2 = 0.0; 
  ProbLogDistr  = 0.0;
  Pt2kink = 0.0;
  MaxNumberOfCollisions = 0.0;  
  ProbOfInelInteraction = 0.0; 
  CofNuclearDestructionPr = 0.0; 
  CofNuclearDestruction = 0.0;
  R2ofNuclearDestruction  = 0.0; 
  ExcitationEnergyPerWoundedNucleon  = 0.0;
  DofNuclearDestruction  = 0.0; 
  Pt2ofNuclearDestruction = 0.0; 
  MaxPt2ofNuclearDestruction = 0.0; 

  for ( G4int i = 0; i < 4; i++ ) {
    for ( G4int j = 0; j < 7; j++ ) {
      ProcParams[i][j] = 0.0;
    }
  }

  return;
}

