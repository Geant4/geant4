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
#include "G4HadronicParameters.hh"

//============================================================================

//#define debugFTFparams

//============================================================================

G4FTFParameters::G4FTFParameters() 
{
  // Set-up alternative sets of FTF parameters (called "tunes").
  // Note that the very first tune (with indexTune == 0) corresponds to the default
  // set of parameters, which does not need to be set-up explicitly: that's why
  // the for loop below starts from 1 and not from 0.
  // The check whether an alternative tune has been switched on is done at the
  // level of the G4FTFParamCollection::SetTune method.
  for ( G4int indexTune = 1; indexTune < G4FTFTunings::sNumberOfTunes; ++indexTune ) {
    fArrayParCollBaryonProj[indexTune].SetTune(indexTune);
    fArrayParCollMesonProj[indexTune].SetTune(indexTune);
    fArrayParCollPionProj[indexTune].SetTune(indexTune);
  }
  
  StringMass = new G4LundStringFragmentation;  // for estimation of min. mass of diffr. states
  Reset();
  csGGinstance = 
    G4CrossSectionDataSetRegistry::Instance()->GetComponentCrossSection("Glauber-Gribov");
  if (!csGGinstance) {
    csGGinstance = new G4ComponentGGHadronNucleusXsc();
  }

  EnableDiffDissociationForBGreater10 = G4HadronicParameters::Instance()->EnableDiffDissociationForBGreater10();
  
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
    
    const G4int indexTune = G4FTFTunings::Instance()->GetIndexTune( particle, KineticEnergy );
    
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
    SetParams( 0, fArrayParCollBaryonProj[indexTune].GetProc0A1(),
	          fArrayParCollBaryonProj[indexTune].GetProc0B1(),          
		  fArrayParCollBaryonProj[indexTune].GetProc0A2(),
	          fArrayParCollBaryonProj[indexTune].GetProc0B2(), 
		  fArrayParCollBaryonProj[indexTune].GetProc0A3(),
	          fArrayParCollBaryonProj[indexTune].GetProc0Atop(),     
		  fArrayParCollBaryonProj[indexTune].GetProc0Ymin() );   // Qexchange without Exc.
    SetParams( 1, fArrayParCollBaryonProj[indexTune].GetProc1A1(),
	          fArrayParCollBaryonProj[indexTune].GetProc1B1(),          
		  fArrayParCollBaryonProj[indexTune].GetProc1A2(),
	          fArrayParCollBaryonProj[indexTune].GetProc1B2(), 
		  fArrayParCollBaryonProj[indexTune].GetProc1A3(),
	          fArrayParCollBaryonProj[indexTune].GetProc1Atop(),     
		  fArrayParCollBaryonProj[indexTune].GetProc1Ymin() );   // Qexchange with Exc.
    if ( Xinel > 0.0 ) {
      SetParams( 2, 6.0/Xinel, 0.0, -6.0/Xinel*16.28, 3.0, 0.0, 0.0, 0.93 );  // Projectile diffraction
      SetParams( 3, 6.0/Xinel, 0.0, -6.0/Xinel*16.28, 3.0, 0.0, 0.0, 0.93 );  // Target diffraction

      SetParams( 4, fArrayParCollBaryonProj[indexTune].GetProc4A1(),
		    fArrayParCollBaryonProj[indexTune].GetProc4B1(),          
		    fArrayParCollBaryonProj[indexTune].GetProc4A2(),
		    fArrayParCollBaryonProj[indexTune].GetProc4B2(), 
		    fArrayParCollBaryonProj[indexTune].GetProc4A3(),
		    fArrayParCollBaryonProj[indexTune].GetProc4Atop(),     
		    fArrayParCollBaryonProj[indexTune].GetProc4Ymin() );  // Qexchange with Exc. Additional multiplier
    } else {  // if Xinel=0., zero everything out (obviously)
      SetParams( 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
      SetParams( 3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
      SetParams( 4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
    }

    if ( (AbsProjectileBaryonNumber > 10 || NumberOfTargetNucleons > 10) && !EnableDiffDissociationForBGreater10 ) {
      // It is not decided what to do with diffraction dissociation in Had-Nucl and Nucl-Nucl interactions
      // For the moment both ProjDiffDisso & TgtDiffDisso for A > 10 are set to false,
      // so both projectile and target diffraction are turned OFF
      if ( ! fArrayParCollBaryonProj[indexTune].IsProjDiffDissociation() )
         SetParams( 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -100.0 );  // Projectile diffraction
      if ( ! fArrayParCollBaryonProj[indexTune].IsTgtDiffDissociation() )
         SetParams( 3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -100.0 );  // Target diffraction
    }

    SetDeltaProbAtQuarkExchange( fArrayParCollBaryonProj[indexTune].GetDeltaProbAtQuarkExchange() );

    if ( NumberOfTargetNucleons > 26 ) {
      SetProbOfSameQuarkExchange( 1.0 );
    } else {
      SetProbOfSameQuarkExchange( fArrayParCollBaryonProj[indexTune].GetProbOfSameQuarkExchange() );
    }

    SetProjMinDiffMass(    fArrayParCollBaryonProj[indexTune].GetProjMinDiffMass() );     // GeV 
    SetProjMinNonDiffMass( fArrayParCollBaryonProj[indexTune].GetProjMinNonDiffMass() );  // GeV 

    SetTarMinDiffMass(     fArrayParCollBaryonProj[indexTune].GetTgtMinDiffMass() );      // GeV
    SetTarMinNonDiffMass(  fArrayParCollBaryonProj[indexTune].GetTgtMinNonDiffMass() );   // GeV 

    SetAveragePt2(         fArrayParCollBaryonProj[indexTune].GetAveragePt2() );          // GeV^2       
    SetProbLogDistrPrD(    fArrayParCollBaryonProj[indexTune].GetProbLogDistrPrD() ); 
    SetProbLogDistr(       fArrayParCollBaryonProj[indexTune].GetProbLogDistr() ); 

  } else if ( ProjectilePDGcode == -2212  ||  ProjectilePDGcode == -2112 ) {  // Projectile is anti_proton or anti_neutron
    
    // Below, in the call to the G4FTFTunings::GetIndexTune method, we pass the proton
    // as projectile, instead of the real one, because for switching on/off diffraction
    // we assume the same treatment for anti_proton/anti_neutron as for proton/neutron,
    // whereas all other parameters for anti_proton/anti_neutron are hardwired.
    const G4int indexTune = G4FTFTunings::Instance()->GetIndexTune( G4Proton::Definition(), KineticEnergy );
    
    //        Proc#   A1      B1            A2       B2   A3   Atop       Ymin
    SetParams( 0,      0.0 , 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,  1000.0  );  // Qexchange without Exc. 
    SetParams( 1,      0.0 , 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,  1000.0  );  // Qexchange with    Exc.
    if ( Xinel > 0.) {
      SetParams( 2, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0 , 0.93 );  // Projectile diffraction
      SetParams( 3, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0 , 0.93 );  // Target diffraction
      SetParams( 4,       1.0, 0.0 ,             0.0, 0.0 , 0.0, 0.0 , 0.93 );  // Qexchange with    Exc. Additional multiply
    } else {
      SetParams( 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
      SetParams( 3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
      SetParams( 4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
    }

    if ( AbsProjectileBaryonNumber > 10  ||  NumberOfTargetNucleons > 10 ) {
      // It is not decided what to do with diffraction dissociation in Had-Nucl and Nucl-Nucl interactions
      // For the moment both ProjDiffDisso & TgtDiffDisso are set to false,
      // so both projectile and target diffraction are turned OFF
      if ( ! fArrayParCollBaryonProj[indexTune].IsProjDiffDissociation() )
         SetParams( 2,       0.0, 0.0 ,           0.0  , 0.0 , 0.0, 0.0   , -100.0  );  // Projectile diffraction
      if ( ! fArrayParCollBaryonProj[indexTune].IsTgtDiffDissociation() )
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
    
    const G4int indexTune = G4FTFTunings::Instance()->GetIndexTune( particle, KineticEnergy );
 
    //        Proc#   A1      B1            A2       B2      A3   Atop        Ymin
    /* --> original code
    SetParams( 0,  150.0,    1.8 ,       -247.3,     2.3,    0.,   1. ,       2.3 ); 
    SetParams( 1,   5.77,    0.6 ,        -5.77,     0.8,    0.,   0. ,       0.0 );
    SetParams( 2,   2.27,    0.5 ,     -98052.0,     4.0,    0.,   0. ,       3.0 ); 
    SetParams( 3,    7.0,    0.9,        -85.28,     1.9,  0.08,   0. ,       2.2 );
    SetParams( 4,    1.0,    0.0 ,       -11.02,     1.0,   0.0,   0. ,       2.4 );  // Qexchange with    Exc. Additional multiply
    */
    //         Proc#
    SetParams( 0, fArrayParCollPionProj[indexTune].GetProc0A1(),
	          fArrayParCollPionProj[indexTune].GetProc0B1(),          
		  fArrayParCollPionProj[indexTune].GetProc0A2(),
	          fArrayParCollPionProj[indexTune].GetProc0B2(), 
		  fArrayParCollPionProj[indexTune].GetProc0A3(),
	          fArrayParCollPionProj[indexTune].GetProc0Atop(),     
		  fArrayParCollPionProj[indexTune].GetProc0Ymin() );   // Qexchange without Exc.
    SetParams( 1, fArrayParCollPionProj[indexTune].GetProc1A1(),
	          fArrayParCollPionProj[indexTune].GetProc1B1(),          
		  fArrayParCollPionProj[indexTune].GetProc1A2(),
	          fArrayParCollPionProj[indexTune].GetProc1B2(), 
		  fArrayParCollPionProj[indexTune].GetProc1A3(),
	          fArrayParCollPionProj[indexTune].GetProc1Atop(),     
		  fArrayParCollPionProj[indexTune].GetProc1Ymin() );   // Qexchange with Exc.
    SetParams( 2, fArrayParCollPionProj[indexTune].GetProc2A1(),
	          fArrayParCollPionProj[indexTune].GetProc2B1(),          
		  fArrayParCollPionProj[indexTune].GetProc2A2(),
	          fArrayParCollPionProj[indexTune].GetProc2B2(), 
		  fArrayParCollPionProj[indexTune].GetProc2A3(),
	          fArrayParCollPionProj[indexTune].GetProc2Atop(),     
		  fArrayParCollPionProj[indexTune].GetProc2Ymin() );   // Projectile diffraction
    SetParams( 3, fArrayParCollPionProj[indexTune].GetProc3A1(),
	          fArrayParCollPionProj[indexTune].GetProc3B1(),          
		  fArrayParCollPionProj[indexTune].GetProc3A2(),
	          fArrayParCollPionProj[indexTune].GetProc3B2(), 
		  fArrayParCollPionProj[indexTune].GetProc3A3(),
	          fArrayParCollPionProj[indexTune].GetProc3Atop(),     
		  fArrayParCollPionProj[indexTune].GetProc3Ymin() );   // Target diffraction
    SetParams( 4, fArrayParCollPionProj[indexTune].GetProc4A1(),
	          fArrayParCollPionProj[indexTune].GetProc4B1(),          
		  fArrayParCollPionProj[indexTune].GetProc4A2(),
	          fArrayParCollPionProj[indexTune].GetProc4B2(), 
		  fArrayParCollPionProj[indexTune].GetProc4A3(),
	          fArrayParCollPionProj[indexTune].GetProc4Atop(),     
		  fArrayParCollPionProj[indexTune].GetProc4Ymin() );   // Qexchange with Exc. Additional multiply

    // NOTE: how can it be |ProjectileBaryonNumber| > 10 if projectile is a pion ??? 
    //
    if ( AbsProjectileBaryonNumber > 10  ||  NumberOfTargetNucleons > 10 ) {
       if ( ! fArrayParCollPionProj[indexTune].IsProjDiffDissociation() )
          SetParams( 2,      0.0 , 0.0 ,           0.0  , 0.0 , 0.0, 0.0  ,  -100.0 );  // Projectile diffraction
       if ( ! fArrayParCollPionProj[indexTune].IsTgtDiffDissociation() )
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
    SetDeltaProbAtQuarkExchange( fArrayParCollPionProj[indexTune].GetDeltaProbAtQuarkExchange() );
    SetProjMinDiffMass( fArrayParCollPionProj[indexTune].GetProjMinDiffMass() );                  // GeV 
    SetProjMinNonDiffMass( fArrayParCollPionProj[indexTune].GetProjMinNonDiffMass() );            // GeV 
    SetTarMinDiffMass( fArrayParCollPionProj[indexTune].GetTgtMinDiffMass() );                    // GeV
    SetTarMinNonDiffMass( fArrayParCollPionProj[indexTune].GetTgtMinNonDiffMass() );              // GeV 
    SetAveragePt2( fArrayParCollPionProj[indexTune].GetAveragePt2() );                            // GeV^2       
    SetProbLogDistrPrD( fArrayParCollPionProj[indexTune].GetProbLogDistrPrD() ); 
    SetProbLogDistr( fArrayParCollPionProj[indexTune].GetProbLogDistr() ); 
    
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

    const G4int indexTune = G4FTFTunings::Instance()->GetIndexTune( particle, KineticEnergy );

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
    double coeff = fArrayParCollMesonProj[indexTune].GetNuclearTgtDestructP1();
    // 
    // NOTE (JVY): Set this switch to false/true on line 138
    //
    if ( fArrayParCollMesonProj[indexTune].IsNuclearTgtDestructP1_ADEP() )       
    {                                                             
      coeff *= G4double(NumberOfTargetNucleons);                  
    }                                                             
    double exfactor = G4Exp( fArrayParCollMesonProj[indexTune].GetNuclearTgtDestructP2()
                           * (Ylab-fArrayParCollMesonProj[indexTune].GetNuclearTgtDestructP3()) );
    coeff *= exfactor;
    coeff /= ( 1.+ exfactor );

    SetCofNuclearDestruction( coeff );

    SetR2ofNuclearDestruction( fArrayParCollMesonProj[indexTune].GetR2ofNuclearDestruct() );
    SetDofNuclearDestruction( fArrayParCollMesonProj[indexTune].GetDofNuclearDestruct() );
    coeff = fArrayParCollMesonProj[indexTune].GetPt2NuclearDestructP2();
    exfactor = G4Exp( fArrayParCollMesonProj[indexTune].GetPt2NuclearDestructP3()
                    * (Ylab-fArrayParCollMesonProj[indexTune].GetPt2NuclearDestructP4()) );
    coeff *= exfactor;
    coeff /= ( 1. + exfactor );
    SetPt2ofNuclearDestruction( (fArrayParCollMesonProj[indexTune].GetPt2NuclearDestructP1()+coeff)*CLHEP::GeV*CLHEP::GeV ); 

    SetMaxPt2ofNuclearDestruction( fArrayParCollMesonProj[indexTune].GetMaxPt2ofNuclearDestruct() );
    SetExcitationEnergyPerWoundedNucleon( fArrayParCollMesonProj[indexTune].GetExciEnergyPerWoundedNucleon() ); 

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

    // Below, in the call to the G4FTFTunings::GetIndexTune method, we pass the proton
    // as projectile, instead of the real one, because for the treatment of nuclear
    // destruction, we assume for this category of hadron projectiles the same treatment
    // as for "baryon".
    const G4int indexTune = G4FTFTunings::Instance()->GetIndexTune( G4Proton::Definition(), KineticEnergy );
    
    // NOTE (JVY) FIXME !!! Will decide later how/if to make this one configurable...
    //
    SetMaxNumberOfCollisions( Plab, 2.0 );

    // projectile destruction - does NOT really matter for particle projectile, only for a nucleus projectile
    //
    double coeff = 0.;
    coeff = fArrayParCollBaryonProj[indexTune].GetNuclearProjDestructP1();
    // 
    // NOTE (JVY): Set this switch to false/true on line 136
    //
    if ( fArrayParCollBaryonProj[indexTune].IsNuclearProjDestructP1_NBRNDEP() )   
    {                                                             
      coeff *= G4double(AbsProjectileBaryonNumber);               
    }                                                             
    double exfactor = G4Exp( fArrayParCollBaryonProj[indexTune].GetNuclearProjDestructP2()*
			     (Ylab-fArrayParCollBaryonProj[indexTune].GetNuclearProjDestructP3()) ); 
    coeff *= exfactor;
    coeff /= ( 1.+ exfactor );
    SetCofNuclearDestructionPr( coeff );

    // target desctruction
    //
    coeff = fArrayParCollBaryonProj[indexTune].GetNuclearTgtDestructP1();
    // 
    // NOTE (JVY): Set this switch to false/true on line 138
    //
    if ( fArrayParCollBaryonProj[indexTune].IsNuclearTgtDestructP1_ADEP() )       
    {                                                             
      coeff *= G4double(NumberOfTargetNucleons);                  
    }                                                             
    exfactor = G4Exp( fArrayParCollBaryonProj[indexTune].GetNuclearTgtDestructP2()*
		      (Ylab-fArrayParCollBaryonProj[indexTune].GetNuclearTgtDestructP3()) );
    coeff *= exfactor;
    coeff /= ( 1.+ exfactor );
    SetCofNuclearDestruction(   coeff );

    SetR2ofNuclearDestruction( fArrayParCollBaryonProj[indexTune].GetR2ofNuclearDestruct() );
    SetDofNuclearDestruction( fArrayParCollBaryonProj[indexTune].GetDofNuclearDestruct() );

    coeff = fArrayParCollBaryonProj[indexTune].GetPt2NuclearDestructP2();
    exfactor = G4Exp( fArrayParCollBaryonProj[indexTune].GetPt2NuclearDestructP3()*
		      (Ylab-fArrayParCollBaryonProj[indexTune].GetPt2NuclearDestructP4())  );
    coeff *= exfactor;
    coeff /= ( 1. + exfactor );
    SetPt2ofNuclearDestruction( (fArrayParCollBaryonProj[indexTune].GetPt2NuclearDestructP1()+coeff)*CLHEP::GeV*CLHEP::GeV ); 

    SetMaxPt2ofNuclearDestruction( fArrayParCollBaryonProj[indexTune].GetMaxPt2ofNuclearDestruct() );
    SetExcitationEnergyPerWoundedNucleon( fArrayParCollBaryonProj[indexTune].GetExciEnergyPerWoundedNucleon() ); 

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

