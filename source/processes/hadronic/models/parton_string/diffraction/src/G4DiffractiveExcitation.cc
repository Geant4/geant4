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
//      ---------------- G4DiffractiveExcitation --------------
//             by Gunter Folger, October 1998.
//        diffractive Excitation used by strings models
//             Take a projectile and a target
//             excite the projectile and target
//  Essential changed by V. Uzhinsky in November - December 2006
//  in order to put it in a correspondence with original FRITIOF
//  model. Variant of FRITIOF with nucleon de-excitation is implemented.
//  Other changes by V.Uzhinsky in May 2007 were introduced to fit
//  meson-nucleon interactions. Additional changes by V. Uzhinsky
//  were introduced in December 2006. They treat diffraction dissociation
//  processes more exactly.
//  Correct treatment of the diffraction dissociation - 2012, V. Uzhinsky
//  Mass distributions for resonances and uu-diquark suppression in protons,
//  and dd-diquarks suppression in neutrons were introduced by V. Uzhinsky, 2014
// ---------------------------------------------------------------------

#include "globals.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4DiffractiveExcitation.hh"
#include "G4FTFParameters.hh"
#include "G4ElasticHNScattering.hh"

#include "G4RotationMatrix.hh"
#include "G4ParticleDefinition.hh" 
#include "G4ParticleTable.hh"
#include "G4SampleResonance.hh"
#include "G4VSplitableHadron.hh"
#include "G4ExcitedString.hh"
#include "G4Neutron.hh"

#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

//#include "G4ios.hh"

//============================================================================

//#define debugFTFexcitation
//#define debug_heavyHadrons

//============================================================================

G4DiffractiveExcitation::G4DiffractiveExcitation() {}


//============================================================================

G4DiffractiveExcitation::~G4DiffractiveExcitation() {}


//============================================================================

G4bool G4DiffractiveExcitation::ExciteParticipants( G4VSplitableHadron*    projectile, 
                                                    G4VSplitableHadron*    target,
                                                    G4FTFParameters*       theParameters,
                                                    G4ElasticHNScattering* theElastic ) const {

  #ifdef debugFTFexcitation
  G4cout << G4endl << "FTF ExciteParticipants --------------" << G4endl;
  #endif

  CommonVariables common;

  // Projectile parameters
  common.Pprojectile = projectile->Get4Momentum();
  if ( common.Pprojectile.z() < 0.0 ) return false;
  common.ProjectilePDGcode = projectile->GetDefinition()->GetPDGEncoding();
  common.absProjectilePDGcode = std::abs( common.ProjectilePDGcode );
  common.M0projectile = projectile->GetDefinition()->GetPDGMass();  //Uzhi Aug.2019  common.Pprojectile.mag();
  G4double ProjectileRapidity = common.Pprojectile.rapidity();

  // Target parameters
  common.Ptarget = target->Get4Momentum();
  common.TargetPDGcode = target->GetDefinition()->GetPDGEncoding();
  common.absTargetPDGcode = std::abs( common.TargetPDGcode );
  common.M0target = target->GetDefinition()->GetPDGMass();  //Uzhi Aug.2019  common.Ptarget.mag();
  G4double TargetRapidity = common.Ptarget.rapidity();

  // Kinematical properties of the interactions
  G4LorentzVector Psum = common.Pprojectile + common.Ptarget;  // 4-momentum in Lab.
  common.S = Psum.mag2();
  common.SqrtS = std::sqrt( common.S ); 

  // Check off-shellness of the participants
  G4bool toBePutOnMassShell = true;  //Uzhi Aug.2019  false;
  common.MminProjectile = common.BrW.GetMinimumMass( projectile->GetDefinition() );
  /* Uzhi Aug.2019
  if ( common.M0projectile < common.MminProjectile ) {
    toBePutOnMassShell = true;
    common.M0projectile = common.BrW.SampleMass( projectile->GetDefinition(), 
                                                 projectile->GetDefinition()->GetPDGMass() 
                                                 + 5.0*projectile->GetDefinition()->GetPDGWidth() );
  }
  */
  common.M0projectile2 = common.M0projectile * common.M0projectile;
  common.ProjectileDiffStateMinMass    = theParameters->GetProjMinDiffMass();
  common.ProjectileNonDiffStateMinMass = theParameters->GetProjMinNonDiffMass();
  if ( common.M0projectile > common.ProjectileDiffStateMinMass ) {
    common.ProjectileDiffStateMinMass    = common.MminProjectile + 220.0*MeV;  //Uzhi Aug.2019  common.M0projectile + 220.0*MeV;
    common.ProjectileNonDiffStateMinMass = common.MminProjectile + 220.0*MeV;  //Uzhi Aug.2019  common.M0projectile + 220.0*MeV;
    if ( common.absProjectilePDGcode > 3000 ) {  // Strange baryon
      common.ProjectileDiffStateMinMass    += 140.0*MeV;
      common.ProjectileNonDiffStateMinMass += 140.0*MeV;
    }
  }
  common.MminTarget = common.BrW.GetMinimumMass( target->GetDefinition() );
  /* Uzhi Aug.2019
  if ( common.M0target < common.MminTarget ) {
    toBePutOnMassShell = true;
    common.M0target = common.BrW.SampleMass( target->GetDefinition(),
                                             target->GetDefinition()->GetPDGMass() 
                                             + 5.0*target->GetDefinition()->GetPDGWidth() );
  }
  */
  common.M0target2 = common.M0target * common.M0target;
  common.TargetDiffStateMinMass    = theParameters->GetTarMinDiffMass();
  common.TargetNonDiffStateMinMass = theParameters->GetTarMinNonDiffMass();
  if ( common.M0target > common.TargetDiffStateMinMass ) {
    common.TargetDiffStateMinMass    = common.MminTarget + 220.0*MeV;  //Uzhi Aug.2019  common.M0target + 220.0*MeV;
    common.TargetNonDiffStateMinMass = common.MminTarget + 220.0*MeV;  //Uzhi Aug.2019  common.M0target + 220.0*MeV;
    if ( common.absTargetPDGcode > 3000 ) {  // Strange baryon
      common.TargetDiffStateMinMass    += 140.0*MeV;
      common.TargetNonDiffStateMinMass += 140.0*MeV;
    }
  };
  #ifdef debugFTFexcitation
  G4cout << "Proj Targ PDGcodes " << common.ProjectilePDGcode << " " << common.TargetPDGcode << G4endl
         << "Mprojectile  Y " << common.Pprojectile.mag() << " " << ProjectileRapidity << G4endl        // Uzhi Aug.2019
         << "M0projectile Y " << common.M0projectile << " " << ProjectileRapidity << G4endl;
  G4cout << "Mtarget      Y " << common.Ptarget.mag() << " " << TargetRapidity << G4endl                // Uzhi Aug.2019
         << "M0target     Y " << common.M0target << " " << TargetRapidity << G4endl;
  G4cout << "Pproj " << common.Pprojectile << G4endl << "Ptarget " << common.Ptarget << G4endl;
  #endif

  // Transform momenta to cms and then rotate parallel to z axis;
  common.toCms = G4LorentzRotation( -1 * Psum.boostVector() );
  G4LorentzVector Ptmp = common.toCms * common.Pprojectile;
  if ( Ptmp.pz() <= 0.0 ) return false;  // "String" moving backwards in  CMS, abort collision!
  common.toCms.rotateZ( -1*Ptmp.phi() );
  common.toCms.rotateY( -1*Ptmp.theta() );
  common.toLab = common.toCms.inverse();
  common.Pprojectile.transform( common.toCms );
  common.Ptarget.transform( common.toCms );

  G4double SumMasses = common.M0projectile + common.M0target;  // + 220.0*MeV;
  #ifdef debugFTFexcitation
  G4cout << "SqrtS     " << common.SqrtS << G4endl << "M0pr M0tr SumM " << common.M0projectile 
         << " " << common.M0target << " " << SumMasses << G4endl;
  #endif
  if ( common.SqrtS < SumMasses ) return false;  // The model cannot work at low energy

  common.PZcms2 = ( sqr( common.S ) + sqr( common.M0projectile2 ) + sqr( common.M0target2 )
                    - 2.0 * ( common.S * ( common.M0projectile2 + common.M0target2 ) 
                              + common.M0projectile2 * common.M0target2 ) ) / 4.0 / common.S;
  #ifdef debugFTFexcitation
  G4cout << "PZcms2 after toBePutOnMassShell " << common.PZcms2 << G4endl;
  #endif
  if ( common.PZcms2 < 0.0 ) return false;  // It can be in an interaction with off-shell nuclear nucleon

  common.PZcms = std::sqrt( common.PZcms2 );
  if ( toBePutOnMassShell ) {
    if ( common.Pprojectile.z() > 0.0 ) {
      common.Pprojectile.setPz(  common.PZcms );
      common.Ptarget.setPz(     -common.PZcms );
    } else {
      common.Pprojectile.setPz( -common.PZcms );
      common.Ptarget.setPz(      common.PZcms );
    }
    common.Pprojectile.setE( std::sqrt( common.M0projectile2 
                                        + common.Pprojectile.x() * common.Pprojectile.x() 
                                        + common.Pprojectile.y() * common.Pprojectile.y() 
                                        + common.PZcms2 ) );
    common.Ptarget.setE( std::sqrt( common.M0target2 
                                    + common.Ptarget.x() * common.Ptarget.x()
                                    + common.Ptarget.y() * common.Ptarget.y() 
                                    + common.PZcms2 ) );
  }
  #ifdef debugFTFexcitation
  G4cout << "Start --------------------" << G4endl << "Proj M0 Mdif Mndif " << common.M0projectile
         << " " << common.ProjectileDiffStateMinMass << "  " << common.ProjectileNonDiffStateMinMass
         << G4endl
         << "Targ M0 Mdif Mndif " << common.M0target << " " << common.TargetDiffStateMinMass 
         << " " << common.TargetNonDiffStateMinMass << G4endl << "SqrtS " << common.SqrtS << G4endl
         << "Proj CMS " << common.Pprojectile << G4endl << "Targ CMS " << common.Ptarget << G4endl;
  #endif

  // Check for possible quark exchange
  ProjectileRapidity = common.Pprojectile.rapidity();
  TargetRapidity = common.Ptarget.rapidity();
  G4double QeNoExc = theParameters->GetProcProb( 0, ProjectileRapidity - TargetRapidity );
  G4double QeExc   = theParameters->GetProcProb( 1, ProjectileRapidity - TargetRapidity ) *
                     theParameters->GetProcProb( 4, ProjectileRapidity - TargetRapidity );
  common.ProbProjectileDiffraction = 
    theParameters->GetProcProb( 2, ProjectileRapidity - TargetRapidity );
  common.ProbTargetDiffraction = 
    theParameters->GetProcProb( 3, ProjectileRapidity - TargetRapidity );
  common.ProbOfDiffraction = common.ProbProjectileDiffraction + common.ProbTargetDiffraction;
  #ifdef debugFTFexcitation
  G4cout << "Proc Probs " << QeNoExc << " " << QeExc << " " 
         << common.ProbProjectileDiffraction << " " << common.ProbTargetDiffraction << G4endl 
         << "ProjectileRapidity " << ProjectileRapidity << G4endl;
  #endif

  if ( QeNoExc + QeExc + common.ProbProjectileDiffraction + common.ProbTargetDiffraction > 1.0 ) {
    QeNoExc = 1.0 - QeExc - common.ProbProjectileDiffraction - common.ProbTargetDiffraction;
  }
  if ( QeExc + QeNoExc != 0.0 ) {
    common.ProbExc = QeExc / ( QeExc + QeNoExc );
  }
  if ( 1.0 - QeExc - QeNoExc > 0.0 ) { 
    common.ProbProjectileDiffraction /= ( 1.0 - QeExc - QeNoExc );
    common.ProbTargetDiffraction     /= ( 1.0 - QeExc - QeNoExc );
  }
  #ifdef debugFTFexcitation
  G4cout << "Proc Probs " << QeNoExc << " " << QeExc << " " 
         << common.ProbProjectileDiffraction << " " << common.ProbTargetDiffraction << G4endl 
         << "ProjectileRapidity " << ProjectileRapidity << G4endl;
  #endif

  // Try out charge exchange
  G4int returnCode = 1;
  if ( G4UniformRand() < QeExc + QeNoExc ) {
    returnCode = ExciteParticipants_doChargeExchange( projectile, target, theParameters, 
                                                      theElastic, common );
  }

  G4bool returnResult = false;
  if ( returnCode == 0 ) {
    returnResult = true;  // Successfully ended: no need of extra work
  } else if ( returnCode == 1 ) {

    common.ProbOfDiffraction = common.ProbProjectileDiffraction + common.ProbTargetDiffraction;
    #ifdef debugFTFexcitation
    G4cout << "Excitation --------------------" << G4endl
           << "Proj M0 MdMin MndMin " << common.M0projectile << " " 
           << common.ProjectileDiffStateMinMass << "  " << common.ProjectileNonDiffStateMinMass
           << G4endl
           << "Targ M0 MdMin MndMin " << common.M0target << " " << common.TargetDiffStateMinMass
           << " " << common.TargetNonDiffStateMinMass << G4endl << "SqrtS " << common.SqrtS 
           << G4endl
           << "Prob: ProjDiff TargDiff + Sum " << common.ProbProjectileDiffraction << " " 
           << common.ProbTargetDiffraction << " " << common.ProbOfDiffraction << G4endl;
    #endif
    if ( common.ProbOfDiffraction != 0.0 ) {
      common.ProbProjectileDiffraction /= common.ProbOfDiffraction;
    } else {
      common.ProbProjectileDiffraction = 0.0;
    }
    #ifdef debugFTFexcitation
    G4cout << "Prob: ProjDiff TargDiff + Sum " << common.ProbProjectileDiffraction << " " 
           << common.ProbTargetDiffraction << " " << common.ProbOfDiffraction << G4endl;
    #endif
    common.ProjectileDiffStateMinMass2    = sqr( common.ProjectileDiffStateMinMass );
    common.ProjectileNonDiffStateMinMass2 = sqr( common.ProjectileNonDiffStateMinMass );
    common.TargetDiffStateMinMass2        = sqr( common.TargetDiffStateMinMass );
    common.TargetNonDiffStateMinMass2     = sqr( common.TargetNonDiffStateMinMass );
    // Choose between diffraction and non-diffraction process
    if ( G4UniformRand() < common.ProbOfDiffraction ) {
      
      returnResult = ExciteParticipants_doDiffraction( projectile, target, theParameters, common );
      
    } else {
      
      returnResult = ExciteParticipants_doNonDiffraction( projectile, target, theParameters, common );
      
    }
    if ( returnResult ) {
      common.Pprojectile += common.Qmomentum;
      common.Ptarget     -= common.Qmomentum;
      // Transform back and update SplitableHadron Participant.
      common.Pprojectile.transform( common.toLab );
      common.Ptarget.transform( common.toLab );
      #ifdef debugFTFexcitation
      G4cout << "Mproj " << common.Pprojectile.mag() << G4endl << "Mtarg " << common.Ptarget.mag()
             << G4endl;
      #endif
      projectile->Set4Momentum( common.Pprojectile );
      target->Set4Momentum( common.Ptarget );
      projectile->IncrementCollisionCount( 1 );
      target->IncrementCollisionCount( 1 );
    }
  }

  return returnResult;
}

//-----------------------------------------------------------------------------

G4int G4DiffractiveExcitation::
ExciteParticipants_doChargeExchange( G4VSplitableHadron*    projectile, 
                                     G4VSplitableHadron*    target,
                                     G4FTFParameters*       theParameters,
                                     G4ElasticHNScattering* theElastic,
                                     G4DiffractiveExcitation::CommonVariables& common ) const {
  // First of the three utility methods used only by ExciteParticipants: 
  // it does the sampling for the charge-exchange case.
  // This method returns an integer code - instead of a boolean, with the following meaning:
  //   "0" : successfully ended and nothing else needs to be done;
  //   "1" : successfully completed, but the work needs to be continued;
  //  "99" : unsuccessfully ended, nothing else can be done.
  G4int returnCode = 99;

  G4double DeltaProbAtQuarkExchange = theParameters->GetDeltaProbAtQuarkExchange();
  G4ParticleDefinition* TestParticle = 0;
  G4double MtestPr = 0.0, MtestTr = 0.0;

  #ifdef debugFTFexcitation
  G4cout << "Q exchange --------------------------" << G4endl;
  #endif

  // The target hadron is always a nucleon (i.e. either a proton or a neutron,
  // never an antinucleon), therefore only a quark (not an antiquark) can be
  // exchanged between the projectile hadron and the target hadron (otherwise
  // we could get a quark-quark-antiquark system which cannot be a bound state).
  // This implies that any projectile meson or anti-meson - given that it has
  // a constituent quark in all cases - can have charge exchange with a target
  // hadron. Instead, any projectile anti-baryon can never have charge exchange
  // with a target hadron (because it has only constituent anti-quarks);
  // projectile baryons, instead can have charge exchange with a target hadron.
  
  G4int NewProjCode = 0, NewTargCode = 0, ProjQ1 = 0, ProjQ2 = 0, ProjQ3  = 0;
  //  Projectile unpacking
  if ( common.absProjectilePDGcode < 1000 ) {  // projectile is meson 
    UnpackMeson(  common.ProjectilePDGcode, ProjQ1, ProjQ2 );  
  } else {                                     // projectile is baryon
    UnpackBaryon( common.ProjectilePDGcode, ProjQ1, ProjQ2, ProjQ3 );
  }
  //  Target unpacking
  G4int TargQ1 = 0, TargQ2 = 0, TargQ3 = 0;
  UnpackBaryon( common.TargetPDGcode, TargQ1, TargQ2, TargQ3 ); 
  #ifdef debugFTFexcitation
  G4cout << "Proj Quarks " << ProjQ1 << " " << ProjQ2 << " " << ProjQ3 << G4endl
         << "Targ Quarks " << TargQ1 << " " << TargQ2 << " " << TargQ3 << G4endl;
  #endif
  
  // Sampling of exchanged quarks
  G4int ProjExchangeQ = 0, TargExchangeQ = 0;
  if ( common.absProjectilePDGcode < 1000 ) {  // projectile is meson 

    G4bool isProjQ1Quark = false;
    ProjExchangeQ = ProjQ2;
    if ( ProjQ1 > 0 ) {  // ProjQ1 is a quark
      isProjQ1Quark = true;
      ProjExchangeQ = ProjQ1;
    }
    // Exchange of non-identical quarks is allowed
    G4int NpossibleStates = 0;
    if ( ProjExchangeQ != TargQ1 ) NpossibleStates++;
    if ( ProjExchangeQ != TargQ2 ) NpossibleStates++;
    if ( ProjExchangeQ != TargQ3 ) NpossibleStates++;  
    G4int Nsampled = (G4int)G4RandFlat::shootInt( G4long( NpossibleStates ) )+1;
    NpossibleStates = 0;
    if ( ProjExchangeQ != TargQ1 ) {
      if ( ++NpossibleStates == Nsampled ) {
        TargExchangeQ = TargQ1; TargQ1 = ProjExchangeQ; 
        isProjQ1Quark ? ProjQ1 = TargExchangeQ : ProjQ2 = TargExchangeQ;
      }
    }
    if ( ProjExchangeQ != TargQ2 ) {
      if ( ++NpossibleStates == Nsampled ) {
        TargExchangeQ = TargQ2; TargQ2 = ProjExchangeQ;
        isProjQ1Quark ? ProjQ1 = TargExchangeQ : ProjQ2 = TargExchangeQ;
      }
    }
    if ( ProjExchangeQ != TargQ3 ) {
      if ( ++NpossibleStates == Nsampled ) {
        TargExchangeQ = TargQ3; TargQ3 = ProjExchangeQ;
        isProjQ1Quark ? ProjQ1 = TargExchangeQ : ProjQ2 = TargExchangeQ;
      }
    }
    #ifdef debugFTFexcitation
    G4cout << "Exchanged Qs in Pr Tr " << ProjExchangeQ << " " << TargExchangeQ << G4endl;
    #endif

    G4int aProjQ1 = std::abs( ProjQ1 ), aProjQ2 = std::abs( ProjQ2 );
    G4bool ProjExcited = false;
    const G4int maxNumberOfAttempts = 50;
    G4int attempts = 0;
    while ( attempts++ < maxNumberOfAttempts ) {  /* Loop checking, 10.08.2015, A.Ribon */

      // Determination of a new projectile ID which satisfies energy-momentum conservation
      G4double ProbSpin0 = 0.5;
      G4double Ksi = G4UniformRand();
      if ( aProjQ1 == aProjQ2 ) {
        if ( G4UniformRand() < ProbSpin0 ) {  // Meson spin = 0 (pseudo-scalar)
          if ( aProjQ1 < 3 ) {
            NewProjCode = 111;                // pi0
            if ( Ksi < 0.5 ) {
              NewProjCode = 221;              // eta
              if ( Ksi < 0.25 ) {
                NewProjCode = 331;            // eta'
              }
            } 
          } else if ( aProjQ1 == 3 ) {
              NewProjCode = 221;              // eta
              if ( Ksi < 0.5 ) {
                NewProjCode = 331;            // eta'
              }
          } else if ( aProjQ1 == 4 ) {
	    NewProjCode = 441;                // eta_c(1S)
	  } else if ( aProjQ1 == 5 ) {
	    NewProjCode = 551;                // eta_b(1S)
	  }
        } else {                              // Meson spin = 1 (vector meson)
          if ( aProjQ1 < 3 ) {
            NewProjCode = 113;                // rho0
            if ( Ksi < 0.5 ) {
              NewProjCode = 223;              // omega
            } 
          } else if ( aProjQ1 == 3 ) {
            NewProjCode = 333;                // phi
          } else if ( aProjQ1 == 4 ) {
	    NewProjCode = 443;                // J/psi(1S)
	  } else if ( aProjQ1 == 5 ) {
	    NewProjCode = 553;                // Upsilon(1S)
	  }
        }
      } else {
        if ( aProjQ1 > aProjQ2 ) {
          NewProjCode = aProjQ1*100 + aProjQ2*10 + 1;
        } else {
          NewProjCode = aProjQ2*100 + aProjQ1*10 + 1;
        }
      }
      #ifdef debugFTFexcitation
      G4cout << "NewProjCode " << NewProjCode << G4endl;
      #endif
      // Decide (with 50% probability) whether the projectile hadrons is excited,
      // but not in the case of charmed and bottom hadrons (because in Geant4
      // there are no excited charmed and bottom states).
      ProjExcited = false;
      if ( aProjQ1 <= 3  &&  aProjQ2 <= 3  &&  G4UniformRand() < 0.5 ) {
        NewProjCode += 2;  // Excited meson (J=1 -> 2*J+1=2+1=3, last digit of the PDG code) 
        ProjExcited = true;
      }

      // So far we have used the absolute values of the PDG codes of the two constituent quarks:
      // now look at their signed values to set properly the signed of the meson's PDG code.
      G4int value = ProjQ1, absValue = aProjQ1, Qquarks = 0;
      for ( G4int iQ = 0; iQ < 2; ++iQ ) {  // 0 : first quark; 1 : second quark
        if ( iQ == 1 ) {
          value = ProjQ2; absValue = aProjQ2;
        }
        if ( absValue == 2  ||  absValue == 4 ) {  // u or c
          Qquarks += 2*value/absValue;  // u, c : positively charged +2 (e/3 unit)
        } else {
          Qquarks -= value/absValue;    // d, s, b : negatively charged -1 (e/3 unit)
        }
      }
      // If Qquarks is positive, the sign of NewProjCode is fine.
      // If Qquarks is negative, then the sign of NewProjCode needs to be reversed.
      // If Qquarks is zero, then we need to distinguish between 3 cases:
      // 1. If aProjQ1 and aProjQ2 are the same, then the sign of NewProjCode is fine
      //    (because the antiparticle is the same as the particle, e.g. eta, eta').
      // 2. If aProjQ1 and aProjQ2 are not the same, given that Qquarks is zero,
      //    we have only two possibilities:
      //    a. aProjQ1 and aProjQ2 are two different down-type quarks, i.e.
      //       (s,d) or (b,d), or (b,s). In this case, the sign of NewProjCode
      //       is fine (because the heaviest of the two down-type quarks has
      //       to be anti-quark belonging to the initial projectile, which
      //       implies a meson with positive PDG code, e.g. B0 (bbar,d), Bs (bbar,s).
      //    b. aProjQ1 and aProjQ2 are two different up-type quarks, i.e. (u,c).
      //       The heaviest of the two (c) has to be an anti-quark (cbar) left
      //       in the projectile, therefore the sign of NewProjCode needs to be
      //       reverse: 421 -> -421 anti_D0 (cbar,u)
      if ( Qquarks < 0  || ( Qquarks == 0  &&  aProjQ1 != aProjQ2  &&  aProjQ1%2 == 0 ) ) {
	NewProjCode *= -1;
      }
      #ifdef debugFTFexcitation
      G4cout << "NewProjCode +2 or 0 " << NewProjCode << G4endl;
      G4cout<<"+++++++++++++++++++++++++++++++++++++++"<<G4endl;
      G4cout<<ProjQ1<<" "<<ProjQ2<<" "<<Qquarks<<G4endl;
      G4cout<<"+++++++++++++++++++++++++++++++++++++++"<<G4endl;
      #endif

      // Proj 
      TestParticle = G4ParticleTable::GetParticleTable()->FindParticle( NewProjCode );
      if ( ! TestParticle ) continue;
      common.MminProjectile = common.BrW.GetMinimumMass( TestParticle );
      if ( common.SqrtS - common.M0target < common.MminProjectile ) continue;
      MtestPr = common.BrW.SampleMass( TestParticle, TestParticle->GetPDGMass() 
                                                     + 5.0*TestParticle->GetPDGWidth() );
      #ifdef debugFTFexcitation
      G4cout << "TestParticle Name " << NewProjCode << " " << TestParticle->GetParticleName()
             << G4endl
             << "MtestPart MtestPart0 "<<MtestPr<<" "<<TestParticle->GetPDGMass()<<G4endl
             << "M0projectile projectile PDGMass " << common.M0projectile << " " 
             << projectile->GetDefinition()->GetPDGMass() << G4endl;
      #endif

      // Targ
      NewTargCode = NewNucleonId( TargQ1, TargQ2, TargQ3 );
      #ifdef debugFTFexcitation
      G4cout << "New TrQ " << TargQ1 << " " << TargQ2 << " " << TargQ3 << G4endl
             << "NewTargCode " << NewTargCode << G4endl;
      #endif
      
      // If the target has not heavy (charm or bottom) constituent quark,
      // see whether a Delta isobar can be created.
      if ( TargQ1 <= 3  &&  TargQ2 <= 3  &&  TargQ3 <= 3 ) {
        if ( TargQ1 != TargQ2  &&  TargQ1 != TargQ3  &&  TargQ2 != TargQ3 ) {  // Lambda or Sigma0 ?
          if ( G4UniformRand() < 0.5 ) {
            NewTargCode += 2;
          } else if ( G4UniformRand() < 0.75 ) {
            NewTargCode = 3122;  // Lambda
          }
        } else if ( TargQ1 == TargQ2  &&  TargQ1 == TargQ3 ) {
          NewTargCode += 2; ProjExcited = true;                         // Create Delta isobar
        } else if ( target->GetDefinition()->GetPDGiIsospin() == 3 ) {  // Delta was the target
          if ( G4UniformRand() > DeltaProbAtQuarkExchange ) {
            NewTargCode += 2; ProjExcited = true;
          }
        } else if ( ! ProjExcited  &&
                    G4UniformRand() < DeltaProbAtQuarkExchange  &&      // Nucleon was the target
                    common.SqrtS > common.M0projectile +                // Delta mass
	  	                   G4ParticleTable::GetParticleTable()->FindParticle( 2224 )->GetPDGMass() ) {
          NewTargCode += 2;  // Create Delta isobar
        }
      }

      // Protection against:
      // - Lambda* (i.e. excited Lambda state)         NOT existing in PDG ,   ->  Lambda
      // - Sigma*  (i.e. excited Sigma hyperon states) NOT existing in Geant4  ->  Sigma
      // - Xi*     (i.e. excited Xi hyperon states)    NOT existing in Geant4  ->  Xi
      if ( NewTargCode == 3124  ||   // Lambda* NOT existing in PDG !
	   NewTargCode == 3224  ||   // Sigma*+ NOT existing in Geant4
	   NewTargCode == 3214  ||   // Sigma*0 NOT existing in Geant4
	   NewTargCode == 3114  ||   // Sigma*- NOT existing in Geant4
	   NewTargCode == 3324  ||   // Xi*0    NOT existing in Geant4
	   NewTargCode == 3314  ) {  // Xi*-    NOT existing in Geant4
	//G4cout << "G4DiffractiveExcitation::ExciteParticipants_doChargeExchange INEXISTING TARGET PDGcode="
	//       << NewTargCode << G4endl;
	NewTargCode -= 2;  // Corresponding ground-state hyperon
      }
      
      // Special treatment for charmed and bottom baryons : in Geant4 there are no Xi_c' and Xi_b'
      // so we need to transform them by hand to Xi_c and Xi_b, respectively.
      #ifdef debug_heavyHadrons
      G4int initialNewTargCode = NewTargCode;
      #endif
      if      ( NewTargCode == 4322 ) NewTargCode = 4232;  // Xi_c'+  ->  Xi_c+
      else if ( NewTargCode == 4312 ) NewTargCode = 4132;  // Xi_c'0  ->  Xi_c0
      else if ( NewTargCode == 5312 ) NewTargCode = 5132;  // Xi_b'-  ->  Xi_b-
      else if ( NewTargCode == 5322 ) NewTargCode = 5232;  // Xi_b'0  ->  Xi_b0
      #ifdef debug_heavyHadrons
      if ( NewTargCode != initialNewTargCode ) {
	G4cout << "G4DiffractiveExcitation::ExciteParticipants_doChargeExchange : forcing (inexisting in G4)" << G4endl
	       << "\t target heavy baryon with pdgCode=" << initialNewTargCode
	       << " into pdgCode=" << NewTargCode << G4endl;
      }
      #endif
      
      TestParticle = G4ParticleTable::GetParticleTable()->FindParticle( NewTargCode );
      if ( ! TestParticle ) continue;
      #ifdef debugFTFexcitation
      G4cout << "New targ " << NewTargCode << " " << TestParticle->GetParticleName() << G4endl;
      #endif
      common.MminTarget = common.BrW.GetMinimumMass( TestParticle );
      if ( common.SqrtS - MtestPr < common.MminTarget ) continue;
      MtestTr = common.BrW.SampleMass( TestParticle, TestParticle->GetPDGMass() 
                                                     + 5.0*TestParticle->GetPDGWidth() );
      if ( common.SqrtS > MtestPr + MtestTr ) break;

    }  // End of while loop
    if ( attempts >= maxNumberOfAttempts ) return returnCode;  // unsuccessfully ended, nothing else can be done

    if ( MtestPr >= common.Pprojectile.mag()  ||  projectile->GetStatus() != 0 ) { 
      common.M0projectile = MtestPr;
    }
    #ifdef debugFTFexcitation
    G4cout << "M0projectile After check " << common.M0projectile << G4endl;
    #endif
    common.M0projectile2 = common.M0projectile * common.M0projectile;
    common.ProjectileDiffStateMinMass    = common.M0projectile + 220.0*MeV;  // 220 MeV=m_pi+80 MeV 
    common.ProjectileNonDiffStateMinMass = common.M0projectile + 220.0*MeV;  // 220 MeV=m_pi+80 MeV
    if ( MtestTr >= common.Ptarget.mag()  ||  target->GetStatus() != 0 ) {
      common.M0target = MtestTr;
    }
    common.M0target2 = common.M0target * common.M0target;
    #ifdef debugFTFexcitation
    G4cout << "New targ M0 M0^2 " << common.M0target << " " << common.M0target2 << G4endl;
    #endif
    common.TargetDiffStateMinMass    = common.M0target + 220.0*MeV;          // 220 MeV=m_pi+80 MeV;    
    common.TargetNonDiffStateMinMass = common.M0target + 220.0*MeV;          // 220 MeV=m_pi+80 MeV; 

  } else {  // of the if ( common.absProjectilePDGcode < 1000 ) ; the projectile is baryon now

    // If it is a projectile anti-baryon, no quark exchange is possible with a target hadron,
    // therefore returns immediately 1 (which means "successfully completed, but the work
    // needs to be continued").
    if ( common.ProjectilePDGcode < 0 ) return 1;

    // Choose randomly, with equal probability, whether to consider the quarks of the 
    // projectile or target hadron for selecting the flavour of the exchanged quark.
    G4bool isProjectileExchangedQ = false;
    G4int firstQ      = TargQ1, secondQ      = TargQ2, thirdQ      = TargQ3;
    G4int otherFirstQ = ProjQ1, otherSecondQ = ProjQ2, otherThirdQ = ProjQ3;
    if ( G4UniformRand() < 0.5 ) {
      isProjectileExchangedQ = true;
      firstQ      = ProjQ1; secondQ      = ProjQ2; thirdQ      = ProjQ3;
      otherFirstQ = TargQ1; otherSecondQ = TargQ2; otherThirdQ = TargQ3;
    }
    // Choose randomly, with equal probability, which of the three quarks of the
    // selected (projectile or target) hadron is the exhanged quark.
    G4int exchangedQ = 0;
    G4double Ksi = G4UniformRand();
    if ( Ksi < 0.333333 ) {
      exchangedQ = firstQ;
    } else if ( 0.333333 <= Ksi  &&  Ksi < 0.666667 ) {
      exchangedQ = secondQ;
    } else {
      exchangedQ = thirdQ;
    }
    #ifdef debugFTFexcitation
    G4cout << "Exchange Qs isProjectile Q " << isProjectileExchangedQ << " " << exchangedQ << " ";
    #endif
    // The exchanged quarks (one of the projectile hadron and one of the target hadron)
    // are always accepted if they have different flavour, else (i.e. when they have the
    // same flavour) they are accepted only with a specified probability.
    G4double probSame = theParameters->GetProbOfSameQuarkExchange();
    const G4int MaxCount = 100;
    G4int count = 0, otherExchangedQ = 0; 
    do {
      if ( exchangedQ != otherFirstQ  ||  G4UniformRand() < probSame ) {
        otherExchangedQ = otherFirstQ; otherFirstQ = exchangedQ; exchangedQ = otherExchangedQ;
      } else {
        if ( exchangedQ != otherSecondQ  ||  G4UniformRand() < probSame ) {
          otherExchangedQ = otherSecondQ; otherSecondQ = exchangedQ; exchangedQ = otherExchangedQ;
        } else {
          if ( exchangedQ != otherThirdQ  ||  G4UniformRand() < probSame ) {
            otherExchangedQ = otherThirdQ; otherThirdQ = exchangedQ; exchangedQ = otherExchangedQ;
          }
        }
      }
    } while ( otherExchangedQ == 0  &&  ++count < MaxCount );
    if ( count >= MaxCount ) return returnCode;  // All attempts failed: unsuccessfully ended, nothing else can be done 
    // Swap (between projectile and target hadron) the two quarks that have been sampled
    // as "exchanged" quarks.
    if ( Ksi < 0.333333 ) {
      firstQ = exchangedQ;
    } else if ( 0.333333 <= Ksi  &&  Ksi < 0.666667 ) {
      secondQ = exchangedQ;
    } else {
      thirdQ = exchangedQ;
    }
    if ( isProjectileExchangedQ ) {
      ProjQ1 = firstQ;      ProjQ2 = secondQ;      ProjQ3 = thirdQ;
      TargQ1 = otherFirstQ; TargQ2 = otherSecondQ; TargQ3 = otherThirdQ;
    } else {
      TargQ1 = firstQ;      TargQ2 = secondQ;      TargQ3 = thirdQ;
      ProjQ1 = otherFirstQ; ProjQ2 = otherSecondQ; ProjQ3 = otherThirdQ;
    }
    #ifdef debugFTFexcitation
    G4cout << "Exchange Qs Pr  Tr " << ( isProjectileExchangedQ ? exchangedQ : otherExchangedQ )
           << " " << ( isProjectileExchangedQ ? otherExchangedQ : exchangedQ ) << G4endl;
    #endif

    NewProjCode = NewNucleonId( ProjQ1, ProjQ2, ProjQ3 );
    NewTargCode = NewNucleonId( TargQ1, TargQ2, TargQ3 );
    // Decide whether the new projectile hadron is a Delta particle; 
    // then decide whether the new target hadron is a Delta particle.
    // Notice that a Delta particle has the last PDG digit "4" (because its spin is 3/2),
    // whereas a nucleon has "2" (because its spin is 1/2).
    // If the new projectile hadron or the new target hadron has a heavy (c or b)
    // constituent quark, then skip this part (because Geant4 does not have
    // excited charm and bottom hadrons).
    for ( G4int iHadron = 0; iHadron < 2; iHadron++ ) {
      // First projectile hadron, then target hadron
      G4int codeQ1 = ProjQ1, codeQ2 = ProjQ2, codeQ3 = ProjQ3, newHadCode = NewProjCode;
      G4double massConstraint = common.M0target;
      G4bool isHadronADelta = ( projectile->GetDefinition()->GetPDGiIsospin() == 3 );
      if ( iHadron == 1 ) {  // Target hadron
        codeQ1 = TargQ1, codeQ2 = TargQ2, codeQ3 = TargQ3, newHadCode = NewTargCode;
        massConstraint = common.M0projectile;
        isHadronADelta = ( target->GetDefinition()->GetPDGiIsospin() == 3 );
      }
      if ( codeQ1 > 3  ||  codeQ2 > 3  ||  codeQ3 > 3 ) continue;  // No excited charm or bottom states in Geant4
      if ( codeQ1 == codeQ2  &&  codeQ1 == codeQ3 ) {  // The three quarks are the same
        newHadCode += 2;  // Delta++ (uuu) or Delta- (ddd) : spin 3/2, last PDG digit = 4 (so +2 wrt p or n)
      } else if ( isHadronADelta ) {  // Hadron (projectile or target) was Delta
        if ( G4UniformRand() > DeltaProbAtQuarkExchange ) { 
          newHadCode += 2;  // Delta+ (uud) or Delta0 (udd) : spin 3/2, last PDG digit = 4 (so +2 wrt p or n)
        } else {
          newHadCode += 0;  // No delta (so the last PDG digit remains 2)
        }
      } else {  // Hadron (projectile or target) was Nucleon
        if ( G4UniformRand() < DeltaProbAtQuarkExchange  &&  
             common.SqrtS > G4ParticleTable::GetParticleTable()->FindParticle( 2224 )->GetPDGMass() 
                            + massConstraint ) {
          newHadCode += 2;  // Delta+ (uud) or Delta0 (udd) : spin 3/2, last PDG digit = 4 (so +2 wrt p or n)
        } else { 
          newHadCode += 0;  // No delta (so the last PDG digit remains 2)
        }
      }
      if ( iHadron == 0 ) {  // Projectile hadron
        NewProjCode = newHadCode;
      } else {               // Target hadron
        NewTargCode = newHadCode;
      }
    }
    #ifdef debugFTFexcitation
    G4cout << "NewProjCode NewTargCode " << NewProjCode << " " << NewTargCode << G4endl;
    #endif

    // Protection against:
    // - Lambda* (i.e. excited Lambda state)         NOT existing in PDG ,   ->  Lambda
    // - Sigma*  (i.e. excited Sigma hyperon states) NOT existing in Geant4  ->  Sigma
    // - Xi*     (i.e. excited Xi hyperon states)    NOT existing in Geant4  ->  Xi
    if ( NewProjCode == 3124  ||   // Lambda* NOT existing in PDG !
	 NewProjCode == 3224  ||   // Sigma*+ NOT existing in Geant4
	 NewProjCode == 3214  ||   // Sigma*0 NOT existing in Geant4
	 NewProjCode == 3114  ||   // Sigma*- NOT existing in Geant4
	 NewProjCode == 3324  ||   // Xi*0    NOT existing in Geant4
	 NewProjCode == 3314  ) {  // Xi*-    NOT existing in Geant4
      //G4cout << "G4DiffractiveExcitation::ExciteParticipants_doChargeExchange INEXISTING PROJECTILE PDGcode="
      //       << NewProjCode << G4endl;
      NewProjCode -= 2;  // Corresponding ground-state hyperon
    }
    if ( NewTargCode == 3124  ||   // Lambda* NOT existing in PDG !
	 NewTargCode == 3224  ||   // Sigma*+ NOT existing in Geant4
	 NewTargCode == 3214  ||   // Sigma*0 NOT existing in Geant4
	 NewTargCode == 3114  ||   // Sigma*- NOT existing in Geant4
	 NewTargCode == 3324  ||   // Xi*0    NOT existing in Geant4
	 NewTargCode == 3314  ) {  // Xi*-    NOT existing in Geant4
      //G4cout << "G4DiffractiveExcitation::ExciteParticipants_doChargeExchange INEXISTING TARGET PDGcode="
      //       << NewTargCode << G4endl;
      NewTargCode -= 2;  // Corresponding ground-state hyperon
    }

    // Special treatment for charmed and bottom baryons : in Geant4 there are no Xi_c' and Xi_b'
    // so we need to transform them by hand to the, respectively, Xi_c and Xi_b.
    #ifdef debug_heavyHadrons
    G4int initialNewProjCode = NewProjCode, initialNewTargCode = NewTargCode;
    #endif
    if      ( NewProjCode == 4322 ) NewProjCode = 4232;  // Xi_c'+  ->  Xi_c+
    else if ( NewProjCode == 4312 ) NewProjCode = 4132;  // Xi_c'0  ->  Xi_c0
    else if ( NewProjCode == 5312 ) NewProjCode = 5132;  // Xi_b'-  ->  Xi_b-
    else if ( NewProjCode == 5322 ) NewProjCode = 5232;  // Xi_b'0  ->  Xi_b0
    if      ( NewTargCode == 4322 ) NewTargCode = 4232;  // Xi_c'+  ->  Xi_c+
    else if ( NewTargCode == 4312 ) NewTargCode = 4132;  // Xi_c'0  ->  Xi_c0
    else if ( NewTargCode == 5312 ) NewTargCode = 5132;  // Xi_b'-  ->  Xi_b-
    else if ( NewTargCode == 5322 ) NewTargCode = 5232;  // Xi_b'0  ->  Xi_b0
    #ifdef debug_heavyHadrons
    if ( NewProjCode != initialNewProjCode  ||  NewTargCode != initialNewTargCode ) {
      G4cout << "G4DiffractiveExcitation::ExciteParticipants_doChargeExchange : forcing (inexisting in G4)" << G4endl
             << "\t heavy baryon into an existing one:" << G4endl;
      if ( NewProjCode != initialNewProjCode ) {
	G4cout << "\t \t projectile baryon with pdgCode=" << initialNewProjCode
	       << " into pdgCode=" << NewProjCode << G4endl;
      }
      if ( NewTargCode != initialNewTargCode ) {
        G4cout << "\t \t target baryon with pdgCode=" << initialNewTargCode
	       << " into pdgCode=" << NewTargCode << G4endl;
      }
    }
    #endif

    // Sampling of the masses of the projectile and target nucleons.
    // Because of energy conservation, the ordering of the sampling matters:
    // randomly, half of the time we sample first the target nucleon mass and
    // then the projectile nucleon mass, and the other half of the time we
    // sample first the projectile nucleon mass and then the target nucleon mass.
    G4VSplitableHadron* firstHadron  = target;
    G4VSplitableHadron* secondHadron = projectile;
    G4int firstHadronCode = NewTargCode, secondHadronCode = NewProjCode;
    G4double massConstraint = common.M0projectile;
    G4bool isFirstTarget = true;
    if ( G4UniformRand() < 0.5 ) {
      // Sample first the projectile nucleon mass, then the target nucleon mass.
      firstHadron  = projectile;      secondHadron = target;
      firstHadronCode = NewProjCode;  secondHadronCode = NewTargCode;
      massConstraint = common.M0target;
      isFirstTarget = false;
    }
    G4double Mtest1st = 0.0, Mtest2nd = 0.0, Mmin1st = 0.0, Mmin2nd = 0.0;
    for ( int iSamplingCase = 0; iSamplingCase < 2; iSamplingCase++ ) {
      G4VSplitableHadron* aHadron = firstHadron;
      G4int aHadronCode = firstHadronCode;
      if ( iSamplingCase == 1 ) {  // Second nucleon mass sampling
        aHadron = secondHadron; aHadronCode = secondHadronCode; massConstraint = Mtest1st;
      }
      G4double MtestHadron = 0.0, MminHadron = 0.0;
      if ( aHadron->GetStatus() == 1  ||  aHadron->GetStatus() == 2 ) { 
        TestParticle = G4ParticleTable::GetParticleTable()->FindParticle( aHadronCode );
        if ( ! TestParticle ) return returnCode;  // Not possible to find such a hadron: unsuccessfully ended, nothing else can be done
        MminHadron = common.BrW.GetMinimumMass( TestParticle );
        if ( common.SqrtS - massConstraint < MminHadron ) return returnCode;  // Kinematically impossible: unsuccessfully ended, nothing else can be done
        if ( TestParticle->GetPDGWidth() == 0.0 ) { 
          MtestHadron = common.BrW.SampleMass( TestParticle, TestParticle->GetPDGMass() );
        } else {
          const G4int maxNumberOfAttempts = 50;
          G4int attempts = 0;
          while ( attempts < maxNumberOfAttempts ) {
            attempts++;
            MtestHadron = common.BrW.SampleMass( TestParticle, TestParticle->GetPDGMass() 
                                                 + 5.0*TestParticle->GetPDGWidth() );
            if ( common.SqrtS < MtestHadron + massConstraint ) {  
              continue;  // Kinematically unacceptable: try again
            } else {  
              break;     // Kinematically acceptable: the mass sampling is successful 
            }
          }
          if ( attempts >= maxNumberOfAttempts ) return returnCode;  // All attempts failed: unsuccessfully ended, nothing else can be done
        }
      }
      if ( iSamplingCase == 0 ) {  
        Mtest1st = MtestHadron;  Mmin1st = MminHadron;
      } else {
        Mtest2nd = MtestHadron;  Mmin2nd = MminHadron;
      }
    }  // End for loop on the two sampling cases (1st and 2nd)
    if ( isFirstTarget ) {
      MtestTr = Mtest1st;    MtestPr = Mtest2nd;
      common.MminTarget = Mmin1st;  common.MminProjectile = Mmin2nd;
    } else {
      MtestTr = Mtest2nd;    MtestPr = Mtest1st;
      common.MminTarget = Mmin2nd;  common.MminProjectile = Mmin1st;
    }
 
    if ( MtestPr != 0.0 ) {
      common.M0projectile = MtestPr;
      common.M0projectile2 = common.M0projectile * common.M0projectile;
      common.ProjectileDiffStateMinMass =    common.M0projectile + 220.0*MeV;  // 220 MeV=m_pi+80 MeV
      common.ProjectileNonDiffStateMinMass = common.M0projectile + 220.0*MeV;  // 220 MeV=m_pi+80 MeV
    }      
    if ( MtestTr != 0.0 ) {
      common.M0target = MtestTr;
      common.M0target2 = common.M0target * common.M0target;
      common.TargetDiffStateMinMass    = common.M0target + 220.0*MeV;          // 220 MeV=m_pi+80 MeV;    
      common.TargetNonDiffStateMinMass = common.M0target + 220.0*MeV;          // 220 MeV=m_pi+80 MeV; 
    }

  }  // End of if ( common.absProjectilePDGcode < 1000 )

  // If we assume that final state hadrons after the charge exchange will be
  // in the ground states, we have to put 
  if ( common.SqrtS < common.M0projectile + common.M0target ) return returnCode;  // unsuccessfully ended, nothing else can be done

  common.PZcms2 = ( sqr( common.S ) + sqr( common.M0projectile2 ) + sqr( common.M0target2 )
                    - 2.0 * ( common.S * ( common.M0projectile2 + common.M0target2 )
                              + common.M0projectile2 * common.M0target2 ) ) / 4.0 / common.S;
  #ifdef debugFTFexcitation
  G4cout << "At the end// NewProjCode " << NewProjCode << G4endl 
         << "At the end// NewTargCode " << NewTargCode << G4endl
         << "M0pr  M0tr  SqS " << common.M0projectile << " " << common.M0target << " " 
         << common.SqrtS << G4endl
         << "M0pr2 M0tr2 SqS " << common.M0projectile2 << " " << common.M0target2 << " "
         << common.SqrtS << G4endl
         << "PZcms2 after the change " << common.PZcms2 << G4endl << G4endl;
  #endif
  if ( common.PZcms2 < 0.0 ) return returnCode; // It can be if energy is not sufficient for Delta;
                                                // unsuccessfully ended, nothing else can be done
  projectile->SetDefinition( G4ParticleTable::GetParticleTable()->FindParticle( NewProjCode ) );
  target->SetDefinition( G4ParticleTable::GetParticleTable()->FindParticle( NewTargCode ) );
  common.PZcms = std::sqrt( common.PZcms2 );
  common.Pprojectile.setPz( common.PZcms );
  common.Pprojectile.setE( std::sqrt( common.M0projectile2 + common.PZcms2 ) );
  common.Ptarget.setPz(    -common.PZcms );
  common.Ptarget.setE( std::sqrt( common.M0target2 + common.PZcms2 ) );
  #ifdef debugFTFexcitation
  G4cout << "Proj Targ and Proj+Targ in CMS" << G4endl << common.Pprojectile << G4endl 
         << common.Ptarget << G4endl << common.Pprojectile + common.Ptarget << G4endl;
  #endif

  if ( projectile->GetStatus() != 0 ) projectile->SetStatus( 2 );
  if ( target->GetStatus()     != 0 ) target->SetStatus( 2 );

  // Check for possible excitation of the participants
  if ( common.SqrtS < common.M0projectile + common.TargetDiffStateMinMass  ||
       common.SqrtS < common.ProjectileDiffStateMinMass + common.M0target  ||
       common.ProbOfDiffraction == 0.0 ) common.ProbExc = 0.0;

  if ( G4UniformRand() > common.ProbExc ) {  // Make elastic scattering
    #ifdef debugFTFexcitation
    G4cout << "Make elastic scattering of new hadrons" << G4endl;
    #endif
    common.Pprojectile.transform( common.toLab );
    common.Ptarget.transform( common.toLab );
    projectile->Set4Momentum( common.Pprojectile );
    target->Set4Momentum( common.Ptarget );
    G4bool Result = theElastic->ElasticScattering( projectile, target, theParameters );
    #ifdef debugFTFexcitation
    G4cout << "Result of el. scatt " << Result << G4endl << "Proj Targ and Proj+Targ in Lab"
           << G4endl << projectile->Get4Momentum() << G4endl << target->Get4Momentum() << G4endl
           << projectile->Get4Momentum() + target->Get4Momentum() << " " 
           << (projectile->Get4Momentum() + target->Get4Momentum()).mag() << G4endl;
    #endif
    if ( Result ) returnCode = 0;  // successfully ended and nothing else needs to be done 
    return returnCode;
  }

  #ifdef debugFTFexcitation
  G4cout << "Make excitation of new hadrons" << G4endl;
  #endif

  // Redefinition of ProbOfDiffraction because the probabilities are changed due to quark exchange
  common.ProbOfDiffraction = common.ProbProjectileDiffraction + common.ProbTargetDiffraction;
  if ( common.ProbOfDiffraction != 0.0 ) {
    common.ProbProjectileDiffraction /= common.ProbOfDiffraction;
    common.ProbTargetDiffraction     /= common.ProbOfDiffraction;
  }

  return returnCode = 1;  // successfully completed, but the work needs to be continued
}

//-----------------------------------------------------------------------------

G4bool G4DiffractiveExcitation::
ExciteParticipants_doDiffraction( G4VSplitableHadron* projectile, G4VSplitableHadron* target,
                                  G4FTFParameters* theParameters, 
                                  G4DiffractiveExcitation::CommonVariables& common ) const {
  // Second of the three utility methods used only by ExciteParticipants: 
  // it does the sampling for the diffraction case, either projectile or target diffraction.

  G4bool isProjectileDiffraction = false;
  if ( G4UniformRand() < common.ProbProjectileDiffraction ) {  // projectile diffraction
    isProjectileDiffraction = true;
    #ifdef debugFTFexcitation
    G4cout << "projectile diffraction" << G4endl;
    #endif
    common.ProjMassT2 = common.ProjectileDiffStateMinMass2;
    common.ProjMassT  = common.ProjectileDiffStateMinMass;
    common.TargMassT2 = common.M0target2;
    common.TargMassT  = common.M0target;
  } else {                                                     // target diffraction
    #ifdef debugFTFexcitation
    G4cout << "Target diffraction" << G4endl;
    #endif
    common.ProjMassT2 = common.M0projectile2;
    common.ProjMassT  = common.M0projectile; 
    common.TargMassT2 = common.TargetDiffStateMinMass2;
    common.TargMassT  = common.TargetDiffStateMinMass;
  }

  // Check whether the interaction is possible
  if ( common.SqrtS < common.ProjMassT + common.TargMassT ) return false;
  
  common.PZcms2 = ( sqr( common.S ) + sqr( common.ProjMassT2 ) + sqr( common.TargMassT2 )
                    - 2.0 * ( common.S * ( common.ProjMassT2 + common.TargMassT2 )
                          + common.ProjMassT2 * common.TargMassT2 ) ) / 4.0 / common.S;
  if ( common.PZcms2 < 0.0 ) return false;
  common.maxPtSquare = common.PZcms2;

  G4double DiffrAveragePt2 = theParameters->GetAvaragePt2ofElasticScattering()*1.2;
  G4bool loopCondition = true;
  G4int whilecount = 0;
  do {  // Generate pt and mass of projectile

    whilecount++;
    if ( whilecount > 1000 ) {
      common.Qmomentum = G4LorentzVector( 0.0, 0.0, 0.0, 0.0 );
      return false;  //  Ignore this interaction
    };

    common.Qmomentum = G4LorentzVector( GaussianPt( DiffrAveragePt2, common.maxPtSquare ), 0 );
    common.Pt2 = G4ThreeVector( common.Qmomentum.vect() ).mag2();
    if ( isProjectileDiffraction ) {  // projectile diffraction
      common.ProjMassT2 = common.ProjectileDiffStateMinMass2 + common.Pt2;
      common.TargMassT2 = common.M0target2 + common.Pt2;
    } else {                          // target diffraction
      common.ProjMassT2 = common.M0projectile2 + common.Pt2;
      common.TargMassT2 = common.TargetDiffStateMinMass2 + common.Pt2;
    }
    common.ProjMassT = std::sqrt( common.ProjMassT2 );
    common.TargMassT = std::sqrt( common.TargMassT2 );
    if ( common.SqrtS < common.ProjMassT + common.TargMassT ) continue;

    common.PZcms2 = ( sqr( common.S ) + sqr( common.ProjMassT2 ) + sqr( common.TargMassT2 )
                      - 2.0 * ( common.S * ( common.ProjMassT2 + common.TargMassT2 )
                                + common.ProjMassT2 * common.TargMassT2 ) ) / 4.0 / common.S;
    if ( common.PZcms2 < 0.0 ) continue;

    common.PZcms = std::sqrt( common.PZcms2 );
    if ( isProjectileDiffraction ) {  // projectile diffraction
      common.PMinusMin = std::sqrt( common.ProjMassT2 + common.PZcms2 ) - common.PZcms; 
      common.PMinusMax = common.SqrtS - common.TargMassT;
      common.PMinusNew = ChooseP( common.PMinusMin, common.PMinusMax );
      common.TMinusNew = common.SqrtS - common.PMinusNew;
      common.Qminus = common.Ptarget.minus() - common.TMinusNew;
      common.TPlusNew = common.TargMassT2 / common.TMinusNew;
      common.Qplus = common.Ptarget.plus() - common.TPlusNew;
      common.Qmomentum.setPz( (common.Qplus - common.Qminus)/2.0 );
      common.Qmomentum.setE(  (common.Qplus + common.Qminus)/2.0 );
      loopCondition = ( ( common.Pprojectile + common.Qmomentum ).mag2() <  
                        common.ProjectileDiffStateMinMass2 ); 
    } else {                          // target diffraction
      common.TPlusMin = std::sqrt( common.TargMassT2 + common.PZcms2 ) - common.PZcms;
      common.TPlusMax = common.SqrtS - common.ProjMassT;
      common.TPlusNew = ChooseP( common.TPlusMin, common.TPlusMax );
      common.PPlusNew = common.SqrtS - common.TPlusNew;
      common.Qplus = common.PPlusNew - common.Pprojectile.plus();
      common.PMinusNew = common.ProjMassT2 / common.PPlusNew;
      common.Qminus = common.PMinusNew - common.Pprojectile.minus();
      common.Qmomentum.setPz( (common.Qplus - common.Qminus)/2.0 );
      common.Qmomentum.setE(  (common.Qplus + common.Qminus)/2.0 );
      loopCondition =  ( ( common.Ptarget - common.Qmomentum ).mag2() < 
                         common.TargetDiffStateMinMass2 );
    }

  } while ( loopCondition );  /* Loop checking, 10.08.2015, A.Ribon */
          // Repeat the sampling because there was not any excitation

  if ( isProjectileDiffraction ) {  // projectile diffraction
    projectile->SetStatus( 0 );
    if ( projectile->GetStatus() == 2 ) projectile->SetStatus( 1 );
    if ( target->GetStatus() == 1  &&  target->GetSoftCollisionCount() == 0 ) target->SetStatus( 2 );
  } else {                          // target diffraction
    target->SetStatus( 0 );
  }

  return true;
}
 
//-----------------------------------------------------------------------------

G4bool G4DiffractiveExcitation::
ExciteParticipants_doNonDiffraction( G4VSplitableHadron* projectile,
                                     G4VSplitableHadron* target,
                                     G4FTFParameters*    theParameters,
                                     G4DiffractiveExcitation::CommonVariables& common ) const {
  // Third of the three utility methods used only by ExciteParticipants: 
  // it does the sampling for the non-diffraction case.

  #ifdef debugFTFexcitation
  G4cout << "Non-diffraction process" << G4endl;
  #endif

  // Check whether the interaction is possible
  common.ProjMassT2 = common.ProjectileNonDiffStateMinMass2;
  common.ProjMassT  = common.ProjectileNonDiffStateMinMass;
  common.TargMassT2 = common.TargetNonDiffStateMinMass2;
  common.TargMassT  = common.TargetNonDiffStateMinMass;
  if ( common.SqrtS < common.ProjMassT + common.TargMassT ) return false;

  common.PZcms2 = ( sqr( common.S ) + sqr( common.ProjMassT2 ) + sqr( common.TargMassT2 )
                  - 2.0 * ( common.S * ( common.ProjMassT2 + common.TargMassT2 )
                        + common.ProjMassT2 * common.TargMassT2 ) ) / 4.0 / common.S;
  common.maxPtSquare = common.PZcms2;
  
  G4int whilecount = 0;
  do {  // Generate pt and masses

    whilecount++;
    if ( whilecount > 1000 ) {
      common.Qmomentum = G4LorentzVector( 0.0, 0.0, 0.0, 0.0 );
      return false;  // Ignore this interaction
    };

    common.Qmomentum = G4LorentzVector( GaussianPt( theParameters->GetAveragePt2(), 
                                                    common.maxPtSquare ), 0 );
    common.Pt2 = G4ThreeVector( common.Qmomentum.vect() ).mag2();
    common.ProjMassT2 = common.ProjectileNonDiffStateMinMass2 + common.Pt2;
    common.ProjMassT  = std::sqrt( common.ProjMassT2 );
    common.TargMassT2 = common.TargetNonDiffStateMinMass2 + common.Pt2;
    common.TargMassT  = std::sqrt( common.TargMassT2 );
    if ( common.SqrtS < common.ProjMassT + common.TargMassT ) continue;

    common.PZcms2 =( sqr( common.S ) + sqr( common.ProjMassT2 ) + sqr( common.TargMassT2 )
                     - 2.0 * ( common.S * ( common.ProjMassT2 + common.TargMassT2 )
                               + common.ProjMassT2 * common.TargMassT2 ) ) / 4.0 / common.S;
    if ( common.PZcms2 < 0.0 ) continue;

    common.PZcms = std::sqrt( common.PZcms2 );
    common.PMinusMin = std::sqrt( common.ProjMassT2 + common.PZcms2 ) - common.PZcms;
    common.PMinusMax = common.SqrtS - common.TargMassT;
    common.TPlusMin = std::sqrt( common.TargMassT2 + common.PZcms2 ) - common.PZcms;
    common.TPlusMax = common.SqrtS - common.ProjMassT;

    if ( G4UniformRand() <= 0.5 ) {  // Random choice projectile or target sampling
      if ( G4UniformRand() < theParameters->GetProbLogDistrPrD() ) {
        common.PMinusNew = ChooseP( common.PMinusMin, common.PMinusMax );
      } else {
        common.PMinusNew = ( common.PMinusMax - common.PMinusMin )*G4UniformRand() + common.PMinusMin;
      }
      if ( G4UniformRand() < theParameters->GetProbLogDistr() ) {
        common.TPlusNew = ChooseP( common.TPlusMin, common.TPlusMax );
      } else {
        common.TPlusNew = ( common.TPlusMax - common.TPlusMin )*G4UniformRand() + common.TPlusMin;
      }
    } else {
      if ( G4UniformRand() < theParameters->GetProbLogDistr() ) {
        common.TPlusNew = ChooseP( common.TPlusMin, common.TPlusMax );
      } else {
        common.TPlusNew = ( common.TPlusMax - common.TPlusMin )*G4UniformRand() + common.TPlusMin;
      }
      if ( G4UniformRand() < theParameters->GetProbLogDistrPrD() ) {
        common.PMinusNew = ChooseP( common.PMinusMin, common.PMinusMax );
      } else {
        common.PMinusNew = ( common.PMinusMax - common.PMinusMin )*G4UniformRand() + common.PMinusMin;
      }
    }

    common.Qminus = common.PMinusNew - common.Pprojectile.minus();
    common.Qplus = -( common.TPlusNew - common.Ptarget.plus() );
    common.Qmomentum.setPz( (common.Qplus - common.Qminus)/2.0 );
    common.Qmomentum.setE(  (common.Qplus + common.Qminus)/2.0 );
    #ifdef debugFTFexcitation
    G4cout <<"Sampled: Mpr, MdifPr, Mtr, MdifTr "<<G4endl
           << ( common.Pprojectile + common.Qmomentum ).mag() << " " 
           << common.ProjectileNonDiffStateMinMass << G4endl 
           << ( common.Ptarget - common.Qmomentum ).mag() << " "
           << common.TargetNonDiffStateMinMass << G4endl;
    #endif

  } while ( ( common.Pprojectile + common.Qmomentum ).mag2() < 
            common.ProjectileNonDiffStateMinMass2  ||  // No double Diffraction
            ( common.Ptarget - common.Qmomentum ).mag2() <  
            common.TargetNonDiffStateMinMass2 );  /* Loop checking, 10.08.2015, A.Ribon */

  projectile->SetStatus( 0 );
  target->SetStatus( 0 );
  return true;
}


//============================================================================

void G4DiffractiveExcitation::CreateStrings( G4VSplitableHadron* hadron, 
                                             G4bool isProjectile,
                                             G4ExcitedString*& FirstString, 
                                             G4ExcitedString*& SecondString,
                                             G4FTFParameters* theParameters ) const {

  //G4cout << "Create Strings SplitUp " << hadron << G4endl
  //       << "Defin " << hadron->GetDefinition() << G4endl
  //       << "Defin " << hadron->GetDefinition()->GetPDGEncoding() << G4endl;

  G4bool HadronIsString = hadron->IsSplit();
  if( ! HadronIsString ) hadron->SplitUp();

  G4Parton* start = hadron->GetNextParton();
  if ( start == nullptr ) { 
    G4cout << " G4FTFModel::String() Error: No start parton found" << G4endl;
    FirstString = 0; SecondString = 0;
    return;
  }

  G4Parton* end = hadron->GetNextParton();
  if ( end == nullptr ) { 
    G4cout << " G4FTFModel::String() Error: No end parton found" <<  G4endl;
    FirstString = 0; SecondString = 0;
    return;
  }

  //G4cout << start << " " << start->GetPDGcode() << " " << end << " " << end->GetPDGcode()
  //       << G4endl
  //       << "Create string " << start->GetPDGcode() << " " << end->GetPDGcode() << G4endl;

  if ( HadronIsString ) {
    if ( isProjectile ) {
      FirstString = new G4ExcitedString( end, start, +1 );
    } else {
      FirstString = new G4ExcitedString( end, start, -1 );
    }
    FirstString->SetTimeOfCreation( hadron->GetTimeOfCreation() );
    FirstString->SetPosition( hadron->GetPosition() );
    SecondString = 0;
    return;
  }

  G4LorentzVector Phadron = hadron->Get4Momentum();
  //G4cout << "String mom " << Phadron << G4endl;
  G4LorentzVector Pstart( 0.0, 0.0, 0.0, 0.0 );
  G4LorentzVector Pend( 0.0, 0.0, 0.0, 0.0 );
  G4LorentzVector Pkink( 0.0, 0.0, 0.0, 0.0 );
  G4LorentzVector PkinkQ1( 0.0, 0.0, 0.0, 0.0 );
  G4LorentzVector PkinkQ2( 0.0, 0.0, 0.0, 0.0 );

  G4int PDGcode_startQ = std::abs( start->GetDefinition()->GetPDGEncoding() );
  G4int PDGcode_endQ   = std::abs(   end->GetDefinition()->GetPDGEncoding() );
  //G4cout << "PDGcode_startQ " << PDGcode_startQ << " PDGcode_endQ   " << PDGcode_endQ << G4endl;

  G4double Wmin( 0.0 );
  if ( isProjectile ) {
    Wmin = theParameters->GetProjMinDiffMass();
  } else {
    Wmin = theParameters->GetTarMinDiffMass();
  }

  G4double W = hadron->Get4Momentum().mag();
  G4double W2 = W*W;
  G4double Pt( 0.0 ), x1( 0.0 ), x3( 0.0 );  // x2( 0.0 ) 
  G4bool Kink = false;

  if ( ! ( ( start->GetDefinition()->GetParticleSubType() == "di_quark"  &&
               end->GetDefinition()->GetParticleSubType() == "di_quark"  )  ||
           ( start->GetDefinition()->GetParticleSubType() == "quark"     &&
               end->GetDefinition()->GetParticleSubType() == "quark"     ) ) ) { 
    // Kinky strings are allowed only for qq-q strings;
    // Kinky strings are impossible for other systems (qq-qqbar, q-qbar)
    // according to the analysis of Pbar P interactions

    if ( W > Wmin ) {  // Kink is possible
      if ( hadron->GetStatus() == 0 ) {
        G4double Pt2kink = theParameters->GetPt2Kink(); // For non-diffractive
        if ( Pt2kink ) {
          Pt = std::sqrt( Pt2kink * ( G4Pow::GetInstance()->powA( W2/16.0/Pt2kink + 1.0, G4UniformRand() ) - 1.0 ) );
        } else {
          Pt = 0.0;
        }
      } else {
        Pt = 0.0;
      }

      if ( Pt > 500.0*MeV ) {
        G4double Ymax = G4Log( W/2.0/Pt + std::sqrt( W2/4.0/Pt/Pt - 1.0 ) );
        G4double Y = Ymax*( 1.0 - 2.0*G4UniformRand() );
        x1 = 1.0 - Pt/W * G4Exp( Y );
        x3 = 1.0 - Pt/W * G4Exp(-Y );
        //x2 = 2.0 - x1 - x3;

        G4double Mass_startQ = 650.0*MeV;
        if ( PDGcode_startQ <  3 ) Mass_startQ =  325.0*MeV;  // For constituent up or down quark
        if ( PDGcode_startQ == 3 ) Mass_startQ =  500.0*MeV;  // For constituent strange quark
        if ( PDGcode_startQ == 4 ) Mass_startQ = 1600.0*MeV;  // For constituent charm quark
        if ( PDGcode_startQ == 5 ) Mass_startQ = 4500.0*MeV;  // For constituent bottom quark
        G4double Mass_endQ = 650.0*MeV;
        if ( PDGcode_endQ <  3 ) Mass_endQ =  325.0*MeV;      // For constituent up or down quark
        if ( PDGcode_endQ == 3 ) Mass_endQ =  500.0*MeV;      // For constituent strange quark
        if ( PDGcode_endQ == 4 ) Mass_endQ = 1600.0*MeV;      // For constituent charm quark
        if ( PDGcode_endQ == 5 ) Mass_endQ = 4500.0*MeV;      // For constituent bottom quark

        G4double P2_1 = W2*x1*x1/4.0 - Mass_endQ*Mass_endQ;
        G4double P2_3 = W2*x3*x3/4.0 - Mass_startQ*Mass_startQ;
        G4double P2_2 = sqr( (2.0 - x1 - x3)*W/2.0 );
        if ( P2_1 <= 0.0  ||  P2_3 <= 0.0 ) { 
          Kink = false;
        } else {
          G4double P_1 = std::sqrt( P2_1 );
          G4double P_2 = std::sqrt( P2_2 );
          G4double P_3 = std::sqrt( P2_3 );
          G4double CosT12 = ( P2_3 - P2_1 - P2_2 ) / (2.0*P_1*P_2);
          G4double CosT13 = ( P2_2 - P2_1 - P2_3 ) / (2.0*P_1*P_3);

          if ( std::abs( CosT12 ) > 1.0  ||  std::abs( CosT13 ) > 1.0 ) {
            Kink = false;
          } else { 
            Kink = true; 
            Pt = P_2 * std::sqrt( 1.0 - CosT12*CosT12 );  // because system was rotated
            Pstart.setPx( -Pt ); Pstart.setPy( 0.0 ); Pstart.setPz( P_3*CosT13 ); 
            Pend.setPx(   0.0 ); Pend.setPy(   0.0 ); Pend.setPz(          P_1 ); 
            Pkink.setPx(   Pt ); Pkink.setPy(  0.0 ); Pkink.setPz(  P_2*CosT12 );
            Pstart.setE( x3*W/2.0 );                
            Pkink.setE( Pkink.vect().mag() );
            Pend.setE( x1*W/2.0 );

            G4double XkQ = GetQuarkFractionOfKink( 0.0, 1.0 );
            if ( Pkink.getZ() > 0.0 ) {
              if ( XkQ > 0.5 ) {
                PkinkQ1 = XkQ*Pkink;
              } else {
                PkinkQ1 = (1.0 - XkQ)*Pkink;
              }
            } else {
              if ( XkQ > 0.5 ) {
                PkinkQ1 = (1.0 - XkQ)*Pkink;
              } else {
                PkinkQ1 = XkQ*Pkink;
              }
            }

            PkinkQ2 = Pkink - PkinkQ1;
            // Minimizing Pt1^2+Pt3^2
            G4double Cos2Psi = ( sqr(x1) - sqr(x3) + 2.0*sqr( x3*CosT13 ) ) /
                               std::sqrt( sqr( sqr(x1) - sqr(x3) ) + sqr( 2.0*x1*x3*CosT13 ) );
            G4double Psi = std::acos( Cos2Psi );

            G4LorentzRotation Rotate;
            if ( isProjectile ) {
              Rotate.rotateY( Psi );
            } else {
              Rotate.rotateY( pi + Psi );
            }                   
            Rotate.rotateZ( twopi * G4UniformRand() );
            Pstart *= Rotate;
            Pkink *= Rotate;
            PkinkQ1 *= Rotate;
            PkinkQ2 *= Rotate;
            Pend *= Rotate;
          }
        }  // End of if ( P2_1 <= 0.0  ||  P2_3 <= 0.0 )
      }  // End of if ( Pt > 500.0*MeV )
    } // End of if ( W > Wmin ) : check for a kink
  }  // end of qq-q string selection

  if ( Kink ) {  // Kink is possible

    //G4cout << "Kink is sampled!" << G4endl;
    std::vector< G4double > QuarkProbabilitiesAtGluonSplitUp = 
        theParameters->GetQuarkProbabilitiesAtGluonSplitUp();

    G4int QuarkInGluon( 1 ); G4double Ksi = G4UniformRand();
    for ( unsigned int Iq = 0; Iq < 3; Iq++ ) {
      //G4cout << "Iq " << Iq << G4endl;
      if ( Ksi > QuarkProbabilitiesAtGluonSplitUp[Iq] ) QuarkInGluon++;
    }
    //G4cout << "Last Iq " << QuarkInGluon << G4endl;
    G4Parton* Gquark = new G4Parton( QuarkInGluon );
    G4Parton* Ganti_quark = new G4Parton( -QuarkInGluon );
    //G4cout << "Lorentz " << G4endl;

    G4LorentzRotation toCMS( -1 * Phadron.boostVector() );
    G4LorentzRotation toLab( toCMS.inverse() );
    //G4cout << "Pstart " << Pstart << G4endl;
    //G4cout << "Pend   " << Pend << G4endl;
    //G4cout << "Kink1  " <<PkinkQ1<<G4endl;
    //G4cout << "Kink2  " <<PkinkQ2<<G4endl;
    //G4cout << "Pstart " << Pstart << G4endl<<G4endl;

    Pstart.transform( toLab );  start->Set4Momentum( Pstart );
    PkinkQ1.transform( toLab );
    PkinkQ2.transform( toLab );
    Pend.transform( toLab );    end->Set4Momentum( Pend );
    //G4cout << "Pstart " << Pstart << G4endl;
    //G4cout << "Pend   " << Pend << G4endl;
    //G4cout << "Defin " << hadron->GetDefinition()<< G4endl;
    //G4cout << "Defin " << hadron->GetDefinition()->GetPDGEncoding()<< G4endl;

    //G4int absPDGcode = std::abs( hadron->GetDefinition()->GetPDGEncoding() );
    G4int absPDGcode = 1500;
    if ( start->GetDefinition()->GetParticleSubType() == "quark"  &&
         end->GetDefinition()->GetParticleSubType() == "quark" ) {
      absPDGcode = 110;
    }
    //G4cout << "absPDGcode " << absPDGcode << G4endl;

    if ( absPDGcode < 1000 ) {  // meson
      if ( isProjectile ) { // Projectile
        if ( end->GetDefinition()->GetPDGEncoding() > 0 ) {  // A quark on the end
          FirstString  = new G4ExcitedString( end   , Ganti_quark, +1 );
          SecondString = new G4ExcitedString( Gquark, start      , +1 );
          Ganti_quark->Set4Momentum( PkinkQ1 );
          Gquark->Set4Momentum( PkinkQ2 );
        } else {  // Anti_Quark on the end
          FirstString  = new G4ExcitedString( end        , Gquark, +1 );
          SecondString = new G4ExcitedString( Ganti_quark, start , +1 );
          Gquark->Set4Momentum( PkinkQ1 );
          Ganti_quark->Set4Momentum( PkinkQ2 );
        }
      } else {  // Target
        if ( end->GetDefinition()->GetPDGEncoding() > 0 ) { // A quark on the end
          FirstString  = new G4ExcitedString( Ganti_quark, end   , -1 );
          SecondString = new G4ExcitedString( start      , Gquark, -1 );
          Ganti_quark->Set4Momentum( PkinkQ2 );
          Gquark->Set4Momentum( PkinkQ1 );
        } else {  // Anti_Quark on the end
          FirstString  = new G4ExcitedString( Gquark, end        , -1 );
          SecondString = new G4ExcitedString( start , Ganti_quark, -1 );
          Gquark->Set4Momentum( PkinkQ2 );
          Ganti_quark->Set4Momentum( PkinkQ1 );
        }
      }
    } else {  // Baryon/AntiBaryon
      if ( isProjectile ) {  // Projectile
        if ( end->GetDefinition()->GetParticleType() == "diquarks"  &&
             end->GetDefinition()->GetPDGEncoding() > 0 ) {  // DiQuark on the end
          FirstString  = new G4ExcitedString( end        , Gquark, +1 );
          SecondString = new G4ExcitedString( Ganti_quark, start , +1 );
          Gquark->Set4Momentum( PkinkQ1 );
          Ganti_quark->Set4Momentum( PkinkQ2 );
        } else {                            // Anti_DiQuark on the end or quark
          FirstString  = new G4ExcitedString( end   , Ganti_quark, +1 );
          SecondString = new G4ExcitedString( Gquark, start      , +1 );
          Ganti_quark->Set4Momentum( PkinkQ1 );
          Gquark->Set4Momentum( PkinkQ2 );
        }
      } else {  // Target
        if ( end->GetDefinition()->GetParticleType() == "diquarks"  &&
             end->GetDefinition()->GetPDGEncoding() > 0 ) {  // DiQuark on the end
          Gquark->Set4Momentum( PkinkQ1 );
          Ganti_quark->Set4Momentum( PkinkQ2 );
          FirstString  = new G4ExcitedString(         end, Gquark, -1 );
          SecondString = new G4ExcitedString( Ganti_quark,  start, -1 );
        } else {  // Anti_DiQuark on the end or Q
          FirstString  = new G4ExcitedString( Ganti_quark, end   , -1 );
          SecondString = new G4ExcitedString( start      , Gquark, -1 );
          Gquark->Set4Momentum( PkinkQ2 );
          Ganti_quark->Set4Momentum( PkinkQ1 );
        }
      }
    }

    FirstString->SetTimeOfCreation( hadron->GetTimeOfCreation() );
    FirstString->SetPosition( hadron->GetPosition() );
    SecondString->SetTimeOfCreation( hadron->GetTimeOfCreation() );
    SecondString->SetPosition( hadron->GetPosition() );
  
  } else {  // Kink is impossible

    if ( isProjectile ) {
      FirstString = new G4ExcitedString( end, start, +1 );
    } else {
      FirstString = new G4ExcitedString( end, start, -1 );
    }
    FirstString->SetTimeOfCreation( hadron->GetTimeOfCreation() );
    FirstString->SetPosition( hadron->GetPosition() );
    SecondString = 0;
    // momenta of string ends
    G4LorentzVector HadronMom = hadron->Get4Momentum();
    G4LorentzVector Pstart1 = G4LorentzVector( HadronMom.px()/2.0, HadronMom.py()/2.0, 0.0, 0.0 );  //    Quark momentum
    G4LorentzVector   Pend1 = G4LorentzVector( HadronMom.px()/2.0, HadronMom.py()/2.0, 0.0, 0.0 );  // Di-quark momentum
    G4double Pz  = HadronMom.pz();
    G4double Eh  = HadronMom.e();
    G4double Pt2 = sqr( HadronMom.px() ) + sqr( HadronMom.py() );
    G4double Mt2 = HadronMom.mt2();
    G4double Exp = std::sqrt( sqr(Pz) + ( sqr(Mt2) - 4.0*sqr(Eh)*Pt2/4.0 )/Mt2 )/2.0;
    G4double Pzq  = Pz/2.0 - Exp;                      Pstart1.setZ( Pzq );
    G4double Eq   = std::sqrt( sqr(Pzq) + Pt2/4.0 );   Pstart1.setE( Eq );
    G4double Pzqq = Pz/2.0 + Exp;                      Pend1.setZ(Pzqq);
    G4double Eqq  = std::sqrt( sqr(Pzqq) + Pt2/4.0 );  Pend1.setE(Eqq);
    start->Set4Momentum( Pstart1 );
    end->Set4Momentum( Pend1 );
    Pstart = Pstart1;  Pend = Pend1;

  }  // End of "if (Kink)"

  //G4cout << "Quarks in the string at creation" << FirstString->GetRightParton()->GetPDGcode()
  //       << " " << FirstString->GetLeftParton()->GetPDGcode() << G4endl
  //       << FirstString << " " << SecondString << G4endl;

  #ifdef G4_FTFDEBUG
  G4cout << " generated string flavors          " << start->GetPDGcode() << " / " 
         << end->GetPDGcode() << G4endl << " generated string momenta:   quark " 
         << start->Get4Momentum() << "mass : " << start->Get4Momentum().mag() << G4endl
         << " generated string momenta: Diquark " << end->Get4Momentum() << "mass : " 
         << end->Get4Momentum().mag() << G4endl << " sum of ends                       "
         << Pstart + Pend << G4endl << " Original                          " 
         << hadron->Get4Momentum() << " "<<hadron->Get4Momentum().mag() << G4endl;
  #endif

  return;
}


//============================================================================

G4double G4DiffractiveExcitation::ChooseP( G4double Pmin, G4double Pmax ) const {
  // Choose an x between Xmin and Xmax with P(x) ~ 1/x . 
  // To be improved...
  G4double range = Pmax - Pmin;                    
  if ( Pmin <= 0.0 || range <= 0.0 ) {
    G4cout << " Pmin, range : " << Pmin << " , " << range << G4endl;
    throw G4HadronicException( __FILE__, __LINE__,
                               "G4DiffractiveExcitation::ChooseP : Invalid arguments " );
  }
  G4double P = Pmin * G4Pow::GetInstance()->powA( Pmax/Pmin, G4UniformRand() ); 
  //G4double P = (Pmax - Pmin) * G4UniformRand() + Pmin;
  return P;
}


//============================================================================

G4ThreeVector G4DiffractiveExcitation::GaussianPt( G4double AveragePt2, G4double maxPtSquare ) const {
  //  @@ this method is used in FTFModel as well. Should go somewhere common!
  G4double Pt2( 0.0 );
  if ( AveragePt2 <= 0.0 ) {
    Pt2 = 0.0;
  } else {
    Pt2 = -AveragePt2 * G4Log( 1.0 + G4UniformRand() * 
                                       ( G4Exp( -maxPtSquare/AveragePt2 ) - 1.0 ) );
  }
  G4double Pt = std::sqrt( Pt2 );
  G4double phi = G4UniformRand() * twopi;
  return G4ThreeVector( Pt * std::cos( phi ), Pt * std::sin( phi ), 0.0 );
}


//============================================================================

G4double G4DiffractiveExcitation::GetQuarkFractionOfKink( G4double zmin, G4double zmax ) const {
  G4double z, yf;
  const G4int maxNumberOfLoops = 10000;
  G4int loopCounter = 0;
  do {
    z = zmin + G4UniformRand() * (zmax - zmin);
    yf = z*z + sqr(1.0 - z);
  } while ( ( G4UniformRand() > yf ) && 
            ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */
  if ( loopCounter >= maxNumberOfLoops ) {
    z = 0.5*(zmin + zmax);  // Just something acceptable, without any physics consideration.
  }
  return z;
}


//============================================================================

void G4DiffractiveExcitation::UnpackMeson( const G4int IdPDG, G4int& Q1, G4int& Q2 ) const {
  G4int absIdPDG = std::abs( IdPDG );
  if ( ! ( absIdPDG == 111 || absIdPDG == 221 || absIdPDG == 331 ||     // Pi0  ,  Eta ,   Eta'
           absIdPDG == 441 || absIdPDG == 443 || absIdPDG == 553 ) ) {  // Etac ,  J/psi , Upsilon 
    // All other projectile mesons (including charmed and bottom ones)
    Q1 =  absIdPDG / 100;
    Q2 = (absIdPDG % 100) / 10;
    G4int anti = 1 - 2 * ( std::max( Q1, Q2 ) % 2 );
    if ( IdPDG < 0 ) anti *= -1;
    Q1 *= anti;
    Q2 *= -1 * anti;
  } else {
    if ( absIdPDG == 441 || absIdPDG == 443 ) {  // Etac , J/psi
      Q1 =  4; Q2 = -4;
    } else if ( absIdPDG == 553 ) {              // Upsilon
      Q1 =  5; Q2 = -5;
    } else {                                     // Pi0 , Eta , Eta'
      if ( G4UniformRand() < 0.5 ) {
	Q1 = 1; Q2 = -1;
      } else {
	Q1 = 2; Q2 = -2;
      }
    }
  }
  return;
}


//============================================================================

void G4DiffractiveExcitation::UnpackBaryon( G4int IdPDG, 
                                            G4int& Q1, G4int& Q2, G4int& Q3 ) const {
  Q1 = IdPDG          / 1000;
  Q2 = (IdPDG % 1000) / 100;
  Q3 = (IdPDG % 100)  / 10;
  return;
}


//============================================================================

G4int G4DiffractiveExcitation::NewNucleonId( G4int Q1, G4int Q2, G4int Q3 ) const {
  // Order the three integers in such a way that Q1 >= Q2 >= Q3
  G4int TmpQ( 0 );
  if ( Q3 > Q2 ) {
    TmpQ = Q2;
    Q2 = Q3;
    Q3 = TmpQ;
  } else if ( Q3 > Q1 ) {
    TmpQ = Q1;
    Q1 = Q3;
    Q3 = TmpQ;
  }
  if ( Q2 > Q1 ) {
    TmpQ = Q1;
    Q1 = Q2;
    Q2 = TmpQ;
  }
  // By now Q1 >= Q2 >= Q3
  G4int NewCode = Q1*1000 + Q2*100 + Q3*10 + 2;
  return NewCode;
}


//============================================================================

G4DiffractiveExcitation::G4DiffractiveExcitation( const G4DiffractiveExcitation& ) {
  throw G4HadronicException( __FILE__, __LINE__, 
                             "G4DiffractiveExcitation copy constructor not meant to be called" );
}


//============================================================================

const G4DiffractiveExcitation & G4DiffractiveExcitation::operator=( const G4DiffractiveExcitation& ) {
  throw G4HadronicException( __FILE__, __LINE__, 
                             "G4DiffractiveExcitation = operator not meant to be called" );
  return *this;
}


//============================================================================

G4bool G4DiffractiveExcitation::operator==( const G4DiffractiveExcitation& ) const {
  throw G4HadronicException( __FILE__, __LINE__, 
                             "G4DiffractiveExcitation == operator not meant to be called" );
}


//============================================================================

G4bool G4DiffractiveExcitation::operator!= ( const G4DiffractiveExcitation& ) const {
  throw G4HadronicException( __FILE__, __LINE__, 
                             "G4DiffractiveExcitation != operator not meant to be called" );
}

