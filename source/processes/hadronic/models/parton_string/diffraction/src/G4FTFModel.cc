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
// $Id: G4FTFModel.cc 81880 2014-06-06 12:50:36Z gcosmo $
// GEANT4 tag $Name:  $
//

// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4FTFModel ----------------
//             by Gunter Folger, May 1998.
//       class implementing the excitation in the FTF Parton String Model
//
//                Vladimir Uzhinsky, November - December 2012
//       simulation of nucleus-nucleus interactions was implemented.
// ------------------------------------------------------------

#include <utility> 

#include "G4FTFModel.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4FTFParameters.hh"
#include "G4FTFParticipants.hh"
#include "G4DiffractiveSplitableHadron.hh"
#include "G4InteractionContent.hh"
#include "G4LorentzRotation.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"


//============================================================================

//#define debugFTFmodel
//#define debugReggeonCascade
//#define debugPutOnMassShell
//#define debugAdjust
//#define debugBuildString


//============================================================================

G4FTFModel::G4FTFModel( const G4String& modelName ) : 
  G4VPartonStringModel( modelName ),
  theExcitation( new G4DiffractiveExcitation() ),
  theElastic( new G4ElasticHNScattering() ),
  theAnnihilation( new G4FTFAnnihilation() )
{
  G4VPartonStringModel::SetThisPointer( this );
  theParameters = 0;
  NumberOfInvolvedNucleonsOfTarget = 0;
  NumberOfInvolvedNucleonsOfProjectile= 0;

  LowEnergyLimit = 2000.0*MeV;
  HighEnergyInter = true;

  G4LorentzVector tmp( 0.0, 0.0, 0.0, 0.0 );
  ProjectileResidual4Momentum        = tmp;
  ProjectileResidualMassNumber       = 0;
  ProjectileResidualCharge           = 0;
  ProjectileResidualExcitationEnergy = 0.0;

  TargetResidual4Momentum            = tmp;
  TargetResidualMassNumber           = 0;
  TargetResidualCharge               = 0;
  TargetResidualExcitationEnergy     = 0.0;

  SetEnergyMomentumCheckLevels( 2.0*perCent, 150.0*MeV );
}


//============================================================================

struct DeleteVSplitableHadron { void operator()( G4VSplitableHadron* aH ){ delete aH; } };


//============================================================================

G4FTFModel::~G4FTFModel() {
   // Because FTF model can be called for various particles
   // theParameters must be erased at the end of each call.
   // Thus the delete is also in G4FTFModel::GetStrings() method.
   if ( theParameters   != 0 ) delete theParameters; 
   if ( theExcitation   != 0 ) delete theExcitation;
   if ( theElastic      != 0 ) delete theElastic; 
   if ( theAnnihilation != 0 ) delete theAnnihilation;

   // Erasing of strings created at annihilation.
   if ( theAdditionalString.size() != 0 ) {
     std::for_each( theAdditionalString.begin(), theAdditionalString.end(), 
                    DeleteVSplitableHadron() );
   }
   theAdditionalString.clear();

   // Erasing of target involved nucleons.
   if ( NumberOfInvolvedNucleonsOfTarget != 0 ) {
     for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; i++ ) {
       G4VSplitableHadron* aNucleon = TheInvolvedNucleonsOfTarget[i]->GetSplitableHadron();
       if ( aNucleon ) delete aNucleon;
     }
   }

   // Erasing of projectile involved nucleons.
   if ( NumberOfInvolvedNucleonsOfProjectile != 0 ) {
     for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfProjectile; i++ ) {
       G4VSplitableHadron* aNucleon = TheInvolvedNucleonsOfProjectile[i]->GetSplitableHadron();
       if ( aNucleon ) delete aNucleon;
     }
   }
}


//============================================================================

void G4FTFModel::Init( const G4Nucleus& aNucleus, const G4DynamicParticle& aProjectile ) {

  theProjectile = aProjectile;  

  G4double PlabPerParticle( 0.0 );  // Laboratory momentum Pz per particle/nucleon

  #ifdef debugFTFmodel
  G4cout << "FTF init Proj Name " << theProjectile.GetDefinition()->GetParticleName() << G4endl
         << "FTF init Proj Mass " << theProjectile.GetMass() 
         << " " << theProjectile.GetMomentum() << G4endl
         << "FTF init Proj B Q  " << theProjectile.GetDefinition()->GetBaryonNumber()
         << " " << (G4int) theProjectile.GetDefinition()->GetPDGCharge() << G4endl 
         << "FTF init Target A Z " << aNucleus.GetA_asInt() 
         << " " << aNucleus.GetZ_asInt() << G4endl;
  #endif

  theParticipants.SetProjectileNucleus( 0 );

  G4LorentzVector tmp( 0.0, 0.0, 0.0, 0.0 );
  ProjectileResidualMassNumber       = 0;
  ProjectileResidualCharge           = 0;
  ProjectileResidualExcitationEnergy = 0.0;
  ProjectileResidual4Momentum        = tmp;

  TargetResidualMassNumber       = aNucleus.GetA_asInt();
  TargetResidualCharge           = aNucleus.GetZ_asInt();
  TargetResidualExcitationEnergy = 0.0;
  TargetResidual4Momentum        = tmp;
  G4double TargetResidualMass = G4ParticleTable::GetParticleTable()->GetIonTable()
                                ->GetIonMass( TargetResidualCharge, TargetResidualMassNumber );

  TargetResidual4Momentum.setE( TargetResidualMass );

  if ( std::abs( theProjectile.GetDefinition()->GetBaryonNumber() ) <= 1 ) { 
    // Projectile is a hadron : meson or baryon
    PlabPerParticle = theProjectile.GetMomentum().z();
    ProjectileResidualMassNumber = std::abs( theProjectile.GetDefinition()->GetBaryonNumber() );
    ProjectileResidualCharge = G4int( theProjectile.GetDefinition()->GetPDGCharge() );
    ProjectileResidualExcitationEnergy = 0.0;
    // G4double ProjectileResidualMass = theProjectile.GetMass();
    ProjectileResidual4Momentum.setVect( theProjectile.GetMomentum() );
    ProjectileResidual4Momentum.setE( theProjectile.GetTotalEnergy() );
    if ( PlabPerParticle < LowEnergyLimit ) {
      HighEnergyInter = false;
    } else {
      HighEnergyInter = true;
    }
  } else {
    if ( theProjectile.GetDefinition()->GetBaryonNumber() > 1 ) { 
      // Projectile is a nucleus
      theParticipants.InitProjectileNucleus(theProjectile.GetDefinition()->GetBaryonNumber(),
                                            G4int(theProjectile.GetDefinition()->GetPDGCharge()));
      ProjectileResidualMassNumber = theProjectile.GetDefinition()->GetBaryonNumber();
      ProjectileResidualCharge = G4int( theProjectile.GetDefinition()->GetPDGCharge() );
      PlabPerParticle = theProjectile.GetMomentum().z() / 
                        theProjectile.GetDefinition()->GetBaryonNumber();
      if ( PlabPerParticle < LowEnergyLimit ) {
        HighEnergyInter = false;
      } else {
        HighEnergyInter = true;
      }
    } else if ( theProjectile.GetDefinition()->GetBaryonNumber() < -1 ) { 
      // Projectile is an anti-nucleus
      theParticipants.InitProjectileNucleus( 
          std::abs( theProjectile.GetDefinition()->GetBaryonNumber() ),
          std::abs( G4int( theProjectile.GetDefinition()->GetPDGCharge() ) ) );
      theParticipants.theProjectileNucleus->StartLoop();
      G4Nucleon* aNucleon;
      while ( ( aNucleon = theParticipants.theProjectileNucleus->GetNextNucleon() ) ) {
        if ( aNucleon->GetDefinition() == G4Proton::Proton() ) {
          aNucleon->SetParticleType( G4AntiProton::AntiProton() ); 
        } else if ( aNucleon->GetDefinition() == G4Neutron::Neutron() ) {
          aNucleon->SetParticleType( G4AntiNeutron::AntiNeutron() );
        } 
      }
      ProjectileResidualMassNumber = std::abs( theProjectile.GetDefinition()->GetBaryonNumber() );
      ProjectileResidualCharge = std::abs( G4int(theProjectile.GetDefinition()->GetPDGCharge()) );
      PlabPerParticle = theProjectile.GetMomentum().z() /
                        std::abs( theProjectile.GetDefinition()->GetBaryonNumber() );
      if ( PlabPerParticle < LowEnergyLimit ) {
        HighEnergyInter = false;
      } else {
        HighEnergyInter = true;
      }
    }
    G4ThreeVector BoostVector = theProjectile.GetMomentum() / theProjectile.GetTotalEnergy();
    theParticipants.theProjectileNucleus->DoLorentzBoost( BoostVector );
    theParticipants.theProjectileNucleus->DoLorentzContraction( BoostVector );
    ProjectileResidualExcitationEnergy = 0.0;
    //G4double ProjectileResidualMass = theProjectile.GetMass();
    ProjectileResidual4Momentum.setVect( theProjectile.GetMomentum() );
    ProjectileResidual4Momentum.setE( theProjectile.GetTotalEnergy() );
  }

  // Init target nucleus
  theParticipants.Init( aNucleus.GetA_asInt(), aNucleus.GetZ_asInt() );

  if ( theParameters != 0 ) delete theParameters;
  theParameters = new G4FTFParameters( theProjectile.GetDefinition(), aNucleus.GetA_asInt(),
                                       aNucleus.GetZ_asInt(), PlabPerParticle );

  if ( theAdditionalString.size() != 0 ) {
    std::for_each( theAdditionalString.begin(), theAdditionalString.end(), 
                   DeleteVSplitableHadron() );
  }
  theAdditionalString.clear();

  #ifdef debugFTFmodel
  G4cout << "FTF end of Init" << G4endl << G4endl;
  #endif

}


//============================================================================

G4ExcitedStringVector* G4FTFModel::GetStrings() { 

  #ifdef debugFTFmodel
  G4cout << "G4FTFModel::GetStrings() " << G4endl;
  #endif

  G4ExcitedStringVector* theStrings( 0 );
  theParticipants.GetList( theProjectile, theParameters );
  StoreInvolvedNucleon(); 

  G4bool Success( true );

  if ( HighEnergyInter ) {
    ReggeonCascade(); 
 
    #ifdef debugFTFmodel
    G4cout << "FTF PutOnMassShell " << G4endl;
    #endif

    Success = PutOnMassShell(); 

    #ifdef debugFTFmodel
    G4cout << "FTF PutOnMassShell Success?  " << Success << G4endl;
    #endif

  }

  #ifdef debugFTFmodel
  G4cout << "FTF ExciteParticipants " << G4endl;
  #endif

  if ( Success ) Success = ExciteParticipants();

  #ifdef debugFTFmodel
  G4cout << "FTF ExciteParticipants Success? " << Success << G4endl;
  #endif

  if ( Success ) {       

    #ifdef debugFTFmodel
    G4cout << "FTF BuildStrings ";
    #endif

    theStrings = BuildStrings();

    #ifdef debugFTFmodel
    G4cout << "FTF BuildStrings " << theStrings << " OK" << G4endl
           << "FTF GetResiduals of Nuclei " << G4endl;
    #endif

    GetResiduals();

    if ( theParameters != 0 ) {
      delete theParameters;
      theParameters = 0;
    }
  } else if ( ! GetProjectileNucleus() ) {
    // Erase the hadron projectile
    std::vector< G4VSplitableHadron* > primaries;
    theParticipants.StartLoop();
    while ( theParticipants.Next() ) {
      const G4InteractionContent& interaction = theParticipants.GetInteraction();
      // Do not allow for duplicates
      if ( primaries.end() == 
           std::find( primaries.begin(), primaries.end(), interaction.GetProjectile() ) ) {
        primaries.push_back( interaction.GetProjectile() );
      }
    }
    std::for_each( primaries.begin(), primaries.end(), DeleteVSplitableHadron() );
    primaries.clear();
  }

  // Cleaning of the memory
  G4VSplitableHadron* aNucleon = 0;

  // Erase the projectile nucleons
  for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfProjectile; i++ ) {
    aNucleon = TheInvolvedNucleonsOfProjectile[i]->GetSplitableHadron();
    if ( aNucleon ) delete aNucleon;
  } 
  NumberOfInvolvedNucleonsOfProjectile = 0;

  // Erase the target nucleons
  for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; i++ ) {
    aNucleon = TheInvolvedNucleonsOfTarget[i]->GetSplitableHadron();
    if ( aNucleon ) delete aNucleon;
  } 
  NumberOfInvolvedNucleonsOfTarget = 0;

  #ifdef debugFTFmodel
  G4cout << "End of FTF. Go to fragmentation" << G4endl
         << "To continue - enter 1, to stop - ^C" << G4endl;
  G4int Uzhi; G4cin >> Uzhi;
  #endif

  return theStrings;
}


//============================================================================

void G4FTFModel::StoreInvolvedNucleon() {
  //To store nucleons involved in the interaction

  NumberOfInvolvedNucleonsOfTarget = 0;

  G4V3DNucleus* theTargetNucleus = GetTargetNucleus();
  theTargetNucleus->StartLoop();

  G4Nucleon* aNucleon;
  while ( ( aNucleon = theTargetNucleus->GetNextNucleon() ) ) {
    if ( aNucleon->AreYouHit() ) {
      TheInvolvedNucleonsOfTarget[NumberOfInvolvedNucleonsOfTarget] = aNucleon;
      NumberOfInvolvedNucleonsOfTarget++;
    }
  }

  #ifdef debugFTFmodel
  G4cout << "G4FTFModel::StoreInvolvedNucleon -------------" << G4endl;
  G4cout << "NumberOfInvolvedNucleonsOfTarget " << NumberOfInvolvedNucleonsOfTarget
         << G4endl << G4endl;
  #endif

  if ( ! GetProjectileNucleus() ) return;  // The projectile is a hadron

  // The projectile is a nucleus or an anti-nucleus.

  NumberOfInvolvedNucleonsOfProjectile = 0;

  G4V3DNucleus* theProjectileNucleus = GetProjectileNucleus();
  theProjectileNucleus->StartLoop();

  G4Nucleon* aProjectileNucleon;
  while ( ( aProjectileNucleon = theProjectileNucleus->GetNextNucleon() ) ) {
    if ( aProjectileNucleon->AreYouHit() ) {  
      // Projectile nucleon was involved in the interaction.
      TheInvolvedNucleonsOfProjectile[NumberOfInvolvedNucleonsOfProjectile] = aProjectileNucleon;
      NumberOfInvolvedNucleonsOfProjectile++;
    }
  } 

  #ifdef debugFTFmodel
  G4cout << "NumberOfInvolvedNucleonsOfProjectile " << NumberOfInvolvedNucleonsOfProjectile
         << G4endl << G4endl;
  #endif

  return;
}                        


//============================================================================

void G4FTFModel::ReggeonCascade() { 
  // Implementation of the reggeon theory inspired model

  G4double ExcitationE = theParameters->GetExcitationEnergyPerWoundedNucleon();

  #ifdef debugReggeonCascade
  G4cout << "G4FTFModel::ReggeonCascade -----------" << G4endl
         << "theProjectile.GetTotalMomentum() " << theProjectile.GetTotalMomentum() << G4endl
         << "theProjectile.GetTotalEnergy() " << theProjectile.GetTotalEnergy() << G4endl
         << "ExcitationE/WN " << ExcitationE << G4endl;
  #endif

  G4int InitNINt = NumberOfInvolvedNucleonsOfTarget;

  // Reggeon cascading in target nucleus
  for ( G4int InvTN = 0; InvTN < InitNINt; InvTN++ ) { 
    G4Nucleon* aTargetNucleon = TheInvolvedNucleonsOfTarget[ InvTN ];
    aTargetNucleon->SetBindingEnergy( ExcitationE );

    G4double CreationTime = aTargetNucleon->GetSplitableHadron()->GetTimeOfCreation();

    G4double XofWoundedNucleon = aTargetNucleon->GetPosition().x();
    G4double YofWoundedNucleon = aTargetNucleon->GetPosition().y();
           
    G4V3DNucleus* theTargetNucleus = GetTargetNucleus();
    theTargetNucleus->StartLoop();

    G4Nucleon* Neighbour(0);
    while ( ( Neighbour = theTargetNucleus->GetNextNucleon() ) ) {
      if ( ! Neighbour->AreYouHit() ) {
        G4double impact2 = sqr( XofWoundedNucleon - Neighbour->GetPosition().x() ) +
                           sqr( YofWoundedNucleon - Neighbour->GetPosition().y() );

        if ( G4UniformRand() < theParameters->GetCofNuclearDestruction() *
                               std::exp( -impact2 / theParameters->GetR2ofNuclearDestruction() )
           ) {  
          // The neighbour nucleon is involved in the reggeon cascade
          TheInvolvedNucleonsOfTarget[ NumberOfInvolvedNucleonsOfTarget ] = Neighbour;
          NumberOfInvolvedNucleonsOfTarget++;

          G4VSplitableHadron* targetSplitable; 
          targetSplitable = new G4DiffractiveSplitableHadron( *Neighbour ); 

          Neighbour->Hit( targetSplitable );
          targetSplitable->SetTimeOfCreation( CreationTime );
          targetSplitable->SetStatus( 2 );     
        }
      }
    }
  }

  #ifdef debugReggeonCascade
  G4cout << "Final NumberOfInvolvedNucleonsOfTarget " 
         << NumberOfInvolvedNucleonsOfTarget << G4endl << G4endl;
  #endif

  if ( ! GetProjectileNucleus() ) return;

  // Nucleus-Nucleus Interaction : Destruction of Projectile
  for ( G4int InvPN = 0; InvPN < NumberOfInvolvedNucleonsOfProjectile; InvPN++ ) { 
    G4Nucleon* aProjectileNucleon = TheInvolvedNucleonsOfProjectile[ InvPN ];
    aProjectileNucleon->SetBindingEnergy( ExcitationE );

    G4double CreationTime = aProjectileNucleon->GetSplitableHadron()->GetTimeOfCreation();

    G4double XofWoundedNucleon = aProjectileNucleon->GetPosition().x();
    G4double YofWoundedNucleon = aProjectileNucleon->GetPosition().y();
           
    G4V3DNucleus* theProjectileNucleus = GetProjectileNucleus();
    theProjectileNucleus->StartLoop();

    G4Nucleon* Neighbour( 0 );
    while ( ( Neighbour = theProjectileNucleus->GetNextNucleon() ) ) {
      if ( ! Neighbour->AreYouHit() ) {
        G4double impact2= sqr( XofWoundedNucleon - Neighbour->GetPosition().x() ) +
                          sqr( YofWoundedNucleon - Neighbour->GetPosition().y() );

        if ( G4UniformRand() < theParameters->GetCofNuclearDestruction() *
                               std::exp( -impact2 / theParameters->GetR2ofNuclearDestruction() )
           ) {
          // The neighbour nucleon is involved in the reggeon cascade
          TheInvolvedNucleonsOfProjectile[ NumberOfInvolvedNucleonsOfProjectile ] = Neighbour;
          NumberOfInvolvedNucleonsOfProjectile++;

          G4VSplitableHadron* projectileSplitable; 
          projectileSplitable = new G4DiffractiveSplitableHadron( *Neighbour ); 

          Neighbour->Hit( projectileSplitable );
          projectileSplitable->SetTimeOfCreation( CreationTime );
          projectileSplitable->SetStatus( 2 );     
        }
      }
    }
  }

  #ifdef debugReggeonCascade
  G4cout << "NumberOfInvolvedNucleonsOfProjectile "
         << NumberOfInvolvedNucleonsOfProjectile << G4endl << G4endl;
  #endif
}   


//============================================================================

G4bool G4FTFModel::PutOnMassShell() {

  #ifdef debugPutOnMassShell
  G4cout << "PutOnMassShell start " << G4endl;
  #endif

  if ( ! GetProjectileNucleus() ) {  // The projectile is hadron or anti-baryon

    G4LorentzVector Pprojectile( theProjectile.GetMomentum(), theProjectile.GetTotalEnergy() );
    if ( Pprojectile.z() < 0.0 ) {
      return false;
    }
    G4double Mprojectile  = Pprojectile.mag();
    G4double M2projectile = Pprojectile.mag2();
    G4LorentzVector Psum  = Pprojectile;
    G4double SumMasses = Mprojectile + 20.0*MeV;  // Separation energy for nuclear nucleon
    // if ( ProjectileIsAntiBaryon ) SumMasses = Mprojectile;

    // Target nucleus 
    G4LorentzVector Ptarget( 0.0, 0.0, 0.0, 0.0 );
    G4LorentzVector PtargetResidual( 0.0, 0.0, 0.0, 0.0 );
    G4double ExcitationEnergyPerWoundedNucleon = 
        theParameters->GetExcitationEnergyPerWoundedNucleon();

    #ifdef debugPutOnMassShell
    G4cout << "ExcitationEnergyPerWoundedNucleon " << ExcitationEnergyPerWoundedNucleon << G4endl;
    #endif

    G4V3DNucleus* theNucleus = GetTargetNucleus();
    theNucleus->StartLoop();
    G4Nucleon* aNucleon( 0 );
    while ( ( aNucleon = theNucleus->GetNextNucleon() ) ) {
      Ptarget += aNucleon->Get4Momentum();
      if ( aNucleon->AreYouHit() ) {  // Involved nucleons
        SumMasses += std::sqrt( sqr( aNucleon->GetDefinition()->GetPDGMass() ) 
                                +  aNucleon->Get4Momentum().perp2() );                     
        SumMasses += 20.0*MeV;  // Separation energy for a nucleon
        TargetResidualExcitationEnergy += ExcitationEnergyPerWoundedNucleon;
        TargetResidualMassNumber--;
        TargetResidualCharge -= G4int( aNucleon->GetDefinition()->GetPDGCharge() );
      } else {   // Spectator nucleons
        PtargetResidual += aNucleon->Get4Momentum();
      }
    }

    #ifdef debugPutOnMassShell
    G4cout << "Target residual: Charge, MassNumber " << TargetResidualCharge << " "
           << TargetResidualMassNumber << G4endl << "Target Initial Momentum " << Ptarget 
           << G4endl << "Target Residual Momentum   " << PtargetResidual << G4endl;
    #endif

    Psum += Ptarget;   
    PtargetResidual.setPz( 0.0 ); PtargetResidual.setE( 0.0 );
    G4double TargetResidualMass( 0.0 );
    if ( TargetResidualMassNumber == 0 ) {
      TargetResidualMass = 0.0;
      TargetResidualExcitationEnergy = 0.0;
    } else {
      TargetResidualMass = G4ParticleTable::GetParticleTable()->GetIonTable()->
                             GetIonMass( TargetResidualCharge, TargetResidualMassNumber );
      if ( TargetResidualMassNumber == 1 ) {
        TargetResidualExcitationEnergy = 0.0;
      }
    }
    SumMasses += std::sqrt( sqr( TargetResidualMass ) + PtargetResidual.perp2() );

    G4double SqrtS = Psum.mag();
    G4double     S = Psum.mag2();

    #ifdef debugPutOnMassShell
    G4cout << "Psum " << Psum/GeV << " GeV" << G4endl << "SqrtS " << SqrtS/GeV << " GeV" << G4endl
           << "SumMasses And TargetResidualMass " << SumMasses/GeV << " " 
           << TargetResidualMass/GeV << " GeV" << G4endl;
    #endif

    if ( SqrtS < SumMasses ) {
      return false;  // It is impossible to simulate after putting nuclear nucleons on mass-shell
    }

    SumMasses -= std::sqrt( sqr( TargetResidualMass ) + PtargetResidual.perp2() );
    SumMasses += std::sqrt( sqr( TargetResidualMass + TargetResidualExcitationEnergy )
                            + PtargetResidual.perp2() ); 

    if ( SqrtS < SumMasses ) {
      SumMasses -= std::sqrt( sqr( TargetResidualMass + TargetResidualExcitationEnergy )
                              + PtargetResidual.perp2() );
      SumMasses += std::sqrt( sqr( TargetResidualMass ) + PtargetResidual.perp2() );
      TargetResidualExcitationEnergy = 0.0;
    }

    TargetResidualMass +=TargetResidualExcitationEnergy;

    #ifdef debugPutOnMassShell
    G4cout << "TargetResidualMass SumMasses TargetResidualExcitationEnergy " 
           << TargetResidualMass/GeV << " " << SumMasses/GeV << " " 
           << TargetResidualExcitationEnergy << " MeV" << G4endl;
    #endif

    // Sampling of nucleons what can transfer to delta-isobars

    //G4cout << "Sampling of nucleons what can transfer to delta-isobars ----" << G4endl
    //       << SqrtS/GeV << " " << SumMasses/GeV << G4endl
    //       << "NumberOfInvolvedNucleonsOfTarget " << NumberOfInvolvedNucleonsOfTarget << G4endl;

    G4int MaxNumberOfDeltas = G4int( (SqrtS - SumMasses)/(400.0*MeV) );
    G4int NumberOfDeltas( 0 );

    if ( theNucleus->GetMassNumber() != 1 ) {
      //G4double ProbDeltaIsobar( 0.05 );  // Uzhi 6.07.2012
      //G4double ProbDeltaIsobar( 0.25 );  // Uzhi 13.06.2013 
      G4double ProbDeltaIsobar( 0.10 );  // A.R. 07.08.2013
      for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; i++ ) {
        //G4cout << "i MaxNumberOfDeltas ProbDeltaIsobar " << i << " " << MaxNumberOfDeltas
        //       << " " << ProbDeltaIsobar << G4endl;
        if ( G4UniformRand() < ProbDeltaIsobar  &&  NumberOfDeltas < MaxNumberOfDeltas ) {
          NumberOfDeltas++;
          G4VSplitableHadron* targetSplitable = 
              TheInvolvedNucleonsOfTarget[i]->GetSplitableHadron();
          G4double MassNuc = std::sqrt( sqr( targetSplitable->GetDefinition()->GetPDGMass() )
                                       + targetSplitable->Get4Momentum().perp2() );
          G4int PDGcode = targetSplitable->GetDefinition()->GetPDGEncoding();
          G4ParticleDefinition* Old_def = targetSplitable->GetDefinition();
          G4int newPDGcode = PDGcode/10; newPDGcode = newPDGcode*10 + 4; // Delta
          G4ParticleDefinition* ptr = 
              G4ParticleTable::GetParticleTable()->FindParticle( newPDGcode );
          targetSplitable->SetDefinition( ptr );
          G4double MassDel = std::sqrt( sqr( targetSplitable->GetDefinition()->GetPDGMass() )
                                       + targetSplitable->Get4Momentum().perp2() );
          //G4cout << i << " " << SqrtS/GeV << " " << SumMasses/GeV << " " << MassDel/GeV
          //       << " " << MassNuc << G4endl;
          if ( SqrtS < SumMasses + MassDel - MassNuc ) {  // Uzhi 12.06.2012
            // Change cannot be accepted!
            targetSplitable->SetDefinition( Old_def );
            ProbDeltaIsobar = 0.0;
          } else {  // Change is accepted
            SumMasses += (MassDel - MassNuc);
          }
        } 
      }
    }
    //G4cout << "MaxNumberOfDeltas NumberOfDeltas " << MaxNumberOfDeltas << " " 
    //       << NumberOfDeltas << G4endl;

    G4LorentzRotation toCms( -1*Psum.boostVector() );
    G4LorentzVector Ptmp = toCms*Pprojectile;
    if ( Ptmp.pz() <= 0.0 ) {  // "String" moving backwards in  CMS, abort collision!
      //G4cout << " abort ColliDeleteVSplitableHadronsion! " << G4endl;
      return false; 
    }

    G4LorentzRotation toLab( toCms.inverse() );
    Ptmp = toCms*Ptarget;                      
    G4double YtargetNucleus = Ptmp.rapidity();

    // Ascribing of the involved nucleons Pt and Xminus
    G4double Dcor = theParameters->GetDofNuclearDestruction() / theNucleus->GetMassNumber();
    G4double AveragePt2  = theParameters->GetPt2ofNuclearDestruction();
    G4double maxPtSquare = theParameters->GetMaxPt2ofNuclearDestruction();

    #ifdef debugPutOnMassShell
    G4cout << "Y targetNucleus " << YtargetNucleus << G4endl
           << "Dcor " << theParameters->GetDofNuclearDestruction() << " DcorA " << Dcor
           << " AveragePt2 " << AveragePt2 << G4endl;
    #endif

    G4double M2target( 0.0 );
    G4double WminusTarget( 0.0 );
    G4double WplusProjectile( 0.0 );

    G4int NumberOfTries( 0 );
    G4double ScaleFactor( 1.0 );
    G4bool OuterSuccess( true );

    do {  //while ( ! OuterSuccess )
 
      OuterSuccess = true;

      do {  // while ( SqrtS < Mprojectile + std::sqrt( M2target ) )

        NumberOfTries++;

        if ( NumberOfTries == 100*(NumberOfTries/100) ) {
          // At large number of tries it would be better to reduce the values
          ScaleFactor /= 2.0;
          Dcor       *= ScaleFactor;
          AveragePt2 *= ScaleFactor;
        }

        G4ThreeVector PtSum( 0.0, 0.0, 0.0 );
        G4double XminusSum( 0.0 );
        G4double Xminus( 0.0 );
        G4bool InerSuccess = true;

        do {  // while ( ! InerSuccess )

          InerSuccess = true;
          PtSum = G4ThreeVector( 0.0, 0.0, 0.0 );
          XminusSum = 0.0;

          for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; i++ ) {
            aNucleon = TheInvolvedNucleonsOfTarget[i];
            G4ThreeVector tmpPt = GaussianPt( AveragePt2, maxPtSquare );
            PtSum += tmpPt;
            G4ThreeVector tmpX = GaussianPt( Dcor*Dcor, 1.0 );
            Xminus = tmpX.x();
            XminusSum += Xminus;
            G4LorentzVector tmp( tmpPt.x(), tmpPt.y(), Xminus, aNucleon->Get4Momentum().e() );
            aNucleon->SetMomentum( tmp );
          }

          G4double DeltaX = ( PtSum.x() - PtargetResidual.x() )/NumberOfInvolvedNucleonsOfTarget;
          G4double DeltaY = ( PtSum.y() - PtargetResidual.y() )/NumberOfInvolvedNucleonsOfTarget;
          G4double DeltaXminus( 0.0 );
          if ( TargetResidualMassNumber == 0 ) {
            DeltaXminus = (XminusSum - 1.0) / NumberOfInvolvedNucleonsOfTarget;
          } else {
            DeltaXminus = -1.0 / theNucleus->GetMassNumber();
          }

          XminusSum = 1.0;
          M2target = 0.0;

          for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; i++ ) {
            aNucleon = TheInvolvedNucleonsOfTarget[i];
            Xminus = aNucleon->Get4Momentum().pz() - DeltaXminus;
            XminusSum -= Xminus;               
            if ( TargetResidualMassNumber == 0 ) {
              if ( Xminus <= 0.0  ||  Xminus > 1.0 ) {
                InerSuccess = false; 
                break;
              }
            } else {
              if ( Xminus <= 0.0  ||  Xminus > 1.0  ||  XminusSum <= 0.0  ||  XminusSum > 1.0 ) {
                InerSuccess = false; 
                break;
              }
            }                                          
            G4double Px = aNucleon->Get4Momentum().px() - DeltaX;
            G4double Py = aNucleon->Get4Momentum().py() - DeltaY;
            M2target += ( sqr( aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass() )
                          + sqr( Px ) + sqr( Py ) ) / Xminus;
            G4LorentzVector tmp( Px, Py, Xminus, aNucleon->Get4Momentum().e() ); // 6.12.2010
            aNucleon->SetMomentum( tmp );
          }

          if ( InerSuccess  &&  TargetResidualMassNumber != 0 ) {
            M2target += ( sqr( TargetResidualMass ) + PtargetResidual.perp2() ) / XminusSum;
          }

          #ifdef debugPutOnMassShell
          G4cout << "InerSuccess " << InerSuccess << G4endl << "SqrtS, Mp+Mt, Mt " << SqrtS/GeV
                 << " " << ( Mprojectile + std::sqrt( M2target ) )/GeV << " " 
                 << std::sqrt( M2target )/GeV << G4endl 
                 << "To continue - enter 1, to stop - ^C" << G4endl;
          G4int Uzhi; G4cin >> Uzhi;
          #endif

        } while ( ! InerSuccess );

      } while ( SqrtS < Mprojectile + std::sqrt( M2target ) );

      G4double DecayMomentum2 = sqr( S ) + sqr( M2projectile ) + sqr( M2target )
                                - 2.0*S*M2projectile - 2.0*S*M2target - 2.0*M2projectile*M2target;
      WminusTarget = ( S - M2projectile + M2target + std::sqrt( DecayMomentum2 ) ) / 2.0 / SqrtS;
      WplusProjectile = SqrtS - M2target / WminusTarget;
      G4double Pzprojectile = WplusProjectile/2.0 - M2projectile/2.0/WplusProjectile;
      G4double Eprojectile  = WplusProjectile/2.0 + M2projectile/2.0/WplusProjectile;
      G4double Yprojectile = 0.5 * std::log( (Eprojectile + Pzprojectile)/
                                             (Eprojectile - Pzprojectile) );

      #ifdef debugPutOnMassShell
      G4cout << "DecayMomentum2 " << DecayMomentum2 << G4endl
             << "WminusTarget WplusProjectile " << WminusTarget << " " << WplusProjectile
             << G4endl << "Yprojectile " << Yprojectile << G4endl;
      #endif

      // Now all is O.K
      for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; i++ ) {
        aNucleon = TheInvolvedNucleonsOfTarget[i];
        G4LorentzVector tmp = aNucleon->Get4Momentum();
        G4double Mt2 = sqr( tmp.x() ) + sqr( tmp.y() ) +
                       sqr( aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass() );
        G4double Xminus = tmp.z();
        G4double Pz = -WminusTarget*Xminus/2.0 + Mt2/(2.0*WminusTarget*Xminus);
        G4double E  =  WminusTarget*Xminus/2.0 + Mt2/(2.0*WminusTarget*Xminus);
        G4double YtargetNucleon = 0.5 * std::log( (E + Pz)/(E - Pz) ); 

        #ifdef debugPutOnMassShell
        G4cout << "i YtN Ypr YtN-YtA " << i << " " << YtargetNucleon << " " << YtargetNucleus
               << " " << YtargetNucleon - YtargetNucleus << G4endl;
        #endif

        if ( std::abs( YtargetNucleon - YtargetNucleus ) > 2  ||  Yprojectile < YtargetNucleon ) {
          OuterSuccess = false; 
          break;
        } 
      }  

    } while ( ! OuterSuccess );

    G4double Pzprojectile = WplusProjectile/2.0 - M2projectile/2.0/WplusProjectile;
    G4double Eprojectile  = WplusProjectile/2.0 + M2projectile/2.0/WplusProjectile;
    Pprojectile.setPz( Pzprojectile ); Pprojectile.setE( Eprojectile );

    #ifdef debugPutOnMassShell
    G4cout << "Proj after in CMS " << Pprojectile << G4endl;
    #endif

    // The work with the projectile is finished at the moment.
    Pprojectile.transform( toLab );  
    theProjectile.SetMomentum( Pprojectile.vect() );
    theProjectile.SetTotalEnergy( Pprojectile.e() );

    theParticipants.StartLoop();
    theParticipants.Next();
    G4VSplitableHadron* primary = theParticipants.GetInteraction().GetProjectile();
    primary->Set4Momentum( Pprojectile ); 

    #ifdef debugPutOnMassShell
    G4cout << "Final proj. mom in Lab. " << primary->Get4Momentum() << G4endl;
    #endif

    G4ThreeVector TargetResidual3Momentum( 0.0, 0.0, 1.0 );

    for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; i++ ) {
      aNucleon = TheInvolvedNucleonsOfTarget[i];
      G4LorentzVector tmp = aNucleon->Get4Momentum();
      TargetResidual3Momentum -= tmp.vect();
      G4double Mt2 = sqr( tmp.x() ) + sqr( tmp.y() ) +
                     sqr( aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass() );
      G4double Xminus = tmp.z();
      G4double Pz = -WminusTarget*Xminus/2.0 + Mt2/(2.0*WminusTarget*Xminus);
      G4double E  =  WminusTarget*Xminus/2.0 + Mt2/(2.0*WminusTarget*Xminus);
      tmp.setPz( Pz ); 
      tmp.setE( E );
      tmp.transform( toLab );
      aNucleon->SetMomentum( tmp );
      G4VSplitableHadron* targetSplitable = aNucleon->GetSplitableHadron();
      targetSplitable->Set4Momentum( tmp );
    }

    G4double Mt2Residual = sqr( TargetResidualMass ) + sqr( TargetResidual3Momentum.x() )
                         + sqr( TargetResidual3Momentum.y() );

    #ifdef debugPutOnMassShell
    G4cout << "WminusTarget TargetResidual3Momentum.z() " << WminusTarget << " "
           << TargetResidual3Momentum.z() << G4endl;
    #endif

    G4double PzResidual = 0.0;
    G4double EResidual  = 0.0;
    if ( TargetResidualMassNumber != 0 ) {
      PzResidual = -WminusTarget*TargetResidual3Momentum.z()/2.0 + 
                    Mt2Residual/(2.0*WminusTarget*TargetResidual3Momentum.z());
      EResidual  =  WminusTarget*TargetResidual3Momentum.z()/2.0 + 
                    Mt2Residual/(2.0*WminusTarget*TargetResidual3Momentum.z());
    }

    TargetResidual4Momentum.setPx( TargetResidual3Momentum.x() );
    TargetResidual4Momentum.setPy( TargetResidual3Momentum.y() );
    TargetResidual4Momentum.setPz( PzResidual ); 
    TargetResidual4Momentum.setE( EResidual );

    #ifdef debugPutOnMassShell
    G4cout << "Target Residual4Momentum in CMS" << TargetResidual4Momentum << G4endl;
    #endif

    TargetResidual4Momentum.transform( toLab );

    #ifdef debugPutOnMassShell
    G4cout << "Target Residual4Momentum in Lab " << TargetResidual4Momentum << G4endl
           << "To continue enter - 1, to break - ^C" << G4endl;
    G4int Uzhi; G4cin >> Uzhi;
    #endif

    return true;

  }  // end if ( ! GetProjectileNucleus() )

  // The projectile is a nucleus or an anti-nucleus

  #ifdef debugPutOnMassShell
  G4cout << "PutOnMassShell for Nucleus_Nucleus " << G4endl;
  #endif

  G4LorentzVector Pprojectile( theProjectile.GetMomentum(), theProjectile.GetTotalEnergy() );

  if ( Pprojectile.z() < 0.0 ) {
    return false;
  }

  G4double ExcitationEnergyPerWoundedNucleon = theParameters->GetExcitationEnergyPerWoundedNucleon();

  #ifdef debugPutOnMassShell
  G4cout << "ExcitationEnergyPerWoundedNucleon " << ExcitationEnergyPerWoundedNucleon << G4endl;
  #endif

  G4LorentzVector Psum = Pprojectile;
  G4double SumMasses( 0.0 );

  // Projectile nucleus
  G4V3DNucleus* thePrNucleus = GetProjectileNucleus();
  //G4int PrResidualMassNumber = thePrNucleus->GetMassNumber();
  //G4int PrResidualCharge     = thePrNucleus->GetCharge();
  //ProjectileResidualExcitationEnergy = 0.0;
  G4LorentzVector Pproj( 0.0, 0.0, 0.0, 0.0 );
  G4LorentzVector PprojResidual( 0.0, 0.0, 0.0, 0.0 );
  thePrNucleus->StartLoop();
  G4Nucleon* aNucleon( 0 );
  while ( ( aNucleon = thePrNucleus->GetNextNucleon() ) ) {
    Pproj += aNucleon->Get4Momentum();
    if ( aNucleon->AreYouHit() ) {  // Involved nucleons
      SumMasses += std::sqrt( sqr( aNucleon->GetDefinition()->GetPDGMass() ) + 
                              aNucleon->Get4Momentum().perp2() );                     
      SumMasses += 20.0*MeV;  // Separation energy for a nucleon
      ProjectileResidualExcitationEnergy += ExcitationEnergyPerWoundedNucleon;
      ProjectileResidualMassNumber--;
      ProjectileResidualCharge -= std::abs( G4int( aNucleon->GetDefinition()->GetPDGCharge() ) );
    } else {  // Spectator nucleons
      PprojResidual += aNucleon->Get4Momentum();
    }
  }

  #ifdef debugPutOnMassShell
  G4cout << "Projectile residual: Charge, MassNumber " << ProjectileResidualCharge << " "
         << ProjectileResidualMassNumber << G4endl << "Initial Momentum   " << Pproj << G4endl
         << "Residual Momentum " << PprojResidual << G4endl;
  #endif

  // Target nucleus
  G4V3DNucleus* theNucleus = GetTargetNucleus();
  //TargetResidualMassNumber = theNucleus->GetMassNumber();
  //TargetResidualCharge     = theNucleus->GetCharge();
  //TargetResidualExcitationEnergy = 0.0;
  G4LorentzVector Ptarget( 0.0, 0.0, 0.0, 0.0 );
  G4LorentzVector PtargetResidual( 0.0, 0.0, 0.0, 0.0 );
  theNucleus->StartLoop();
  aNucleon = 0;
  while ( ( aNucleon = theNucleus->GetNextNucleon() ) ) {
    Ptarget += aNucleon->Get4Momentum();
    if ( aNucleon->AreYouHit() ) {  // Involved nucleons
      SumMasses += std::sqrt( sqr( aNucleon->GetDefinition()->GetPDGMass() ) + 
                              aNucleon->Get4Momentum().perp2() );                     
      SumMasses += 20.0*MeV;  // Separation energy for a nucleon
      TargetResidualExcitationEnergy += ExcitationEnergyPerWoundedNucleon;
      TargetResidualMassNumber--;
      TargetResidualCharge -= G4int( aNucleon->GetDefinition()->GetPDGCharge() );
    } else {  // Spectator nucleons
      PtargetResidual += aNucleon->Get4Momentum();
    }
  }

  #ifdef debugPutOnMassShell
  G4cout << "Target residual: Charge, MassNumber " << TargetResidualCharge << " "
         << TargetResidualMassNumber << G4endl << "Initial Momentum   " << Ptarget << G4endl
         << "Residual Momentum " << PtargetResidual << G4endl;
  #endif

  Psum += Ptarget;   

  PprojResidual.setPz( 0.0 );   PprojResidual.setE( 0.0 );
  PtargetResidual.setPz( 0.0 ); PtargetResidual.setE( 0.0 );

  G4double PrResidualMass( 0.0 );
  if ( ProjectileResidualMassNumber == 0 ) {
    PrResidualMass = 0.0;
    ProjectileResidualExcitationEnergy = 0.0;
  } else {
    PrResidualMass = G4ParticleTable::GetParticleTable()->GetIonTable()
                         ->GetIonMass( ProjectileResidualCharge, ProjectileResidualMassNumber );
    if ( ProjectileResidualMassNumber == 1 ) {
      ProjectileResidualExcitationEnergy = 0.0;
    }
  }

  SumMasses += std::sqrt( sqr( PrResidualMass ) + PprojResidual.vect().perp2() );

  G4double TargetResidualMass( 0.0 );
  if ( TargetResidualMassNumber == 0 ) {
    TargetResidualMass = 0.0;
    TargetResidualExcitationEnergy = 0.0;
  } else {
    TargetResidualMass = G4ParticleTable::GetParticleTable()->GetIonTable()
                           ->GetIonMass( TargetResidualCharge, TargetResidualMassNumber );
    if ( TargetResidualMassNumber == 1 ) {
      TargetResidualExcitationEnergy = 0.0;
    }
  }

  SumMasses += std::sqrt( sqr( TargetResidualMass ) + PtargetResidual.perp2() );

  G4double SqrtS = Psum.mag();
  G4double     S = Psum.mag2();

  #ifdef debugPutOnMassShell
  G4cout << "Psum " << Psum/GeV << " GeV" << G4endl << "SqrtS " << SqrtS/GeV << " GeV" << G4endl
         << "SumMasses, PrResidualMass and TargetResidualMass " << SumMasses/GeV << " " 
         << PrResidualMass/GeV << " " << TargetResidualMass/GeV << " GeV" << G4endl;
  #endif

  if ( SqrtS < SumMasses ) {
    return false; // It is impossible to simulate after putting nuclear nucleons on mass-shell. 
  }

  SumMasses -= std::sqrt( sqr( PrResidualMass ) + PprojResidual.perp2() );
  SumMasses += std::sqrt( sqr( PrResidualMass + ProjectileResidualExcitationEnergy ) + 
                          PprojResidual.perp2() ); 
  SumMasses -= std::sqrt( sqr( TargetResidualMass ) + PtargetResidual.perp2() );
  SumMasses += std::sqrt( sqr( TargetResidualMass + TargetResidualExcitationEnergy ) +
                          PtargetResidual.perp2() ); 

  if ( SqrtS < SumMasses ) {
    SumMasses -= std::sqrt( sqr( PrResidualMass + ProjectileResidualExcitationEnergy ) +
                            PprojResidual.perp2() );
    SumMasses += std::sqrt( sqr( PrResidualMass ) + PprojResidual.perp2() );
    ProjectileResidualExcitationEnergy = 0.0;
    SumMasses -= std::sqrt( sqr( TargetResidualMass + TargetResidualExcitationEnergy ) +
                            PtargetResidual.perp2() );
    SumMasses += std::sqrt( sqr( TargetResidualMass ) + PtargetResidual.perp2() );
    TargetResidualExcitationEnergy = 0.0;
  }

  PrResidualMass     += ProjectileResidualExcitationEnergy;
  TargetResidualMass += TargetResidualExcitationEnergy;

  #ifdef debugPutOnMassShell
  G4cout << "PrResidualMass ProjResidualExcitationEnergy " << PrResidualMass/GeV << " "
         << ProjectileResidualExcitationEnergy << " MeV" << G4endl
         << "TargetResidualMass TargetResidualExcitationEnergy " << TargetResidualMass/GeV << " "
         << TargetResidualExcitationEnergy << " MeV" << G4endl
         << "Sum masses " << SumMasses/GeV << G4endl;
  #endif

  // Sampling of projectile nucleons what can transfer to delta-isobars

  G4int MaxNumberOfDeltas = G4int( (SqrtS - SumMasses)/(400.0*MeV) );
  G4int NumberOfDeltas( 0 );

  //G4double ProbDeltaIsobar( 0.05 );  // Uzhi 6.07.2012
  //G4double ProbDeltaIsobar( 0.25 );  // Uzhi 13.06.2013 
  G4double ProbDeltaIsobar( 0.10 );  // A.R. 07.08.2013
  if ( thePrNucleus->GetMassNumber() != 1 ) {
    for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfProjectile; i++ ) {
      if ( G4UniformRand() < ProbDeltaIsobar  &&  NumberOfDeltas < MaxNumberOfDeltas ) {
        NumberOfDeltas++;
        G4VSplitableHadron* projectileSplitable = 
            TheInvolvedNucleonsOfProjectile[i]->GetSplitableHadron();
        G4double MassNuc = std::sqrt( sqr( projectileSplitable->GetDefinition()->GetPDGMass() ) +
                                      projectileSplitable->Get4Momentum().perp2() );
        G4int PDGcode = std::abs( projectileSplitable->GetDefinition()->GetPDGEncoding() );
        G4ParticleDefinition* Old_def = projectileSplitable->GetDefinition();
        G4int newPDGcode = PDGcode/10; newPDGcode = newPDGcode*10 + 4; // Delta
        if ( projectileSplitable->GetDefinition()->GetPDGEncoding() < 0 ) newPDGcode *= -1;
        G4ParticleDefinition* ptr = G4ParticleTable::GetParticleTable()->FindParticle(newPDGcode);
        projectileSplitable->SetDefinition( ptr );
        G4double MassDel = std::sqrt( sqr( projectileSplitable->GetDefinition()->GetPDGMass() ) +
                                      projectileSplitable->Get4Momentum().perp2() );
        if ( SqrtS < SumMasses + MassDel - MassNuc ) { // Change cannot be accepted!
          projectileSplitable->SetDefinition( Old_def );
          ProbDeltaIsobar = 0.0;
        } else { // Change is accepted.
          SumMasses += (MassDel - MassNuc);
        }
      } 
    }
  }

  // Sampling of target nucleons what can transfer to delta-isobars
  if ( theNucleus->GetMassNumber() != 1 ) {
    //G4double ProbDeltaIsobar( 0.05 );  // Uzhi 6.07.2012
    for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; i++ ) {
      if ( G4UniformRand() < ProbDeltaIsobar  &&  NumberOfDeltas < MaxNumberOfDeltas ) {
        NumberOfDeltas++;
        G4VSplitableHadron* targetSplitable = 
            TheInvolvedNucleonsOfTarget[i]->GetSplitableHadron();
        G4double MassNuc = std::sqrt( sqr( targetSplitable->GetDefinition()->GetPDGMass() ) +
                                           targetSplitable->Get4Momentum().perp2() );
        G4int PDGcode = targetSplitable->GetDefinition()->GetPDGEncoding();
        G4ParticleDefinition* Old_def = targetSplitable->GetDefinition();
        G4int newPDGcode = PDGcode/10; newPDGcode=newPDGcode*10 + 4;  // Delta
        G4ParticleDefinition* ptr = G4ParticleTable::GetParticleTable()->FindParticle(newPDGcode);
        targetSplitable->SetDefinition( ptr );
        G4double MassDel = std::sqrt( sqr( targetSplitable->GetDefinition()->GetPDGMass() ) +
                                      targetSplitable->Get4Momentum().perp2() );
        if ( SqrtS < SumMasses + MassDel - MassNuc ) {  // Uzhi 12.06.2012
          targetSplitable->SetDefinition( Old_def );  // Change cannot be accepted!
          ProbDeltaIsobar = 0.0;
        } else {  // Change is accepted.
          SumMasses += (MassDel - MassNuc);
        }
      } 
    }
  }

  G4LorentzRotation toCms( -1*Psum.boostVector() );
  G4LorentzVector Ptmp = toCms*Pprojectile;
  if ( Ptmp.pz() <= 0.0 ) {  // "String" moving backwards in CMS, abort collision !
    //G4cout << " abort ColliDeleteVSplitableHadronsion! " << G4endl;
    return false; 
  }

  G4LorentzRotation toLab( toCms.inverse() );

  Ptmp = toCms*Pproj;                      
  G4double YprojectileNucleus = Ptmp.rapidity();
  Ptmp = toCms*Ptarget;                      
  G4double YtargetNucleus = Ptmp.rapidity();

  // Ascribing of the involved nucleons Pt and Xminus
  G4double DcorP       = theParameters->GetDofNuclearDestruction()/thePrNucleus->GetMassNumber();
  G4double DcorT       = theParameters->GetDofNuclearDestruction()/theNucleus->GetMassNumber();
  G4double AveragePt2  = theParameters->GetPt2ofNuclearDestruction();
  G4double maxPtSquare = theParameters->GetMaxPt2ofNuclearDestruction();

  #ifdef debugPutOnMassShell
  G4cout << "Y projectileNucleus " << YprojectileNucleus << G4endl << "Y targetNucleus     " 
         << YtargetNucleus << G4endl << "Dcor " << theParameters->GetDofNuclearDestruction()
         << " DcorP DcorT " << DcorP << " " << DcorT << " AveragePt2 " << AveragePt2 << G4endl;
  #endif

  G4double M2proj( 0.0 );
  G4double WplusProjectile( 0.0 );
  G4double M2target( 0.0 );
  G4double WminusTarget( 0.0 );
  G4int NumberOfTries( 0 );
  G4double ScaleFactor( 1.0 );
  G4bool OuterSuccess( true );
  
  do {  //  while ( ! OuterSuccess )
  
    OuterSuccess = true;

    do {  // while ( SqrtS < std::sqrt( M2proj ) + std::sqrt( M2target ) )

      NumberOfTries++;

      if ( NumberOfTries == 100*(NumberOfTries/100) ) {
        // At large number of tries it would be better to reduce the values
        ScaleFactor /= 2.0;
        DcorP       *= ScaleFactor;
        DcorT       *= ScaleFactor;
        AveragePt2  *= ScaleFactor;
      }

      // Sampling of kinematical properties of projectile nucleons
      G4ThreeVector PtSum( 0.0, 0.0, 0.0 );
      G4double XplusSum( 0.0 );
      G4double Xplus( 0.0 );
      G4bool InerSuccess = true;

      do {  // while ( ! InerSuccess )
           
        InerSuccess = true;

        PtSum = G4ThreeVector( 0.0, 0.0, 0.0 );
        XplusSum = 0.0;

        for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfProjectile; i++ ) {
          aNucleon = TheInvolvedNucleonsOfProjectile[i];
          G4ThreeVector tmpPt = GaussianPt( AveragePt2, maxPtSquare );
          PtSum += tmpPt;
          G4ThreeVector tmpX = GaussianPt( DcorP*DcorP, 1.0 );
          Xplus = tmpX.x();
          XplusSum += Xplus;
          G4LorentzVector tmp( tmpPt.x(), tmpPt.y(), Xplus, aNucleon->Get4Momentum().e() );
          aNucleon->SetMomentum( tmp );
        }

        G4double DeltaX = (PtSum.x() - PprojResidual.x()) / NumberOfInvolvedNucleonsOfProjectile;
        G4double DeltaY = (PtSum.y() - PprojResidual.y()) / NumberOfInvolvedNucleonsOfProjectile;
        G4double DeltaXplus( 0.0 );
        if ( ProjectileResidualMassNumber == 0 ) {
          DeltaXplus = (XplusSum - 1.0) / NumberOfInvolvedNucleonsOfProjectile;
        } else {
          DeltaXplus = -1.0 / thePrNucleus->GetMassNumber();
        }
        XplusSum = 1.0;
        M2proj= 0.0;

        for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfProjectile; i++ ) {
          aNucleon = TheInvolvedNucleonsOfProjectile[i];
          Xplus = aNucleon->Get4Momentum().pz() - DeltaXplus;
          XplusSum -= Xplus;               
          if ( ProjectileResidualMassNumber == 0 ) {
            if ( Xplus <= 0.0  ||  Xplus > 1.0 ) {
              InerSuccess = false; 
              break;
            }
          } else {
            if ( Xplus <= 0.0  ||  Xplus > 1.0  ||  XplusSum <= 0.0  ||  XplusSum > 1.0 ) {
              InerSuccess = false; 
              break;
            }
          }                                          
          G4double Px = aNucleon->Get4Momentum().px() - DeltaX;
          G4double Py = aNucleon->Get4Momentum().py() - DeltaY;
          M2proj +=( sqr( aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass() ) +
                     sqr( Px ) + sqr( Py ) ) / Xplus;
          G4LorentzVector tmp( Px, Py, Xplus, aNucleon->Get4Momentum().e() ); // 6.12.2010
               aNucleon->SetMomentum(tmp);
        }

        if ( InerSuccess  &&  ProjectileResidualMassNumber != 0 ) {
          M2proj += ( sqr( PrResidualMass ) + PprojResidual.perp2() ) / XplusSum;
        }

      } while ( ! InerSuccess );

      // Sampling of kinematical properties of target nucleons

      G4double XminusSum( 0.0 );
      G4double Xminus( 0.0);

      do {  // while ( ! InerSuccess )

        InerSuccess = true;

        PtSum = G4ThreeVector( 0.0, 0.0, 0.0 );
        XminusSum = 0.0;

        for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; i++ ) {
          aNucleon = TheInvolvedNucleonsOfTarget[i];
          G4ThreeVector tmpPt = GaussianPt( AveragePt2, maxPtSquare );
          PtSum += tmpPt;
          G4ThreeVector tmpX = GaussianPt( DcorT*DcorT, 1.0 );
          Xminus = tmpX.x();
          XminusSum += Xminus;
          G4LorentzVector tmp( tmpPt.x(), tmpPt.y(), Xminus, aNucleon->Get4Momentum().e() );
          aNucleon->SetMomentum( tmp );
        }

        G4double DeltaX = (PtSum.x() - PtargetResidual.x()) / NumberOfInvolvedNucleonsOfTarget;
        G4double DeltaY = (PtSum.y() - PtargetResidual.y()) / NumberOfInvolvedNucleonsOfTarget;
        G4double DeltaXminus( 0.0 );
        if ( TargetResidualMassNumber == 0 ) {
          DeltaXminus = (XminusSum - 1.0) / NumberOfInvolvedNucleonsOfTarget;
        } else {
          DeltaXminus = -1.0 / theNucleus->GetMassNumber();
        }
        
        XminusSum = 1.0;
        M2target = 0.0;

        for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; i++ ) {
          aNucleon = TheInvolvedNucleonsOfTarget[i];
          Xminus = aNucleon->Get4Momentum().pz() - DeltaXminus;
          XminusSum -= Xminus;               
          if ( TargetResidualMassNumber == 0 ) {
            if ( Xminus <= 0.0  ||  Xminus > 1.0 ) {
              InerSuccess = false; 
              break;
            }
          } else {
            if ( Xminus <= 0.0  ||  Xminus > 1.0  ||  XminusSum <= 0.0  ||  XminusSum > 1.0 ) {
              InerSuccess = false; 
              break;
            }
          }                                          
          G4double Px = aNucleon->Get4Momentum().px() - DeltaX;
          G4double Py = aNucleon->Get4Momentum().py() - DeltaY;
          M2target += ( sqr( aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass() ) +
                        sqr( Px ) + sqr( Py ) ) / Xminus;
          G4LorentzVector tmp( Px, Py, Xminus, aNucleon->Get4Momentum().e() ); // 6.12.2010
          aNucleon->SetMomentum( tmp );
        }

        if ( InerSuccess  &&  TargetResidualMassNumber != 0 ) {
          M2target += ( sqr( TargetResidualMass ) + PtargetResidual.perp2() ) / XminusSum;
        }
          
      } while ( ! InerSuccess );

      #ifdef debugPutOnMassShell
      G4cout << "SqrtS, Mp+Mt, Mp, Mt " << SqrtS/GeV << " " 
             << ( std::sqrt( M2proj ) + std::sqrt( M2target) )/GeV << " "
             << std::sqrt( M2proj )/GeV << " " << std::sqrt( M2target )/GeV << G4endl
             << "To continue - enter 1, to stop - ^C" << G4endl;
      G4int Uzhi; G4cin >> Uzhi;
      #endif

    } while ( SqrtS < std::sqrt( M2proj ) + std::sqrt( M2target ) );

    G4double DecayMomentum2 = sqr( S ) + sqr( M2proj ) + sqr( M2target )
                              - 2.0*S*M2proj - 2.0*S*M2target - 2.0*M2proj*M2target;
    WminusTarget = ( S - M2proj + M2target + std::sqrt( DecayMomentum2 ) ) / 2.0 / SqrtS;
    WplusProjectile = SqrtS - M2target/WminusTarget;
    G4double Pzprojectile = WplusProjectile/2.0 - M2proj/2.0/WplusProjectile;
    G4double Eprojectile  = WplusProjectile/2.0 + M2proj/2.0/WplusProjectile;
    G4double Yprojectile  = 0.5 * std::log( (Eprojectile + Pzprojectile)/
                                            (Eprojectile - Pzprojectile) );
    G4double Pztarget = -WminusTarget/2.0 + M2target/2.0/WminusTarget;
    G4double Etarget  =  WminusTarget/2.0 + M2target/2.0/WminusTarget;
    G4double Ytarget  = 0.5 * std::log( (Etarget + Pztarget)/(Etarget - Pztarget) );

    #ifdef debugPutOnMassShell
    G4cout << "DecayMomentum2 " << DecayMomentum2 << G4endl << "WminusTarget WplusProjectile "
           << WminusTarget << " " << WplusProjectile << G4endl
           << "Yprojectile " << Yprojectile << G4endl;
    #endif

    // Now all is O.K.

    for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfProjectile; i++ ) {
      aNucleon = TheInvolvedNucleonsOfProjectile[i];
      G4LorentzVector tmp = aNucleon->Get4Momentum();
      G4double Mt2 = sqr( tmp.x() ) + sqr( tmp.y() ) +
                     sqr( aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass() );
      G4double Xplus = tmp.z();
      G4double Pz = WplusProjectile*Xplus/2.0 - Mt2/(2.0*WplusProjectile*Xplus);
      G4double E =  WplusProjectile*Xplus/2.0 + Mt2/(2.0*WplusProjectile*Xplus);
      G4double YprojNucleon = 0.5 * std::log( (E + Pz)/(E - Pz) ); 

      #ifdef debugPutOnMassShell
      G4cout << "i YpN Ypr YpN-YtA Ypr0 " << i << " " << YprojNucleon << " " << Yprojectile
             << " " << YprojNucleon - Yprojectile << " " << YprojectileNucleus << G4endl;
      #endif

      if ( std::abs( YprojNucleon - YprojectileNucleus ) > 2  ||  Ytarget > YprojNucleon ) {
        OuterSuccess = false; 
        break;
      } 
      // Yprojectile YprojectileNucleus ???? 
    }

    for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; i++ ) {
      aNucleon = TheInvolvedNucleonsOfTarget[i];
      G4LorentzVector tmp = aNucleon->Get4Momentum();
      G4double Mt2 = sqr( tmp.x() ) + sqr( tmp.y() ) +
                     sqr( aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass() );
      G4double Xminus = tmp.z();
      G4double Pz = -WminusTarget*Xminus/2.0 + Mt2/(2.0*WminusTarget*Xminus);
      G4double E =   WminusTarget*Xminus/2.0 + Mt2/(2.0*WminusTarget*Xminus);
      G4double YtargetNucleon = 0.5 * std::log( (E + Pz)/(E - Pz) ); 

      #ifdef debugPutOnMassShell
      G4cout << "i YtN Ypr YtN-YtA " << i << " " << YtargetNucleon << " " << Yprojectile
             << " " << YtargetNucleon - YtargetNucleus << G4endl;
      #endif

      if ( std::abs( YtargetNucleon - YtargetNucleus ) > 2  ||  Yprojectile < YtargetNucleon ) {
        OuterSuccess = false; 
        break;
      } 
    }

  } while ( ! OuterSuccess );

  G4ThreeVector ProjectileResidual3Momentum( 0.0, 0.0, 1.0 );

  for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfProjectile; i++ ) {
    aNucleon = TheInvolvedNucleonsOfProjectile[i];
    G4LorentzVector tmp = aNucleon->Get4Momentum();
    ProjectileResidual3Momentum -= tmp.vect();
    G4double Mt2 = sqr( tmp.x() ) + sqr( tmp.y() ) +
                   sqr( aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass() );
    G4double Xplus = tmp.z();
    G4double Pz = WplusProjectile*Xplus/2.0 - Mt2/(2.0*WplusProjectile*Xplus);
    G4double E =  WplusProjectile*Xplus/2.0 + Mt2/(2.0*WplusProjectile*Xplus);
    tmp.setPz( Pz ); 
    tmp.setE( E );
    tmp.transform( toLab );
    aNucleon->SetMomentum( tmp );
    G4VSplitableHadron* projectileSplitable = aNucleon->GetSplitableHadron();
    projectileSplitable->Set4Momentum( tmp );
  }

  G4double Mt2PrResidual = sqr( PrResidualMass ) + sqr( ProjectileResidual3Momentum.x() ) +
                           sqr( ProjectileResidual3Momentum.y() );

  #ifdef debugPutOnMassShell
  G4cout << "WplusProjectile ProjectileResidual3Momentum.z() " << WplusProjectile << " "
         << ProjectileResidual3Momentum.z() << G4endl;
  #endif

  G4double PzPrResidual = 0.0;
  G4double EPrResidual  = 0.0;
  if ( ProjectileResidualMassNumber != 0 ) {
    PzPrResidual = WplusProjectile*ProjectileResidual3Momentum.z()/2.0 - 
                   Mt2PrResidual/( 2.0*WplusProjectile*ProjectileResidual3Momentum.z() );
    EPrResidual  = WplusProjectile*ProjectileResidual3Momentum.z()/2.0 + 
                   Mt2PrResidual/( 2.0*WplusProjectile*ProjectileResidual3Momentum.z() );
  }
  ProjectileResidual4Momentum.setPx( ProjectileResidual3Momentum.x());
  ProjectileResidual4Momentum.setPy( ProjectileResidual3Momentum.y());
  ProjectileResidual4Momentum.setPz( PzPrResidual ); 
  ProjectileResidual4Momentum.setE( EPrResidual );

  #ifdef debugPutOnMassShell
  G4cout << "Projectile Residual4Momentum in CMS" << ProjectileResidual4Momentum << G4endl;
  #endif

  ProjectileResidual4Momentum.transform( toLab );

  #ifdef debugPutOnMassShell
  G4cout << "Projectile Residual4Momentum in Lab " << ProjectileResidual4Momentum << G4endl;
  #endif

  G4ThreeVector TargetResidual3Momentum( 0.0, 0.0, 1.0 );
  for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; i++ ) {
    aNucleon = TheInvolvedNucleonsOfTarget[i];
    G4LorentzVector tmp = aNucleon->Get4Momentum();
    TargetResidual3Momentum -= tmp.vect();
    G4double Mt2 = sqr( tmp.x() ) + sqr( tmp.y() ) +
                   sqr( aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass() );
    G4double Xminus = tmp.z();
    G4double Pz = -WminusTarget*Xminus/2.0 + Mt2/(2.0*WminusTarget*Xminus);
    G4double E  =  WminusTarget*Xminus/2.0 + Mt2/(2.0*WminusTarget*Xminus);
    tmp.setPz( Pz ); 
    tmp.setE( E );
    tmp.transform( toLab );
    aNucleon->SetMomentum( tmp );
    G4VSplitableHadron* targetSplitable = aNucleon->GetSplitableHadron();
    targetSplitable->Set4Momentum( tmp );
  }

  G4double Mt2TrResidual = sqr( TargetResidualMass ) + sqr( TargetResidual3Momentum.x() ) +
                           sqr( TargetResidual3Momentum.y() );

  #ifdef debugPutOnMassShell
  G4cout << "WminusTarget TargetResidual3Momentum.z() " << WminusTarget
         << " " << TargetResidual3Momentum.z() << G4endl;
  #endif

  G4double PzTrResidual = 0.0;
  G4double ETrResidual  = 0.0;
  if ( TargetResidualMassNumber != 0 ) {
    PzTrResidual = -WminusTarget*TargetResidual3Momentum.z()/2.0 + 
                    Mt2TrResidual/( 2.0*WminusTarget*TargetResidual3Momentum.z() );
    ETrResidual  =  WminusTarget*TargetResidual3Momentum.z()/2.0 + 
                    Mt2TrResidual/( 2.0*WminusTarget*TargetResidual3Momentum.z() );
  }

  TargetResidual4Momentum.setPx( TargetResidual3Momentum.x() );
  TargetResidual4Momentum.setPy( TargetResidual3Momentum.y() );
  TargetResidual4Momentum.setPz( PzTrResidual ); 
  TargetResidual4Momentum.setE( ETrResidual );

  #ifdef debugPutOnMassShell
  G4cout << "Target Residual4Momentum in CMS" << TargetResidual4Momentum << G4endl;
  #endif

  TargetResidual4Momentum.transform( toLab );

  #ifdef debugPutOnMassShell
  G4cout << "Target Residual4Momentum in Lab " << TargetResidual4Momentum << G4endl
         << "To continue enter - 1, to break - ^C" << G4endl;
  G4int Uzhi; G4cin >> Uzhi;
  #endif

  return true;
}  


//============================================================================

G4bool G4FTFModel::ExciteParticipants() {

  #ifdef debugBuildString
  G4cout << "G4FTFModel::ExciteParticipants() " << G4endl;
  #endif

  G4bool Successfull( true );  
  G4int MaxNumOfInelCollisions = G4int( theParameters->GetMaxNumberOfCollisions() );
  if ( MaxNumOfInelCollisions > 0 ) {  //  Plab > Pbound, normal application of FTF is possible
    G4double ProbMaxNumber = theParameters->GetMaxNumberOfCollisions() - MaxNumOfInelCollisions;
    if ( G4UniformRand() < ProbMaxNumber ) MaxNumOfInelCollisions++;
  } else { 
    // Plab < Pbound, normal application of FTF is impossible,low energy corrections applied
    MaxNumOfInelCollisions = 1;
  }

  #ifdef debugBuildString
  G4cout << "MaxNumOfInelCollisions MaxNumOfInelCollisions " << MaxNumOfInelCollisions << G4endl;
  #endif

  G4int CurrentInteraction( 0 );
  theParticipants.StartLoop();

  while ( theParticipants.Next() ) {   

    CurrentInteraction++;
    const G4InteractionContent& collision = theParticipants.GetInteraction();
    G4VSplitableHadron* projectile = collision.GetProjectile();
    G4Nucleon* ProjectileNucleon = collision.GetProjectileNucleon();
    G4VSplitableHadron* target = collision.GetTarget();
    G4Nucleon* TargetNucleon = collision.GetTargetNucleon();

    #ifdef debugBuildString
    G4cout << G4endl << "Interaction # Status    " << CurrentInteraction << " " 
           << collision.GetStatus() << G4endl << "Pr* Tr* " << projectile << " "
           << target << G4endl << "projectile->GetStatus target->GetStatus "
           << projectile->GetStatus() << " " << target->GetStatus() << G4endl
           << "projectile->GetSoftC  target->GetSoftC  " << projectile->GetSoftCollisionCount()
           << " " << target->GetSoftCollisionCount() << G4endl;
    #endif

    if ( collision.GetStatus() ) {
      if ( G4UniformRand() < theParameters->GetProbabilityOfElasticScatt() ) { 
        // Elastic scattering

        #ifdef debugBuildString
        G4cout << "Elastic scattering" << G4endl;
        #endif

        if ( ! HighEnergyInter ) {
          G4bool Annihilation = false;
          G4bool Result = AdjustNucleons( projectile, ProjectileNucleon, target,
                                          TargetNucleon, Annihilation );
          if ( ! Result ) continue;
         } 
         Successfull = theElastic->ElasticScattering( projectile, target, theParameters ) 
                       ||  Successfull; 
      } else if ( G4UniformRand() > theParameters->GetProbabilityOfAnnihilation() ) { 
        // Inelastic scattering

        #ifdef debugBuildString
        G4cout << "Inelastic interaction" << G4endl
               << "MaxNumOfInelCollisions " << MaxNumOfInelCollisions << G4endl;
        #endif

        if ( ! HighEnergyInter ) {
          G4bool Annihilation = false;
          G4bool Result = AdjustNucleons( projectile, ProjectileNucleon, target,
                                          TargetNucleon, Annihilation );
          if ( ! Result ) continue;
        } 
        if ( G4UniformRand() < 
             ( 1.0 - projectile->GetSoftCollisionCount() / MaxNumOfInelCollisions ) ) {
          //if ( ! HighEnergyInter ) {
          //  G4bool Annihilation = false;
          //  G4bool Result = AdjustNucleons( projectile, ProjectileNucleon, target,
          //                                  TargetNucleon, Annihilation );
          //  if ( ! Result ) continue;
          //} 
          if (theExcitation->ExciteParticipants( projectile, target, theParameters, theElastic )){

            #ifdef debugBuildString
            G4cout << "FTF excitation Successfull " << G4endl;
            // G4cout << "After  pro " << projectile->Get4Momentum() << " " 
            //        << projectile->Get4Momentum().mag() << G4endl
            //        << "After  tar " << target->Get4Momentum() << " "
            //        << target->Get4Momentum().mag() << G4endl;
            #endif

          } else {

            #ifdef debugBuildString
            G4cout << "FTF excitation Non Successfull -> Elastic scattering " 
                   << Successfull << G4endl;
            #endif

            Successfull = theElastic->ElasticScattering( projectile, target, theParameters )
                          ||  Successfull; 
          }
        } else { // The inelastic interactition was rejected -> elastic scattering

          #ifdef debugBuildString
          G4cout << "Elastic scat. at rejection inelastic scattering" << G4endl;
          #endif

          //if ( ! HighEnergyInter ) {
          //  G4bool Annihilation = false;
          //  G4bool Result = AdjustNucleons( projectile, ProjectileNucleon, target,
          //                                  TargetNucleon, Annihilation );
          //  if ( ! Result) continue;
          //} 
          Successfull = theElastic->ElasticScattering( projectile, target, theParameters )
                        ||  Successfull;
        } 
      } else {  // Annihilation

        #ifdef debugBuildString
        G4cout << "Annihilation" << G4endl;
        #endif
       
        // Skipping possible interactions of the annihilated nucleons 
        while ( theParticipants.Next() ) {   
          G4InteractionContent& acollision = theParticipants.GetInteraction();
          G4VSplitableHadron* NextProjectileNucleon = acollision.GetProjectile();
          G4VSplitableHadron* NextTargetNucleon = acollision.GetTarget();
          if ( projectile == NextProjectileNucleon  ||  target == NextTargetNucleon ) {
            acollision.SetStatus( 0 );
          }
        }

        // Return to the annihilation
        theParticipants.StartLoop(); 
        for ( G4int I = 0; I < CurrentInteraction; I++ ) theParticipants.Next();

        // At last, annihilation
        if ( ! HighEnergyInter ) {
          G4bool Annihilation = true;
          G4bool Result = AdjustNucleons( projectile, ProjectileNucleon, target,
                                          TargetNucleon, Annihilation );
          if ( ! Result ) continue;
        }             
        G4VSplitableHadron* AdditionalString = 0;
        if ( theAnnihilation->Annihilate( projectile, target, AdditionalString, theParameters ) ){
          Successfull = Successfull  ||  true;

          #ifdef debugBuildString
          G4cout << "Annihilation successfull. " << "*AdditionalString " 
                 << AdditionalString << G4endl;
          //G4cout << "After  pro " << projectile->Get4Momentum() << G4endl;
          //G4cout << "After  tar " << target->Get4Momentum() << G4endl;
          #endif

          if ( AdditionalString != 0 ) theAdditionalString.push_back( AdditionalString );
        }
      } 
    }

    #ifdef debugBuildString
    G4cout << "----------------------------- Final properties " << G4endl
           << "projectile->GetStatus target->GetStatus " << projectile->GetStatus() 
           << " " << target->GetStatus() << G4endl << "projectile->GetSoftC  target->GetSoftC  "
           << projectile->GetSoftCollisionCount() << " " << target->GetSoftCollisionCount()
           << G4endl << "ExciteParticipants() Successfull? " << Successfull << G4endl;
    #endif

  }  // end of while ( theParticipants.Next() )

  return Successfull;
}


//============================================================================

G4bool G4FTFModel::AdjustNucleons( G4VSplitableHadron* SelectedAntiBaryon, 
                                   G4Nucleon* ProjectileNucleon,
                                   G4VSplitableHadron* SelectedTargetNucleon,
                                   G4Nucleon* TargetNucleon,
                                   G4bool Annihilation ) {

  #ifdef debugAdjust
  G4cout << "AdjustNucleons ---------------------------------------" << G4endl
         << "Proj is nucleus? " << GetProjectileNucleus() << G4endl
         << "Proj 4mom " << SelectedAntiBaryon->Get4Momentum() << G4endl
         << "Targ 4mom " << SelectedTargetNucleon->Get4Momentum() << G4endl
         << "Pr ResidualMassNumber Pr ResidualCharge Pr ResidualExcitationEnergy "
         << ProjectileResidualMassNumber << " " << ProjectileResidualCharge << " "
         << ProjectileResidualExcitationEnergy << G4endl
         << "Tr ResidualMassNumber Tr ResidualCharge Tr ResidualExcitationEnergy  "
         << TargetResidualMassNumber << " " << TargetResidualCharge << " "
         << TargetResidualExcitationEnergy << G4endl
         << "Collis. pr tr " << SelectedAntiBaryon->GetSoftCollisionCount()
         << SelectedTargetNucleon->GetSoftCollisionCount() << G4endl;
  #endif

  if ( SelectedAntiBaryon->GetSoftCollisionCount() != 0  && 
       SelectedTargetNucleon->GetSoftCollisionCount() != 0 ) {
    return true; // Selected hadrons were adjusted before.   
  }

  // Ascribing of the involved nucleons Pt and X 
  G4double Dcor = theParameters->GetDofNuclearDestruction();

  G4double DcorP( 0.0 ), DcorT( 0.0 );
  if ( ProjectileResidualMassNumber != 0 ) DcorP = Dcor / G4double(ProjectileResidualMassNumber);
  if ( TargetResidualMassNumber != 0 )     DcorT = Dcor / G4double(TargetResidualMassNumber);

  G4double AveragePt2 = theParameters->GetPt2ofNuclearDestruction();
  G4double maxPtSquare = theParameters->GetMaxPt2ofNuclearDestruction();
  G4double ExcitationEnergyPerWoundedNucleon = 
             theParameters->GetExcitationEnergyPerWoundedNucleon();

  if ( ( ! GetProjectileNucleus()  &&
         SelectedAntiBaryon->GetSoftCollisionCount() == 0  &&
         SelectedTargetNucleon->GetSoftCollisionCount() == 0 )
      ||
       ( SelectedAntiBaryon->GetSoftCollisionCount() != 0  && 
         SelectedTargetNucleon->GetSoftCollisionCount() == 0 ) ) {
    // The case of hadron-nucleus interactions, or
    // the case when projectile nuclear nucleon participated in
    // a collision, but target nucleon did not participate.

    #ifdef debugAdjust 
    G4cout << "case 1, hA prcol=0 trcol=0, AA prcol#0 trcol=0" << G4endl;
    #endif

    if ( TargetResidualMassNumber < 1 ) {
      return false;
    }

    if ( SelectedAntiBaryon->Get4Momentum().rapidity() < TargetResidual4Momentum.rapidity() ) {
      return false;
    }

    if ( TargetResidualMassNumber == 1 ) {
      TargetResidualMassNumber       = 0;
      TargetResidualCharge           = 0;
      TargetResidualExcitationEnergy = 0.0;
      SelectedTargetNucleon->Set4Momentum( TargetResidual4Momentum );
      TargetResidual4Momentum = G4LorentzVector( 0.0, 0.0, 0.0, 0.0 );
      return true;
    }

    G4LorentzVector Psum = SelectedAntiBaryon->Get4Momentum() + TargetResidual4Momentum;
 
    #ifdef debugAdjust
    G4cout << "Targ res Init " << TargetResidual4Momentum << G4endl;
    #endif

    // Transform momenta to cms and then rotate parallel to z axis;
    G4LorentzRotation toCms( -1*Psum.boostVector() );
    G4LorentzVector Pprojectile = SelectedAntiBaryon->Get4Momentum();
    G4LorentzVector Ptmp = toCms * Pprojectile;
    toCms.rotateZ( -1*Ptmp.phi() );
    toCms.rotateY( -1*Ptmp.theta() );
    Pprojectile.transform( toCms );
    G4LorentzRotation toLab( toCms.inverse() );

    G4LorentzVector Ptarget( 0.0, 0.0, 0.0, 0.0 );

    G4double SqrtS = Psum.mag();
    G4double S = sqr( SqrtS );

    G4int TResidualMassNumber = TargetResidualMassNumber - 1;
    G4int TResidualCharge = TargetResidualCharge - 
                            G4int( TargetNucleon->GetDefinition()->GetPDGCharge() );
    G4double TResidualExcitationEnergy = TargetResidualExcitationEnergy + 
                                         ExcitationEnergyPerWoundedNucleon;
    if ( TResidualMassNumber <= 1 ) {
      TResidualExcitationEnergy = 0.0;
    }

    G4double TResidualMass( 0.0 );
    if ( TResidualMassNumber != 0 ) {
      TResidualMass = G4ParticleTable::GetParticleTable()->GetIonTable()
                          ->GetIonMass( TResidualCharge, TResidualMassNumber );
    }
 
    G4double TNucleonMass = TargetNucleon->GetDefinition()->GetPDGMass();
    G4double SumMasses = SelectedAntiBaryon->Get4Momentum().mag() + TNucleonMass + TResidualMass;

    G4bool Stopping = false;

    #ifdef debugAdjust
    G4cout << "Annihilation " << Annihilation << G4endl;
    #endif

    if ( ! Annihilation ) {
      if ( SqrtS < SumMasses ) {
        return false;
      } 
      if ( SqrtS < SumMasses + TResidualExcitationEnergy ) { 

        #ifdef debugAdjust
        G4cout << "TResidualExcitationEnergy " << TResidualExcitationEnergy << G4endl;
        #endif

        TResidualExcitationEnergy = SqrtS - SumMasses;

        #ifdef debugAdjust
        G4cout << "TResidualExcitationEnergy " << TResidualExcitationEnergy << G4endl;
        #endif

        Stopping = true; 
        return false;
      }
    }

    if ( Annihilation ) {

      #ifdef debugAdjust
      G4cout << "SqrtS < SumMasses - TNucleonMass " << SqrtS << " " 
             << SumMasses - TNucleonMass << G4endl;
      #endif

      if ( SqrtS < SumMasses - TNucleonMass ) {
        return false;
      } 

      #ifdef debugAdjust
      G4cout << "SqrtS < SumMasses " << SqrtS << " " << SumMasses << G4endl;
      #endif

      if ( SqrtS < SumMasses ) { 
        TNucleonMass = SqrtS - (SumMasses - TNucleonMass) - TResidualExcitationEnergy; 

        #ifdef debugAdjust
        G4cout << "TNucleonMass " << TNucleonMass << G4endl;
        #endif

        SumMasses = SqrtS - TResidualExcitationEnergy; // Uzhi Feb. 2013
        //TResidualExcitationEnergy =0.0;              // Uzhi Feb. 2013
        Stopping = true;
      }

      #ifdef debugAdjust
      G4cout << "SqrtS < SumMasses " << SqrtS << " " << SumMasses << G4endl;
      #endif

      if ( SqrtS < SumMasses + TResidualExcitationEnergy ) { 
        TResidualExcitationEnergy = SqrtS - SumMasses; 
        Stopping = true;
      }
    }

    #ifdef debugAdjust
    G4cout << "Stopping " << Stopping << G4endl;
    #endif

    if ( Stopping ) {
      // All 3-momenta of particles = 0
      // New projectile
      Ptmp.setPx( 0.0 ); Ptmp.setPy( 0.0 ); Ptmp.setPz( 0.0 );
      Ptmp.setE( SelectedAntiBaryon->Get4Momentum().mag() );

      #ifdef debugAdjust
      G4cout << "Proj stop " << Ptmp << G4endl;
      #endif

      Pprojectile = Ptmp; Pprojectile.transform( toLab );
      SelectedAntiBaryon->Set4Momentum( Pprojectile );

      // New target nucleon
      Ptmp.setE( TNucleonMass );

      #ifdef debugAdjust
      G4cout << "Targ stop " << Ptmp << G4endl;
      #endif

      Ptarget = Ptmp; Ptarget.transform( toLab );
      SelectedTargetNucleon->Set4Momentum( Ptarget );

      // New target residual
      TargetResidualMassNumber       = TResidualMassNumber;
      TargetResidualCharge           = TResidualCharge;
      TargetResidualExcitationEnergy = TResidualExcitationEnergy;

      Ptmp.setE( TResidualMass + TargetResidualExcitationEnergy ); 

      #ifdef debugAdjust
      G4cout << "Resi stop " << Ptmp << G4endl;
      #endif

      Ptmp.transform( toLab );          
      TargetResidual4Momentum = Ptmp;

      #ifdef debugAdjust
      G4cout << Pprojectile << G4endl << Ptarget << G4endl << TargetResidual4Momentum << G4endl;
      #endif

      return true;
    }
        
    G4double Mprojectile  = Pprojectile.mag();
    G4double M2projectile = Pprojectile.mag2();
    G4double WplusProjectile( 0.0 );

    G4LorentzVector TResidual4Momentum = toCms * TargetResidual4Momentum;
    G4double YtargetNucleus = TResidual4Momentum.rapidity();

    TResidualMass += TResidualExcitationEnergy;
    G4double M2target( 0.0 );
    G4double WminusTarget( 0.0 );

    G4ThreeVector PtNucleon( 0.0, 0.0, 0.0 );
    G4double XminusNucleon( 0.0 );
    G4ThreeVector PtResidual( 0.0, 0.0, 0.0 );
    G4double XminusResidual( 0.0 );

    G4int NumberOfTries( 0 );
    G4double ScaleFactor( 1.0 );
    G4bool OuterSuccess( true );

    do { // while ( ! OuterSuccess )
      OuterSuccess = true;

      do { // while ( SqrtS < Mprojectile + std::sqrt( M2target) )

        NumberOfTries++;

        if ( NumberOfTries == 100*(NumberOfTries/100) ) { 
          // At large number of tries it would be better to reduce the values
          ScaleFactor /= 2.0;
          DcorT      *= ScaleFactor;
          AveragePt2 *= ScaleFactor;
        }

        //if ( TargetResidualMassNumber > 1 ) {
        //  PtNucleon = GaussianPt( AveragePt2, maxPtSquare );
        //} else {
        //  PtNucleon = G4ThreeVector( 0.0, 0.0, 0.0 );
        //}
        //PtResidual = -PtNucleon;

        G4bool InerSuccess = true;
        if ( TargetResidualMassNumber > 1 ) {
          do {
            InerSuccess = true;

            PtNucleon = GaussianPt( AveragePt2, maxPtSquare );
            PtResidual = -PtNucleon;

            G4double Mtarget = std::sqrt( sqr( TNucleonMass ) + PtNucleon.mag2() ) + 
                               std::sqrt( sqr( TResidualMass ) + PtResidual.mag2() );
            if ( SqrtS < Mprojectile + Mtarget ) {
              InerSuccess = false; 
              continue;
            }

            G4ThreeVector tmpX = GaussianPt( DcorT*DcorT, 1.0 );
            G4double Xcenter = std::sqrt( sqr( TNucleonMass )  + PtNucleon.mag2() ) / Mtarget;
            XminusNucleon = Xcenter + tmpX.x();
            if ( XminusNucleon <= 0.0  ||  XminusNucleon >= 1.0 ) {
              InerSuccess = false; 
              continue;
            }

            XminusResidual = 1.0 - XminusNucleon;
          } while ( ! InerSuccess );
        } else {
          XminusNucleon  = 1.0;
          XminusResidual = 1.0;  // It must be 0, but in the case calculation of Pz,
                                 // E is problematic.
        }

        M2target = ( sqr( TNucleonMass ) + PtNucleon.mag2() ) / XminusNucleon + 
                   ( sqr( TResidualMass ) + PtResidual.mag2() ) / XminusResidual;

      } while ( SqrtS < Mprojectile + std::sqrt( M2target) );

      G4double DecayMomentum2 = sqr( S ) + sqr( M2projectile ) + sqr( M2target )
                                - 2.0*S*M2projectile - 2.0*S*M2target - 2.0*M2projectile*M2target;

      WminusTarget = ( S - M2projectile + M2target + std::sqrt( DecayMomentum2 ) ) / 2.0 / SqrtS;
      WplusProjectile = SqrtS - M2target / WminusTarget;

      G4double Pzprojectile = WplusProjectile/2.0 - M2projectile/2.0/WplusProjectile;
      G4double Eprojectile  = WplusProjectile/2.0 + M2projectile/2.0/WplusProjectile;
      G4double Yprojectile  = 0.5 * std::log( (Eprojectile + Pzprojectile) /
                                              (Eprojectile - Pzprojectile) );

      #ifdef debugAdjust
      G4cout << "DecayMomentum2 " << DecayMomentum2 << G4endl
             << "WminusTarget WplusProjectile " << WminusTarget << " " << WplusProjectile
             << G4endl << "Yprojectile " << Yprojectile << G4endl;
      #endif

      G4double Mt2 = sqr( TNucleonMass ) + PtNucleon.mag2();
      G4double Pz = -WminusTarget*XminusNucleon/2.0 + Mt2/(2.0*WminusTarget*XminusNucleon);
      G4double E  =  WminusTarget*XminusNucleon/2.0 + Mt2/(2.0*WminusTarget*XminusNucleon);
      G4double YtargetNucleon = 0.5 * std::log( (E + Pz)/(E - Pz) ); 

      #ifdef debugAdjust
      G4cout << "YtN Ytr YtN-Ytr " << " " << YtargetNucleon << " " << YtargetNucleus << " "
             << YtargetNucleon - YtargetNucleus << G4endl
             << "YtN Ypr YtN-Ypr " << " " << YtargetNucleon << " " << Yprojectile
             << " " << YtargetNucleon - Yprojectile << G4endl;
      #endif

      if ( std::abs( YtargetNucleon - YtargetNucleus ) > 2 || Yprojectile < YtargetNucleon ) {
        OuterSuccess = false; 
        continue;
      } 

    } while ( ! OuterSuccess );

    G4double Pzprojectile = WplusProjectile/2.0 - M2projectile/2.0/WplusProjectile;
    G4double Eprojectile  = WplusProjectile/2.0 + M2projectile/2.0/WplusProjectile;
    Pprojectile.setPz( Pzprojectile ); Pprojectile.setE( Eprojectile );

    #ifdef debugAdjust
    G4cout << "Proj after in CMS " << Pprojectile << G4endl;
    #endif

    Pprojectile.transform( toLab ); // The work with the projectile is finished at the moment.

    SelectedAntiBaryon->Set4Momentum( Pprojectile );

    #ifdef debugAdjust
    G4cout << "New proj4M " << Pprojectile << G4endl;
    #endif

    G4double Mt2 = sqr( TNucleonMass ) + PtNucleon.mag2();
    G4double Pz = -WminusTarget*XminusNucleon/2.0 + Mt2/(2.0*WminusTarget*XminusNucleon);
    G4double E  =  WminusTarget*XminusNucleon/2.0 + Mt2/(2.0*WminusTarget*XminusNucleon);

    Ptarget.setPx( PtNucleon.x() ); Ptarget.setPy( PtNucleon.y() );
    Ptarget.setPz( Pz );            Ptarget.setE( E ); 
    Ptarget.transform( toLab );
    SelectedTargetNucleon->Set4Momentum( Ptarget );

    #ifdef debugAdjust
    G4cout << "New targ4M " << Ptarget << G4endl;
    #endif

    // New target residual
    TargetResidualMassNumber       = TResidualMassNumber;
    TargetResidualCharge           = TResidualCharge;
    TargetResidualExcitationEnergy = TResidualExcitationEnergy;

    #ifdef debugAdjust
    G4cout << "TargetResidualMassNumber TargetResidualCharge TargetResidualExcitationEnergy "
           << TargetResidualMassNumber << " " << TargetResidualCharge << " "
           << TargetResidualExcitationEnergy << G4endl;
    #endif

    if ( TargetResidualMassNumber != 0 ) {
      Mt2 = sqr( TResidualMass ) + PtResidual.mag2();
      Pz = -WminusTarget*XminusResidual/2.0 + Mt2/(2.0*WminusTarget*XminusResidual);
      E =   WminusTarget*XminusResidual/2.0 + Mt2/(2.0*WminusTarget*XminusResidual);

      TargetResidual4Momentum.setPx( PtResidual.x() );
      TargetResidual4Momentum.setPy( PtResidual.y() );
      TargetResidual4Momentum.setPz( Pz );
      TargetResidual4Momentum.setE( E );

      #ifdef debugAdjust
      G4cout << "New Residu " << TargetResidual4Momentum << " CMS" << G4endl;
      #endif

      TargetResidual4Momentum.transform( toLab );

      #ifdef debugAdjust
      G4cout << "New Residu " << TargetResidual4Momentum << " Lab" << G4endl;
      #endif

    } else {
      TargetResidual4Momentum = G4LorentzVector( 0.0, 0.0, 0.0, 0.0 );
    }
    return true;

  } else if ( SelectedAntiBaryon->GetSoftCollisionCount() == 0  && 
              SelectedTargetNucleon->GetSoftCollisionCount() != 0 ) {
    // It is assumed that in the case there is ProjectileResidualNucleus

    #ifdef debugAdjust
    G4cout << "case 2,  prcol=0 trcol#0" << G4endl;
    #endif

    if ( ProjectileResidualMassNumber < 1 ) return false;

    if ( ProjectileResidual4Momentum.rapidity() <= 
         SelectedTargetNucleon->Get4Momentum().rapidity() ) {
      return false;
    }

    if ( ProjectileResidualMassNumber == 1 ) {
      ProjectileResidualMassNumber       = 0;
      ProjectileResidualCharge           = 0;
      ProjectileResidualExcitationEnergy = 0.0;
      SelectedAntiBaryon->Set4Momentum( ProjectileResidual4Momentum );        
      ProjectileResidual4Momentum = G4LorentzVector( 0.0, 0.0, 0.0, 0.0 );
      return true;
    }

    G4LorentzVector Psum = ProjectileResidual4Momentum + SelectedTargetNucleon->Get4Momentum();

    // Transform momenta to cms and then rotate parallel to z axis;
    G4LorentzRotation toCms( -1*Psum.boostVector() );
    G4LorentzVector Pprojectile = ProjectileResidual4Momentum;
    G4LorentzVector Ptmp = toCms * Pprojectile;
    toCms.rotateZ( -1*Ptmp.phi() );
    toCms.rotateY( -1*Ptmp.theta() );
    G4LorentzRotation toLab( toCms.inverse() );
    G4LorentzVector Ptarget = toCms * SelectedTargetNucleon->Get4Momentum();
    Pprojectile.transform( toCms );

    G4double SqrtS = Psum.mag();
    G4double S = sqr( SqrtS );

    G4int TResidualMassNumber = ProjectileResidualMassNumber - 1;
    G4int TResidualCharge = ProjectileResidualCharge 
                          - std::abs( G4int(ProjectileNucleon->GetDefinition()->GetPDGCharge()) );
    G4double TResidualExcitationEnergy = ProjectileResidualExcitationEnergy + 
                                         ExcitationEnergyPerWoundedNucleon;
    if ( TResidualMassNumber <= 1 ) {
      TResidualExcitationEnergy = 0.0;
    }

    G4double TResidualMass( 0.0 );
    if ( TResidualMassNumber != 0 ) {
      TResidualMass = G4ParticleTable::GetParticleTable()->GetIonTable()
                        ->GetIonMass( TResidualCharge , TResidualMassNumber );
    }
 
    G4double TNucleonMass = ProjectileNucleon->GetDefinition()->GetPDGMass();

    G4double SumMasses = SelectedTargetNucleon->Get4Momentum().mag() +
                         TNucleonMass + TResidualMass;

    #ifdef debugAdjust
    G4cout << "SelectedTN.mag() PNMass + PResidualMass " 
           << SelectedTargetNucleon->Get4Momentum().mag() << " " 
           << TNucleonMass << " " << TResidualMass << G4endl;
    #endif

    G4bool Stopping = false;

    if ( ! Annihilation ) {
      if ( SqrtS < SumMasses ) {
        return false;
      } 
      if ( SqrtS < SumMasses + TResidualExcitationEnergy ) { 
        TResidualExcitationEnergy = SqrtS - SumMasses;
        Stopping = true; 
        return false;
      }
    }

    if ( Annihilation ) {
      if ( SqrtS < SumMasses - TNucleonMass ) {
        return false;
      } 
      if ( SqrtS < SumMasses ) { 
        TNucleonMass = SqrtS - (SumMasses - TNucleonMass);
        SumMasses = SqrtS;
        TResidualExcitationEnergy = 0.0;
        Stopping = true;
      }

      if ( SqrtS < SumMasses + TResidualExcitationEnergy ) { 
        TResidualExcitationEnergy = SqrtS - SumMasses; 
        Stopping=true;
      }
    }

    #ifdef debugAdjust
    G4cout << "Stopping " << Stopping << G4endl;
    #endif

    if ( Stopping ) { 
      // All 3-momenta of particles = 0
      // New target nucleon
      Ptmp.setPx( 0.0 ); Ptmp.setPy( 0.0 ); Ptmp.setPz( 0.0 );
      Ptmp.setE( SelectedTargetNucleon->Get4Momentum().mag() );
      Ptarget = Ptmp; Ptarget.transform( toLab );
      SelectedTargetNucleon->Set4Momentum( Ptarget );

      // New projectile nucleon
      Ptmp.setE( TNucleonMass );
      Pprojectile = Ptmp; Pprojectile.transform( toLab );
      SelectedAntiBaryon->Set4Momentum( Pprojectile );

      // New projectile residual
      ProjectileResidualMassNumber       = TResidualMassNumber;
      ProjectileResidualCharge           = TResidualCharge;
      ProjectileResidualExcitationEnergy = TResidualExcitationEnergy;

      Ptmp.setE( TResidualMass + ProjectileResidualExcitationEnergy ); 
      Ptmp.transform( toLab );          
      ProjectileResidual4Momentum = Ptmp;

      return true;
    }

    G4double Mtarget  = Ptarget.mag();
    G4double M2target = Ptarget.mag2();

    G4LorentzVector TResidual4Momentum = toCms * ProjectileResidual4Momentum;
    G4double YprojectileNucleus = TResidual4Momentum.rapidity();

    TResidualMass += TResidualExcitationEnergy;

    G4double M2projectile( 0.0 );
    G4double WminusTarget( 0.0 );
    G4double WplusProjectile( 0.0 );
    G4ThreeVector PtNucleon( 0.0, 0.0, 0.0 );
    G4double XplusNucleon( 0.0 );
    G4ThreeVector PtResidual( 0.0, 0.0, 0.0 );
    G4double XplusResidual( 0.0 );
    G4int NumberOfTries( 0 );
    G4double ScaleFactor( 1.0 );
    G4bool OuterSuccess( true );

    do { // while ( ! OuterSuccess )
  
      OuterSuccess = true;

      do { // while ( SqrtS < Mtarget + std::sqrt( M2projectile ) )

        NumberOfTries++;

        if ( NumberOfTries == 100*(NumberOfTries/100) ) { 
          // At large number of tries it would be better to reduce the values
          ScaleFactor /= 2.0;
          DcorP      *= ScaleFactor;
          AveragePt2 *= ScaleFactor;
        }

        #ifdef debugAdjust
        G4cout << "ProjectileResidualMassNumber " << ProjectileResidualMassNumber << G4endl;
        #endif

        if ( ProjectileResidualMassNumber > 1 ) {
          PtNucleon = GaussianPt( AveragePt2, maxPtSquare );
        } else {
          PtNucleon = G4ThreeVector( 0.0, 0.0, 0.0 );
        }
        PtResidual = -PtNucleon;

        G4double Mprojectile = std::sqrt( sqr( TNucleonMass )  + PtNucleon.mag2() ) + 
                               std::sqrt( sqr( TResidualMass ) + PtResidual.mag2() );

        #ifdef debugAdjust
        G4cout << "SqrtS < Mtarget + Mprojectile " << SqrtS << "  " << Mtarget
               << " " << Mprojectile << " " << Mtarget + Mprojectile << G4endl;
        #endif

        M2projectile = sqr( Mprojectile );  // Uzhi 31.08.13
        if ( SqrtS < Mtarget + Mprojectile ) {
          OuterSuccess = false; 
          continue;
        }

        G4double Xcenter = std::sqrt( sqr( TNucleonMass ) + PtNucleon.mag2() ) / Mprojectile;

        G4bool InerSuccess = true;
        if ( ProjectileResidualMassNumber > 1 ) {
          do {
            InerSuccess = true;
            G4ThreeVector tmpX = GaussianPt( DcorP*DcorP, 1.0 );
            XplusNucleon = Xcenter + tmpX.x();
            if ( XplusNucleon <= 0.0  ||  XplusNucleon >= 1.0 ) { 
              InerSuccess = false; 
              continue;
            }
            XplusResidual = 1.0 - XplusNucleon;
          } while ( ! InerSuccess );
        } else {
          XplusNucleon  = 1.0;
          XplusResidual = 1.0; // It must be 0, but in the case determination
                               // of Pz and E will be problematic.
        }

        #ifdef debugAdjust
        G4cout << "TNucleonMass PtNucleon XplusNucleon " << TNucleonMass << " " << PtNucleon
               << " " << XplusNucleon << G4endl
               << "TResidualMass PtResidual XplusResidual " << TResidualMass << " " << PtResidual
               << "  " << XplusResidual << G4endl;
        #endif

        M2projectile = ( sqr( TNucleonMass )  + PtNucleon.mag2() ) / XplusNucleon + 
                       ( sqr( TResidualMass ) + PtResidual.mag2() ) / XplusResidual;

        #ifdef debugAdjust
        G4cout << "SqrtS < Mtarget + std::sqrt(M2projectile) " << SqrtS << "  " << Mtarget
               << " " << std::sqrt( M2projectile ) << " " << Mtarget + std::sqrt( M2projectile ) 
               << G4endl;
        #endif

      } while ( SqrtS < Mtarget + std::sqrt( M2projectile ) );

      G4double DecayMomentum2 = sqr( S ) + sqr( M2projectile ) + sqr( M2target )
                                - 2.0*S*M2projectile - 2.0*S*M2target - 2.0*M2projectile*M2target;

      WplusProjectile = ( S + M2projectile - M2target + std::sqrt( DecayMomentum2 ) )/2.0/SqrtS;
      WminusTarget = SqrtS - M2projectile/WplusProjectile;

      G4double Pztarget = -WminusTarget/2.0 + M2target/2.0/WminusTarget;
      G4double Etarget =   WminusTarget/2.0 + M2target/2.0/WminusTarget;
      G4double Ytarget = 0.5 * std::log( (Etarget + Pztarget)/(Etarget - Pztarget) );

      #ifdef debugAdjust
      G4cout << "DecayMomentum2 " << DecayMomentum2 << G4endl
             << "WminusTarget WplusProjectile " << WminusTarget << " " << WplusProjectile 
             << G4endl << "YtargetNucleon " << Ytarget << G4endl;
      #endif

      G4double Mt2 = sqr( TNucleonMass ) + PtNucleon.mag2();
      G4double Pz = WplusProjectile*XplusNucleon/2.0 - Mt2/(2.0*WplusProjectile*XplusNucleon);
      G4double E =  WplusProjectile*XplusNucleon/2.0 + Mt2/(2.0*WplusProjectile*XplusNucleon);
      G4double YprojectileNucleon = 0.5 * std::log( (E + Pz)/(E - Pz) ); 

      #ifdef debugAdjust
      G4cout << "YpN Ypr YpN-Ypr " << " " << YprojectileNucleon << " " << YprojectileNucleus
             << " " << YprojectileNucleon - YprojectileNucleus << G4endl
             << "YpN Ytr YpN-Ytr " << " " << YprojectileNucleon << " " << Ytarget 
             << " " << YprojectileNucleon - Ytarget << G4endl;
      #endif

      if ( std::abs( YprojectileNucleon - YprojectileNucleus ) > 2  ||
           Ytarget  > YprojectileNucleon ) {
        OuterSuccess = false; 
        continue;
      }
 
    } while ( ! OuterSuccess );

    // New target
    G4double Pztarget = -WminusTarget/2.0 + M2target/2.0/WminusTarget;
    G4double Etarget  =  WminusTarget/2.0 + M2target/2.0/WminusTarget;
    Ptarget.setPz( Pztarget ); Ptarget.setE( Etarget );
    Ptarget.transform( toLab ); // The work with the target nucleon is finished at the moment.
    SelectedTargetNucleon->Set4Momentum( Ptarget );

    #ifdef debugAdjust
    G4cout << "Targ after in Lab " << Ptarget << G4endl;
    #endif

    // New projectile
    G4double Mt2 = sqr( TNucleonMass ) + PtNucleon.mag2();
    G4double Pz = WplusProjectile*XplusNucleon/2.0 - Mt2/(2.0*WplusProjectile*XplusNucleon);
    G4double E =  WplusProjectile*XplusNucleon/2.0 + Mt2/(2.0*WplusProjectile*XplusNucleon);
    Pprojectile.setPx( PtNucleon.x() ); Pprojectile.setPy( PtNucleon.y() );
    Pprojectile.setPz( Pz );            Pprojectile.setE( E ); 
    Pprojectile.transform( toLab );
    SelectedAntiBaryon->Set4Momentum( Pprojectile );

    #ifdef debugAdjust
    G4cout << "Proj after in Lab " << Pprojectile << G4endl;
    #endif

    // New projectile residual
    ProjectileResidualMassNumber       = TResidualMassNumber;
    ProjectileResidualCharge           = TResidualCharge;
    ProjectileResidualExcitationEnergy = TResidualExcitationEnergy;

    if ( ProjectileResidualMassNumber != 0 ) {
      Mt2 = sqr( TResidualMass ) + PtResidual.mag2();
      Pz = WplusProjectile*XplusResidual/2.0 - Mt2/(2.0*WplusProjectile*XplusResidual);
      E  = WplusProjectile*XplusResidual/2.0 + Mt2/(2.0*WplusProjectile*XplusResidual);
      ProjectileResidual4Momentum.setPx( PtResidual.x() );
      ProjectileResidual4Momentum.setPy( PtResidual.y() );
      ProjectileResidual4Momentum.setPz( Pz );
      ProjectileResidual4Momentum.setE( E );
      ProjectileResidual4Momentum.transform( toLab );
    } else {
      ProjectileResidual4Momentum = G4LorentzVector( 0.0, 0.0, 0.0, 0.0 );
    }
    return true;

  } else { // if ( SelectedAntiBaryon->GetSoftCollisionCount() == 0 && 
           //      SelectedTargetNucleon->GetSoftCollisionCount() == 0 )

    //    It can be in the case of nucleus-nucleus interaction only!

    #ifdef debugAdjust
    G4cout << "case 3,  prcol=0 trcol=0" << G4endl;
    #endif

    if ( ! GetProjectileNucleus() ) return false;

    #ifdef debugAdjust
    G4cout << "Proj res Init " << ProjectileResidual4Momentum << G4endl
           << "Targ res Init " << TargetResidual4Momentum << G4endl
           << "ProjectileResidualMassNumber ProjectileResidualCharge " 
           << ProjectileResidualMassNumber << " " << ProjectileResidualCharge << G4endl
           << "TargetResidualMassNumber TargetResidualCharge " << TargetResidualMassNumber
           << " " << TargetResidualCharge << G4endl;
    #endif

    G4LorentzVector Psum = ProjectileResidual4Momentum + TargetResidual4Momentum; 

    // Transform momenta to cms and then rotate parallel to z axis;
    G4LorentzRotation toCms( -1*Psum.boostVector() );
    G4LorentzVector Pprojectile = ProjectileResidual4Momentum;
    G4LorentzVector Ptmp = toCms * Pprojectile;
    toCms.rotateZ( -1*Ptmp.phi() );
    toCms.rotateY( -1*Ptmp.theta() );
    G4LorentzRotation toLab( toCms.inverse() );
    Pprojectile.transform( toCms );
    G4LorentzVector Ptarget = toCms * TargetResidual4Momentum;

    G4double SqrtS = Psum.mag();
    G4double S = sqr( SqrtS );

    G4int PResidualMassNumber = ProjectileResidualMassNumber - 1;
    G4int PResidualCharge = ProjectileResidualCharge - 
                            std::abs( G4int(ProjectileNucleon->GetDefinition()->GetPDGCharge()) );
    G4double PResidualExcitationEnergy = ProjectileResidualExcitationEnergy +
                                         ExcitationEnergyPerWoundedNucleon;
    if ( PResidualMassNumber <= 1 ) {
      PResidualExcitationEnergy = 0.0;
    }

    G4double PResidualMass( 0.0 );
    if ( PResidualMassNumber != 0 ) {
      PResidualMass = G4ParticleTable::GetParticleTable()->GetIonTable()
                        ->GetIonMass( PResidualCharge, PResidualMassNumber );
    }
 
    G4double PNucleonMass = ProjectileNucleon->GetDefinition()->GetPDGMass();

    G4int TResidualMassNumber = TargetResidualMassNumber - 1;
    G4int TResidualCharge = TargetResidualCharge - 
                            G4int( TargetNucleon->GetDefinition()->GetPDGCharge() );
    G4double TResidualExcitationEnergy = TargetResidualExcitationEnergy +
                                         ExcitationEnergyPerWoundedNucleon;
    if ( TResidualMassNumber <= 1 ) {
      TResidualExcitationEnergy = 0.0;
    }
    G4double TResidualMass( 0.0 );
    if ( TResidualMassNumber != 0 ) {
      TResidualMass = G4ParticleTable::GetParticleTable()->GetIonTable()
                        ->GetIonMass( TResidualCharge, TResidualMassNumber );
    }

    G4double TNucleonMass = TargetNucleon->GetDefinition()->GetPDGMass();

    G4double SumMasses = PNucleonMass + PResidualMass + TNucleonMass + TResidualMass;

    #ifdef debugAdjust
    G4cout << "PNucleonMass PResidualMass TNucleonMass TResidualMass " << PNucleonMass 
           << " " << PResidualMass << " " << TNucleonMass << " " << TResidualMass << G4endl
           << "PResidualExcitationEnergy " << PResidualExcitationEnergy << G4endl
           << "TResidualExcitationEnergy " << TResidualExcitationEnergy << G4endl;
    #endif

    G4bool Stopping = false;

    if ( ! Annihilation ) {

      #ifdef debugAdjust
      G4cout << "SqrtS < SumMasses " << SqrtS << " " << SumMasses << G4endl;
      #endif

      if ( SqrtS < SumMasses ) {
        return false;
      } 

      #ifdef debugAdjust
      G4cout << "SqrtS < SumMasses + PResidualExcitationEnergy + TResidualExcitationEnergy "
             << SqrtS << " " << SumMasses + PResidualExcitationEnergy + TResidualExcitationEnergy
             << G4endl;
      #endif

      if ( SqrtS < SumMasses + PResidualExcitationEnergy + TResidualExcitationEnergy ) { 
        Stopping = true; 
        //AR-14Aug2013  return false;
        if ( PResidualExcitationEnergy <= 0.0 ) {
          TResidualExcitationEnergy = SqrtS - SumMasses;
        } else if ( TResidualExcitationEnergy <= 0.0 ) {
          PResidualExcitationEnergy = SqrtS - SumMasses;
        } else {
          G4double Fraction = (SqrtS - SumMasses) / 
                              (PResidualExcitationEnergy + TResidualExcitationEnergy);
          PResidualExcitationEnergy *= Fraction;
          TResidualExcitationEnergy *= Fraction;
        }
      }
    }

    #ifdef debugAdjust
    G4cout << "Stopping " << Stopping << G4endl;
    #endif

    if ( Annihilation ) {
      if ( SqrtS < SumMasses - TNucleonMass ) {
        return false;
      } 
      if ( SqrtS < SumMasses ) { 
        Stopping = true;
        TNucleonMass = SqrtS - (SumMasses - TNucleonMass);
        SumMasses = SqrtS;
        TResidualExcitationEnergy = 0.0;
      }
      if ( SqrtS < SumMasses + PResidualExcitationEnergy + TResidualExcitationEnergy ) { 
        Stopping = true;
        if ( PResidualExcitationEnergy <= 0.0 ) {
          TResidualExcitationEnergy = SqrtS - SumMasses;
        } else if ( TResidualExcitationEnergy <= 0.0 ) {
          PResidualExcitationEnergy = SqrtS - SumMasses;
        } else {
          G4double Fraction = (SqrtS - SumMasses) / 
                              (PResidualExcitationEnergy + TResidualExcitationEnergy);
          PResidualExcitationEnergy *= Fraction;
          TResidualExcitationEnergy *= Fraction;
        }
      }
    }

    if ( Stopping ) {
      // All 3-momenta of particles = 0
      // New projectile
      Ptmp.setPx( 0.0 ); Ptmp.setPy( 0.0 ); Ptmp.setPz( 0.0 );
      Ptmp.setE( PNucleonMass );
      Pprojectile = Ptmp; Pprojectile.transform( toLab );
      SelectedAntiBaryon->Set4Momentum( Pprojectile );

      // New projectile residual
      ProjectileResidualMassNumber       = PResidualMassNumber;
      ProjectileResidualCharge           = PResidualCharge;
      ProjectileResidualExcitationEnergy = PResidualExcitationEnergy;

      Ptmp.setE( PResidualMass + ProjectileResidualExcitationEnergy ); 
      Ptmp.transform( toLab );          
      ProjectileResidual4Momentum = Ptmp;

      // New target nucleon
      Ptmp.setE( TNucleonMass );
      Ptarget = Ptmp; Ptarget.transform( toLab );
      SelectedTargetNucleon->Set4Momentum( Ptarget );

      // New target residual
      TargetResidualMassNumber       = TResidualMassNumber;
      TargetResidualCharge           = TResidualCharge;
      TargetResidualExcitationEnergy = TResidualExcitationEnergy;

      Ptmp.setE( TResidualMass + TargetResidualExcitationEnergy ); 
      Ptmp.transform( toLab );          
      TargetResidual4Momentum = Ptmp;

      return true;
    }

    G4LorentzVector PResidual4Momentum = toCms * ProjectileResidual4Momentum;
    G4double YprojectileNucleus = PResidual4Momentum.rapidity();

    #ifdef debugAdjust
    G4cout << "YprojectileNucleus XcenterP " << YprojectileNucleus << G4endl;
    #endif

    G4LorentzVector TResidual4Momentum = toCms*TargetResidual4Momentum;
    G4double YtargetNucleus = TResidual4Momentum.rapidity();

    PResidualMass += PResidualExcitationEnergy;
    TResidualMass += TResidualExcitationEnergy;

    G4double M2projectile( 0.0 );
    G4double M2target( 0.0 );
    G4double WminusTarget( 0.0 );
    G4double WplusProjectile( 0.0 );

    G4ThreeVector PtNucleonP( 0.0, 0.0, 0.0 );
    G4double XplusNucleon( 0.0 );
    G4ThreeVector PtResidualP( 0.0, 0.0, 0.0 );
    G4double XplusResidual( 0.0 );

    G4ThreeVector PtNucleonT( 0.0, 0.0, 0.0 );
    G4double XminusNucleon( 0.0 );
    G4ThreeVector PtResidualT( 0.0, 0.0, 0.0 );
    G4double XminusResidual( 0.0 );

    G4int NumberOfTries( 0 );
    G4double ScaleFactor( 1.0 );
    G4bool OuterSuccess( true );

    do { // while ( ! OuterSuccess )

      OuterSuccess = true;

      do { // while ( SqrtS < std::sqrt( M2projectile ) + std::sqrt( M2target ) )

        NumberOfTries++;

        if ( NumberOfTries == 100*(NumberOfTries/100) ) { 
          // At large number of tries it would be better to reduce the values
          ScaleFactor /= 2.0;
          DcorP      *= ScaleFactor;
          DcorT      *= ScaleFactor;
          AveragePt2 *= ScaleFactor;
        }

        #ifdef debugAdjust
        //G4cout << "NumberOfTries ScaleFactor " << NumberOfTries << " " << ScaleFactor << G4endl;
        #endif

        if ( ProjectileResidualMassNumber > 1 ) {            
          PtNucleonP = GaussianPt( AveragePt2, maxPtSquare );
        } else {
          PtNucleonP = G4ThreeVector( 0.0, 0.0, 0.0 );
        }
        PtResidualP = -PtNucleonP;

        if ( TargetResidualMassNumber > 1 ) { 
          PtNucleonT = GaussianPt( AveragePt2, maxPtSquare );
        } else {
          PtNucleonT = G4ThreeVector( 0.0, 0.0, 0.0 );
        }
        PtResidualT = -PtNucleonT;

        G4double Mprojectile = std::sqrt( sqr( PNucleonMass )  + PtNucleonP.mag2() ) + 
                               std::sqrt( sqr( PResidualMass ) + PtResidualP.mag2() );
        M2projectile = sqr( Mprojectile );  // Uzhi 31.08.13

        G4double Mtarget = std::sqrt( sqr( TNucleonMass )  + PtNucleonT.mag2() ) + 
                           std::sqrt( sqr( TResidualMass ) + PtResidualT.mag2() );
        M2target = sqr( Mtarget );  // Uzhi 31.08.13        

        if ( SqrtS < Mprojectile + Mtarget ) {
          OuterSuccess = false; 
          continue;
        }

        G4bool InerSuccess = true;

        if ( ProjectileResidualMassNumber > 1 ) { 
          do {
            InerSuccess = true;
            G4ThreeVector tmpX = GaussianPt( DcorP*DcorP, 1.0 );
            G4double XcenterP = std::sqrt( sqr( PNucleonMass ) + PtNucleonP.mag2() ) / Mprojectile;
            XplusNucleon = XcenterP + tmpX.x();

            #ifdef debugAdjust
            //G4cout << "XplusNucleon 1 " << XplusNucleon << G4endl;
            //{ G4int Uzhi; G4cin >> Uzhi; }
            #endif

            if ( XplusNucleon <= 0.0 || XplusNucleon >= 1.0 ) {
              InerSuccess = false; 
              continue;
            }
            XplusResidual = 1.0 - XplusNucleon;
          } while ( ! InerSuccess );

          #ifdef debugAdjust
          //G4cout << "XplusNucleon XplusResidual 2 " << XplusNucleon 
          //       << " " << XplusResidual << G4endl;
          //{ G4int Uzhi; G4cin >> Uzhi; }
          #endif

        } else {
          XplusNucleon  = 1.0;
          XplusResidual = 1.0; // It must be 0
        }

        if ( TargetResidualMassNumber > 1 ) { 
          do {
            InerSuccess = true;

            G4ThreeVector tmpX = GaussianPt( DcorT*DcorT, 1.0 );
            G4double XcenterT = std::sqrt( sqr( TNucleonMass )  + PtNucleonT.mag2() ) / Mtarget;
            XminusNucleon = XcenterT + tmpX.x();
            if ( XminusNucleon <= 0.0 || XminusNucleon >= 1.0 ) {
              InerSuccess = false; 
              continue;
            }
            XminusResidual = 1.0 - XminusNucleon;
          } while ( ! InerSuccess );
        } else {
          XminusNucleon  = 1.0;
          XminusResidual = 1.0;  // It must be 0
        }

        #ifdef debugAdjust
        G4cout << "PtNucleonP " << PtNucleonP << " " << PtResidualP << G4endl
               << "XplusNucleon XplusResidual " << XplusNucleon << " " << XplusResidual << G4endl
               << "PtNucleonT " << PtNucleonT << " " << PtResidualT << G4endl
               << "XminusNucleon XminusResidual " << XminusNucleon << " " << XminusResidual 
               << G4endl;
        #endif

        M2projectile = ( sqr( PNucleonMass ) + PtNucleonP.mag2() )  / XplusNucleon + 
                       ( sqr( PResidualMass) + PtResidualP.mag2() ) / XplusResidual;
        M2target = ( sqr( TNucleonMass )  + PtNucleonT.mag2() )  / XminusNucleon + 
                   ( sqr( TResidualMass ) + PtResidualT.mag2() ) / XminusResidual;

      } while ( SqrtS < std::sqrt( M2projectile ) + std::sqrt( M2target ) );

      G4double DecayMomentum2 = sqr( S ) + sqr( M2projectile ) + sqr( M2target )
                                - 2.0*S*M2projectile - 2.0*S*M2target - 2.0*M2projectile*M2target;

      WplusProjectile = ( S + M2projectile - M2target + std::sqrt( DecayMomentum2 ) )/2.0/SqrtS;
      WminusTarget = SqrtS - M2projectile/WplusProjectile;

      G4double Mt2 = sqr( PNucleonMass ) + PtNucleonP.mag2();
      G4double Pz = WplusProjectile*XplusNucleon/2.0 - Mt2/(2.0*WplusProjectile*XplusNucleon);
      G4double E =  WplusProjectile*XplusNucleon/2.0 + Mt2/(2.0*WplusProjectile*XplusNucleon);
      G4double YprojectileNucleon = 0.5 * std::log( (E + Pz)/(E - Pz) );

      Mt2 = sqr( TNucleonMass ) + PtNucleonT.mag2();
      Pz = -WminusTarget*XminusNucleon/2.0 + Mt2/(2.0*WminusTarget*XminusNucleon);
      E =   WminusTarget*XminusNucleon/2.0 + Mt2/(2.0*WminusTarget*XminusNucleon);
      G4double YtargetNucleon = 0.5 * std::log( (E + Pz)/(E - Pz) ); 

      if ( std::abs( YtargetNucleon - YtargetNucleus ) > 2         || 
           std::abs( YprojectileNucleon - YprojectileNucleus ) > 2 ||
           YprojectileNucleon < YtargetNucleon ) {        
        OuterSuccess = false;
        continue;
      } 

    } while ( ! OuterSuccess );

    #ifdef debugAdjust
    G4cout << "PtNucleonP " << PtNucleonP << G4endl;
    #endif

    G4double Mt2 = sqr( PNucleonMass ) + PtNucleonP.mag2();
    G4double Pz = WplusProjectile*XplusNucleon/2.0 - Mt2/(2.0*WplusProjectile*XplusNucleon);
    G4double E =  WplusProjectile*XplusNucleon/2.0 + Mt2/(2.0*WplusProjectile*XplusNucleon);

    Pprojectile.setPx( PtNucleonP.x() ); Pprojectile.setPy( PtNucleonP.y() );
    Pprojectile.setPz( Pz );             Pprojectile.setE( E ); 
    Pprojectile.transform( toLab );
    SelectedAntiBaryon->Set4Momentum( Pprojectile );

    // New projectile residual
    ProjectileResidualMassNumber       = PResidualMassNumber;
    ProjectileResidualCharge           = PResidualCharge;
    ProjectileResidualExcitationEnergy = PResidualExcitationEnergy;

    #ifdef debugAdjust
    G4cout << "PResidualMass PtResidualP " << PResidualMass << " " << PtResidualP << G4endl;
    #endif

    if ( ProjectileResidualMassNumber != 0 ) {
      Mt2 = sqr( PResidualMass ) + PtResidualP.mag2();
      Pz = WplusProjectile*XplusResidual/2.0 - Mt2/(2.0*WplusProjectile*XplusResidual);
      E  = WplusProjectile*XplusResidual/2.0 + Mt2/(2.0*WplusProjectile*XplusResidual);
      ProjectileResidual4Momentum.setPx( PtResidualP.x() );
      ProjectileResidual4Momentum.setPy( PtResidualP.y() );
      ProjectileResidual4Momentum.setPz( Pz );
      ProjectileResidual4Momentum.setE( E );
      ProjectileResidual4Momentum.transform( toLab );
    } else {
      ProjectileResidual4Momentum = G4LorentzVector( 0.0, 0.0, 0.0, 0.0 );
    } 

    #ifdef debugAdjust
    G4cout << "Pr N R " << Pprojectile << G4endl << "       " 
           << ProjectileResidual4Momentum << G4endl;
    #endif

    Mt2 = sqr( TNucleonMass ) + PtNucleonT.mag2();
    Pz = -WminusTarget*XminusNucleon/2.0 + Mt2/(2.0*WminusTarget*XminusNucleon);
    E =   WminusTarget*XminusNucleon/2.0 + Mt2/(2.0*WminusTarget*XminusNucleon);

    Ptarget.setPx( PtNucleonT.x() ); Ptarget.setPy( PtNucleonT.y() );
    Ptarget.setPz( Pz );             Ptarget.setE( E ); 
    Ptarget.transform( toLab );
    SelectedTargetNucleon->Set4Momentum( Ptarget );

    // New target residual
    TargetResidualMassNumber       = TResidualMassNumber;
    TargetResidualCharge           = TResidualCharge;
    TargetResidualExcitationEnergy = TResidualExcitationEnergy;

    if ( TargetResidualMassNumber != 0 ) {
      Mt2 = sqr( TResidualMass ) + PtResidualT.mag2();
      Pz = -WminusTarget*XminusResidual/2.0 + Mt2/(2.0*WminusTarget*XminusResidual);
      E =   WminusTarget*XminusResidual/2.0 + Mt2/(2.0*WminusTarget*XminusResidual);

      TargetResidual4Momentum.setPx( PtResidualT.x() );
      TargetResidual4Momentum.setPy( PtResidualT.y() );
      TargetResidual4Momentum.setPz( Pz );
      TargetResidual4Momentum.setE( E) ;
      TargetResidual4Momentum.transform( toLab );
    } else {
      TargetResidual4Momentum = G4LorentzVector( 0.0, 0.0, 0.0, 0.0 );
    } 

    #ifdef debugAdjust
    G4cout << "Tr N R " << Ptarget << G4endl << "       " << TargetResidual4Momentum << G4endl;
    #endif

    return true;

  }

}


//============================================================================

G4ExcitedStringVector* G4FTFModel::BuildStrings() {
  // Loop over all collisions; find all primaries, and all targets 
  // (targets may be duplicate in the List (to unique G4VSplitableHadrons) ).

  G4ExcitedStringVector* strings = new G4ExcitedStringVector();
  G4ExcitedString* FirstString( 0 );     // If there will be a kink,
  G4ExcitedString* SecondString( 0 );    // two strings will be produced.

  if ( ! GetProjectileNucleus() ) {

    std::vector< G4VSplitableHadron* > primaries;
    theParticipants.StartLoop();
    while ( theParticipants.Next() ) {
      const G4InteractionContent& interaction = theParticipants.GetInteraction();
      //  do not allow for duplicates ...
      if ( interaction.GetStatus() ) {
        if ( primaries.end() == std::find( primaries.begin(), primaries.end(), 
                                           interaction.GetProjectile() ) ) {
          primaries.push_back( interaction.GetProjectile() );
        }
      }
    }

    #ifdef debugBuildString
    G4cout << "G4FTFModel::BuildStrings()" << G4endl
           << "Number of projectile strings " << primaries.size() << G4endl;
    #endif

    for ( unsigned int ahadron = 0; ahadron < primaries.size(); ahadron++ ) {
      G4bool isProjectile( true );
      //G4cout << "primaries[ahadron] " << primaries[ahadron] << G4endl;
      //if ( primaries[ahadron]->GetStatus() <= 1 ) isProjectile=true;
      FirstString = 0; SecondString = 0;
      theExcitation->CreateStrings( primaries[ ahadron ], isProjectile, 
                                    FirstString, SecondString, theParameters );
      if ( FirstString  != 0 ) strings->push_back( FirstString );
      if ( SecondString != 0 ) strings->push_back( SecondString );

      #ifdef debugBuildString
      G4cout << "FirstString & SecondString? " << FirstString << " " << SecondString << G4endl
             << "Quarks on the FirstString ends " << FirstString->GetRightParton()->GetPDGcode()
             << " " << FirstString->GetLeftParton()->GetPDGcode() << G4endl;
      #endif

    }

    #ifdef debugBuildString
    G4cout << "Check 1 string " << strings->operator[](0)->GetRightParton()->GetPDGcode() 
           << " " << strings->operator[](0)->GetLeftParton()->GetPDGcode() << G4endl << G4endl;
    #endif

    std::for_each( primaries.begin(), primaries.end(), DeleteVSplitableHadron() );
    primaries.clear();

  } else {  // Projectile is a nucleus

    #ifdef debugBuildString
    G4cout << "Building of projectile-like strings" << G4endl;
    #endif

    G4bool isProjectile = true;
    for ( G4int ahadron = 0; ahadron < NumberOfInvolvedNucleonsOfProjectile; ahadron++ ) {

      #ifdef debugBuildString
      G4cout << "Nucleon #, status, intCount " << ahadron << " "
             << TheInvolvedNucleonsOfProjectile[ ahadron ]->GetSplitableHadron()->GetStatus() 
             << " " << TheInvolvedNucleonsOfProjectile[ ahadron ]->GetSplitableHadron()
                       ->GetSoftCollisionCount();
      #endif

      G4VSplitableHadron* aProjectile = 
          TheInvolvedNucleonsOfProjectile[ ahadron ]->GetSplitableHadron();

      #ifdef debugBuildString
      G4cout << G4endl << "ahadron aProjectile Status " << ahadron << " " << aProjectile
             << " " << aProjectile->GetStatus() << G4endl;
      #endif

      if ( aProjectile->GetStatus() == 0 ) { // A nucleon took part in non-diffractive interaction

        #ifdef debugBuildString
        G4cout << "Case1 aProjectile->GetStatus() == 0 " << G4endl;
        #endif

        FirstString = 0; SecondString = 0;
        theExcitation->CreateStrings( 
                           TheInvolvedNucleonsOfProjectile[ ahadron ]->GetSplitableHadron(),
                           isProjectile, FirstString, SecondString, theParameters );
        if ( FirstString  != 0 ) strings->push_back( FirstString );
        if ( SecondString != 0 ) strings->push_back( SecondString );
      } else if ( aProjectile->GetStatus() == 1 && aProjectile->GetSoftCollisionCount() != 0 ) { 
        // Nucleon took part in diffractive interaction

        #ifdef debugBuildString
        G4cout << "Case2 aProjectile->GetStatus() !=0 St==1 SoftCol!=0" << G4endl;
        #endif

        FirstString = 0; SecondString = 0;
        theExcitation->CreateStrings( 
                           TheInvolvedNucleonsOfProjectile[ ahadron ]->GetSplitableHadron(),
                           isProjectile, FirstString, SecondString, theParameters );
        if ( FirstString  != 0 ) strings->push_back( FirstString );
        if ( SecondString != 0 ) strings->push_back( SecondString );
      } else if ( aProjectile->GetStatus() == 1  &&  aProjectile->GetSoftCollisionCount() == 0  &&
                  HighEnergyInter ) {
        // Nucleon was considered as a paricipant of an interaction,
        // but the interaction was skipped due to annihilation.
        // It is now considered as an involved nucleon at high energies.

        #ifdef debugBuildString
        G4cout << "Case3 aProjectile->GetStatus() !=0 St==1 SoftCol==0" << G4endl;
        #endif

        FirstString = 0; SecondString = 0;
        theExcitation->CreateStrings( 
                           TheInvolvedNucleonsOfProjectile[ ahadron ]->GetSplitableHadron(),
                           isProjectile, FirstString, SecondString, theParameters );
        if ( FirstString  != 0 ) strings->push_back( FirstString );
        if ( SecondString != 0 ) strings->push_back( SecondString );

        #ifdef debugBuildString
        G4cout << " Strings are built for nucleon marked for an interaction, but"
               << " the interaction was skipped." << G4endl;
        #endif

      } else if ( aProjectile->GetStatus() == 2 ) {
        // Nucleon which was involved in the Reggeon cascading

        #ifdef debugBuildString
        G4cout << "Case4 aProjectile->GetStatus() !=0 St==2 " << G4endl;
        #endif

        FirstString = 0; SecondString = 0;
        theExcitation->CreateStrings( 
                         TheInvolvedNucleonsOfProjectile[ ahadron ]->GetSplitableHadron(),
                         isProjectile, FirstString, SecondString, theParameters );
        if ( FirstString  != 0 ) strings->push_back( FirstString );
        if ( SecondString != 0 ) strings->push_back( SecondString );

        #ifdef debugBuildString
        G4cout << " Strings are build for involved nucleon." << G4endl;
        #endif

      } else {

        #ifdef debugBuildString
        G4cout << "Case5 " << G4endl;
        #endif

        //TheInvolvedNucleonsOfProjectile[ ahadron ]->Hit( 0 );
        //G4cout << TheInvolvedNucleonsOfProjectile[ ahadron ]->GetSplitableHadron() << G4endl;

        #ifdef debugBuildString
        G4cout << " No string" << G4endl;
        #endif

      }
    }
  } 

  #ifdef debugBuildString
  G4cout << "Building of target-like strings" << G4endl;
  #endif

  G4bool isProjectile = false;
  for ( G4int ahadron = 0; ahadron < NumberOfInvolvedNucleonsOfTarget; ahadron++ ) {
    G4VSplitableHadron* aNucleon = TheInvolvedNucleonsOfTarget[ ahadron ]->GetSplitableHadron();

    #ifdef debugBuildString
    G4cout << "Nucleon #, status, intCount " << aNucleon << " " << ahadron << " "
           << aNucleon->GetStatus() << " " << aNucleon->GetSoftCollisionCount();
    #endif

    if ( aNucleon->GetStatus() == 0 ) { // A nucleon took part in non-diffractive interaction
      FirstString = 0 ; SecondString = 0;
      theExcitation->CreateStrings( aNucleon, isProjectile, 
                                    FirstString, SecondString, theParameters );
      if ( FirstString  != 0 ) strings->push_back( FirstString );
      if ( SecondString != 0 ) strings->push_back( SecondString );

      #ifdef debugBuildString
      G4cout << " 1 case A string is build" << G4endl;
      #endif

    } else if ( aNucleon->GetStatus() == 1  &&  aNucleon->GetSoftCollisionCount() != 0 ) {
      // A nucleon took part in diffractive interaction
      FirstString = 0; SecondString = 0;
      theExcitation->CreateStrings( aNucleon, isProjectile,
                                    FirstString, SecondString, theParameters );
      if ( FirstString  != 0 ) strings->push_back( FirstString );
      if ( SecondString != 0 ) strings->push_back( SecondString );

      #ifdef debugBuildString
      G4cout << "2 case A string is build, nucleon was excited." << G4endl;
      #endif

    } else if ( aNucleon->GetStatus() == 1  &&  aNucleon->GetSoftCollisionCount() == 0  &&
                HighEnergyInter ) {
      // A nucleon was considered as a participant but due to annihilation
      // its interactions were skipped. It will be considered as involved one
      // at high energies.
      FirstString = 0; SecondString = 0;
      theExcitation->CreateStrings( aNucleon, isProjectile,
                                    FirstString, SecondString, theParameters );
      if ( FirstString  != 0 ) strings->push_back( FirstString );
      if ( SecondString != 0 ) strings->push_back( SecondString );

      #ifdef debugBuildString
      G4cout << "3 case A string is build" << G4endl;
      #endif

    } else if ( aNucleon->GetStatus() == 1  &&  aNucleon->GetSoftCollisionCount() == 0  &&
                ! HighEnergyInter ) {
      // A nucleon was considered as a participant but due to annihilation
      // its interactions were skipped. It will be returned to nucleus
      // at low energies energies.
      aNucleon->SetStatus( 4 );
      // ????????? delete aNucleon;

      #ifdef debugBuildString
      G4cout << "4 case A string is not build" << G4endl;
      #endif

    } else if ( aNucleon->GetStatus() == 2 ) { // A nucleon was involved in Reggeon cascading
      FirstString = 0; SecondString = 0;
      theExcitation->CreateStrings( aNucleon, isProjectile, 
                                    FirstString, SecondString, theParameters );
      if ( FirstString  != 0 ) strings->push_back( FirstString );
      if ( SecondString != 0 ) strings->push_back( SecondString );

      #ifdef debugBuildString
      G4cout << "5 case A string is build" << G4endl;
      #endif

    } else {

      #ifdef debugBuildString
      G4cout << "6 case No string" << G4endl;
      #endif

    }
  }

  #ifdef debugBuildString
  G4cout << G4endl << "theAdditionalString.size() " << theAdditionalString.size() 
         << G4endl << G4endl;
  #endif

  isProjectile = true;
  if ( theAdditionalString.size() != 0 ) {
    for ( unsigned int  ahadron = 0; ahadron < theAdditionalString.size(); ahadron++ ) {
      // if ( theAdditionalString[ ahadron ]->GetStatus() <= 1 ) isProjectile = true; 
      FirstString = 0; SecondString = 0;
      theExcitation->CreateStrings( theAdditionalString[ ahadron ], isProjectile,
                                    FirstString, SecondString, theParameters );
      if ( FirstString  != 0 ) strings->push_back( FirstString );
      if ( SecondString != 0 ) strings->push_back( SecondString );
    }
  }

  //for ( unsigned int ahadron = 0; ahadron < strings->size(); ahadron++ ) {
  //  G4cout << ahadron << " " << strings->operator[]( ahadron )->GetRightParton()->GetPDGcode()
  //         << " " << strings->operator[]( ahadron )->GetLeftParton()->GetPDGcode() << G4endl;
  //}
  //G4cout << "------------------------" << G4endl;

  return strings;
}


//============================================================================

void G4FTFModel::GetResiduals() {
  // This method is needed for the correct application of G4PrecompoundModelInterface

  #ifdef debugFTFmodel
  G4cout << "GetResiduals(): HighEnergyInter? GetProjectileNucleus()?"
         << HighEnergyInter << " " << GetProjectileNucleus() << G4endl;
  #endif

  if ( HighEnergyInter ) {

    #ifdef debugFTFmodel
    G4cout << "NumberOfInvolvedNucleonsOfTarget "<< NumberOfInvolvedNucleonsOfTarget << G4endl;
    #endif

    G4double DeltaExcitationE = TargetResidualExcitationEnergy / 
                                G4double( NumberOfInvolvedNucleonsOfTarget );
    G4LorentzVector DeltaPResidualNucleus = TargetResidual4Momentum /
                                            G4double( NumberOfInvolvedNucleonsOfTarget );

    for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; i++ ) {
      G4Nucleon* aNucleon = TheInvolvedNucleonsOfTarget[i];

      #ifdef debugFTFmodel
      G4VSplitableHadron* targetSplitable = aNucleon->GetSplitableHadron();
      G4cout << i << " Hit? " << aNucleon->AreYouHit() << " " << targetSplitable << G4endl;
      if ( targetSplitable ) G4cout << i << "Status " << targetSplitable->GetStatus() << G4endl;
      #endif

      G4LorentzVector tmp = -DeltaPResidualNucleus;
      aNucleon->SetMomentum( tmp );
      aNucleon->SetBindingEnergy( DeltaExcitationE );
    }

    if ( ! GetProjectileNucleus() ) return; // The projectile is a hadron

    #ifdef debugFTFmodel
    G4cout << "NumberOfInvolvedNucleonsOfProjectile " << NumberOfInvolvedNucleonsOfProjectile
           << G4endl << "ProjectileResidualExcitationEnergy ProjectileResidual4Momentum "
           << ProjectileResidualExcitationEnergy << "  " << ProjectileResidual4Momentum << G4endl;
    #endif

    DeltaExcitationE = ProjectileResidualExcitationEnergy /
                       G4double( NumberOfInvolvedNucleonsOfProjectile );
    DeltaPResidualNucleus = ProjectileResidual4Momentum /
                            G4double( NumberOfInvolvedNucleonsOfProjectile );

    for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfProjectile; i++ ) {
      G4Nucleon* aNucleon = TheInvolvedNucleonsOfProjectile[i];

      #ifdef debugFTFmodel
      G4VSplitableHadron* projSplitable = aNucleon->GetSplitableHadron();
      G4cout << i << " Hit? " << aNucleon->AreYouHit() << " " << projSplitable << G4endl;
      if ( projSplitable ) G4cout << i << "Status " << projSplitable->GetStatus() << G4endl;
      #endif

      G4LorentzVector tmp = -DeltaPResidualNucleus;
      aNucleon->SetMomentum( tmp );
      aNucleon->SetBindingEnergy( DeltaExcitationE );
    }
  
    #ifdef debugFTFmodel
    G4cout << "End projectile" << G4endl;
    #endif
   
  } else {

    #ifdef debugFTFmodel
    G4cout << "Low energy interaction: Target nucleus --------------" << G4endl
           << "Tr ResidualMassNumber Tr ResidualCharge Tr ResidualExcitationEnergy  "
           << TargetResidualMassNumber << " " << TargetResidualCharge << " "
           << TargetResidualExcitationEnergy << G4endl;
    #endif

    G4int NumberOfTargetParticipant( 0 );
    for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; i++ ) {
      G4Nucleon* aNucleon = TheInvolvedNucleonsOfTarget[i];
      G4VSplitableHadron* targetSplitable = aNucleon->GetSplitableHadron();
      if ( targetSplitable->GetSoftCollisionCount() != 0 ) NumberOfTargetParticipant++;
    } 

    G4double DeltaExcitationE( 0.0 );
    G4LorentzVector DeltaPResidualNucleus = G4LorentzVector( 0.0, 0.0, 0.0, 0.0 );

    if ( NumberOfTargetParticipant != 0 ) {
      DeltaExcitationE = TargetResidualExcitationEnergy / G4double( NumberOfTargetParticipant );
      DeltaPResidualNucleus = TargetResidual4Momentum / G4double( NumberOfTargetParticipant );
    }

    for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; i++ ) {
      G4Nucleon* aNucleon = TheInvolvedNucleonsOfTarget[i];
      G4VSplitableHadron* targetSplitable = aNucleon->GetSplitableHadron();
      if ( targetSplitable->GetSoftCollisionCount() != 0 ) {
        G4LorentzVector tmp = -DeltaPResidualNucleus;
        aNucleon->SetMomentum( tmp );
        aNucleon->SetBindingEnergy( DeltaExcitationE );
      } else {
        delete targetSplitable;  
        targetSplitable = 0; 
        aNucleon->Hit( targetSplitable );
        aNucleon->SetBindingEnergy( 0.0 );
      }
    } 

    #ifdef debugFTFmodel
    G4cout << "NumberOfTargetParticipant " << NumberOfTargetParticipant << G4endl
           << "TargetResidual4Momentum  " << TargetResidual4Momentum << G4endl;
    #endif

    if ( ! GetProjectileNucleus() ) return; // The projectile is a hadron

    #ifdef debugFTFmodel
    G4cout << "Low energy interaction: Projectile nucleus --------------" << G4endl
           << "Pr ResidualMassNumber Pr ResidualCharge Pr ResidualExcitationEnergy "
           << ProjectileResidualMassNumber << " " << ProjectileResidualCharge << " "
           << ProjectileResidualExcitationEnergy << G4endl;
    #endif

    G4int NumberOfProjectileParticipant( 0 );
    for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfProjectile; i++ ) {
      G4Nucleon* aNucleon = TheInvolvedNucleonsOfProjectile[i];
      G4VSplitableHadron* projectileSplitable = aNucleon->GetSplitableHadron();
      if ( projectileSplitable->GetSoftCollisionCount() != 0 ) 
        NumberOfProjectileParticipant++;
      }
 
      #ifdef debugFTFmodel
      G4cout << "NumberOfProjectileParticipant" << G4endl;
      #endif

      DeltaExcitationE = 0.0;
      DeltaPResidualNucleus = G4LorentzVector( 0.0, 0.0, 0.0, 0.0 );

      if ( NumberOfProjectileParticipant != 0 ) {
        DeltaExcitationE = ProjectileResidualExcitationEnergy / 
                           G4double( NumberOfProjectileParticipant );
        DeltaPResidualNucleus = ProjectileResidual4Momentum /
                                G4double( NumberOfProjectileParticipant );
      }
      //G4cout << "DeltaExcitationE DeltaPResidualNucleus " << DeltaExcitationE
      //       << " " << DeltaPResidualNucleus << G4endl;
      for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfProjectile; i++ ) {
        G4Nucleon* aNucleon = TheInvolvedNucleonsOfProjectile[i];
        G4VSplitableHadron* projectileSplitable = aNucleon->GetSplitableHadron();
        if ( projectileSplitable->GetSoftCollisionCount() != 0 ) {
          G4LorentzVector tmp = -DeltaPResidualNucleus;
          aNucleon->SetMomentum( tmp );
          aNucleon->SetBindingEnergy( DeltaExcitationE );
        } else {
          delete projectileSplitable;  
          projectileSplitable = 0; 
          aNucleon->Hit( projectileSplitable );
          aNucleon->SetBindingEnergy( 0.0 );
        } 
      } 

      #ifdef debugFTFmodel
      G4cout << "NumberOfProjectileParticipant " << NumberOfProjectileParticipant << G4endl
             << "ProjectileResidual4Momentum " << ProjectileResidual4Momentum << G4endl;
      #endif

    }

    #ifdef debugFTFmodel
    G4cout << "End GetResiduals -----------------" << G4endl;
    #endif

}


//============================================================================

G4ThreeVector G4FTFModel::GaussianPt( G4double AveragePt2, G4double maxPtSquare ) const {
  //  @@ this method is used in FTFModel as well. Should go somewhere common!

  G4double Pt2( 0.0 );
  if ( AveragePt2 <= 0.0 ) {
    Pt2 = 0.0;
  } else {
    Pt2 = -AveragePt2 * std::log( 1.0 + G4UniformRand() * 
                                        ( std::exp( -maxPtSquare/AveragePt2 ) -1.0 ) );
  }
  G4double Pt = std::sqrt( Pt2 );
  G4double phi = G4UniformRand() * twopi;

  return G4ThreeVector( Pt*std::cos(phi), Pt*std::sin(phi), 0.0 );    
}


//============================================================================

void G4FTFModel::ModelDescription( std::ostream& desc ) const {
  desc << "please add description here" << G4endl;
}
