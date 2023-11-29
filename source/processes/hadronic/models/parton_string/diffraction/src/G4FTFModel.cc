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
#include "G4KineticTrack.hh"
#include "G4HyperNucleiProperties.hh"
#include "G4Exp.hh"
#include "G4Log.hh"

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
  // ---> JVY theParameters = 0;
  theParameters = new G4FTFParameters();
  //
  NumberOfInvolvedNucleonsOfTarget = 0;
  NumberOfInvolvedNucleonsOfProjectile= 0;
  for ( G4int i = 0; i < 250; ++i ) {
    TheInvolvedNucleonsOfTarget[i] = 0;
    TheInvolvedNucleonsOfProjectile[i] = 0;
  }

  //LowEnergyLimit = 2000.0*MeV;
  LowEnergyLimit = 1000.0*MeV;

  HighEnergyInter = true;

  G4LorentzVector tmp( 0.0, 0.0, 0.0, 0.0 );
  ProjectileResidual4Momentum        = tmp;
  ProjectileResidualMassNumber       = 0;
  ProjectileResidualCharge           = 0;
  ProjectileResidualLambdaNumber     = 0;
  ProjectileResidualExcitationEnergy = 0.0;

  TargetResidual4Momentum            = tmp;
  TargetResidualMassNumber           = 0;
  TargetResidualCharge               = 0;
  TargetResidualExcitationEnergy     = 0.0;

  Bimpact = -1.0;
  BinInterval = false;
  Bmin = 0.0; 
  Bmax = 0.0;
  NumberOfProjectileSpectatorNucleons = 0;
  NumberOfTargetSpectatorNucleons = 0;
  NumberOfNNcollisions = 0;

  SetEnergyMomentumCheckLevels( 2.0*perCent, 150.0*MeV );
}


//============================================================================

struct DeleteVSplitableHadron { void operator()( G4VSplitableHadron* aH ) { delete aH; } };


//============================================================================

G4FTFModel::~G4FTFModel() {
   // Because FTF model can be called for various particles
   //
   // ---> NOTE (JVY): This statement below is no longer true !!! 
   // theParameters must be erased at the end of each call.
   // Thus the delete is also in G4FTFModel::GetStrings() method.
   // ---> JVY
   //
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
     for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; ++i ) {
       G4VSplitableHadron* aNucleon = TheInvolvedNucleonsOfTarget[i]->GetSplitableHadron();
       if ( aNucleon ) delete aNucleon;
     }
   }

   // Erasing of projectile involved nucleons.
   if ( NumberOfInvolvedNucleonsOfProjectile != 0 ) {
     for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfProjectile; ++i ) {
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

  theParticipants.Clean();

  theParticipants.SetProjectileNucleus( 0 );

  G4LorentzVector tmp( 0.0, 0.0, 0.0, 0.0 );
  ProjectileResidualMassNumber       = 0;
  ProjectileResidualCharge           = 0;
  ProjectileResidualLambdaNumber     = 0;
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
    ProjectileResidualMassNumber = std::abs( theProjectile.GetDefinition()->GetBaryonNumber() );
    ProjectileResidualCharge = G4int( theProjectile.GetDefinition()->GetPDGCharge() );
    PlabPerParticle = theProjectile.GetMomentum().z();
    ProjectileResidualExcitationEnergy = 0.0;
    //G4double ProjectileResidualMass = theProjectile.GetMass();
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
      ProjectileResidualMassNumber = theProjectile.GetDefinition()->GetBaryonNumber();
      ProjectileResidualCharge = G4int( theProjectile.GetDefinition()->GetPDGCharge() );
      ProjectileResidualLambdaNumber = theProjectile.GetDefinition()->GetNumberOfLambdasInHypernucleus();
      PlabPerParticle = theProjectile.GetMomentum().z() / ProjectileResidualMassNumber;
      if ( PlabPerParticle < LowEnergyLimit ) {
        HighEnergyInter = false;
      } else {
        HighEnergyInter = true;
      }
      theParticipants.InitProjectileNucleus( ProjectileResidualMassNumber, ProjectileResidualCharge,
					     ProjectileResidualLambdaNumber );
    } else if ( theProjectile.GetDefinition()->GetBaryonNumber() < -1 ) { 
      // Projectile is an anti-nucleus
      ProjectileResidualMassNumber = std::abs( theProjectile.GetDefinition()->GetBaryonNumber() );
      ProjectileResidualCharge = std::abs( G4int( theProjectile.GetDefinition()->GetPDGCharge() ) );
      ProjectileResidualLambdaNumber = theProjectile.GetDefinition()->GetNumberOfAntiLambdasInAntiHypernucleus();
      PlabPerParticle = theProjectile.GetMomentum().z() / ProjectileResidualMassNumber;
      if ( PlabPerParticle < LowEnergyLimit ) {
        HighEnergyInter = false;
      } else {
        HighEnergyInter = true;
      }
      theParticipants.InitProjectileNucleus( ProjectileResidualMassNumber, ProjectileResidualCharge,
                                             ProjectileResidualLambdaNumber );
      theParticipants.GetProjectileNucleus()->StartLoop();
      G4Nucleon* aNucleon;
      while ( ( aNucleon = theParticipants.GetProjectileNucleus()->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
        if ( aNucleon->GetDefinition() == G4Proton::Definition() ) {
          aNucleon->SetParticleType( G4AntiProton::Definition() ); 
        } else if ( aNucleon->GetDefinition() == G4Neutron::Definition() ) {
          aNucleon->SetParticleType( G4AntiNeutron::Definition() );
        } else if ( aNucleon->GetDefinition() == G4Lambda::Definition() ) {
	  aNucleon->SetParticleType( G4AntiLambda::Definition() );
	}
      }
    }

    G4ThreeVector BoostVector = theProjectile.GetMomentum() / theProjectile.GetTotalEnergy();
    theParticipants.GetProjectileNucleus()->DoLorentzBoost( BoostVector );
    theParticipants.GetProjectileNucleus()->DoLorentzContraction( BoostVector );
    ProjectileResidualExcitationEnergy = 0.0;
    //G4double ProjectileResidualMass = theProjectile.GetMass();
    ProjectileResidual4Momentum.setVect( theProjectile.GetMomentum() );
    ProjectileResidual4Momentum.setE( theProjectile.GetTotalEnergy() );
  }

  // Init target nucleus (assumed to be never a hypernucleus)
  theParticipants.Init( aNucleus.GetA_asInt(), aNucleus.GetZ_asInt() );

  NumberOfProjectileSpectatorNucleons = std::abs( theProjectile.GetDefinition()->GetBaryonNumber() );
  NumberOfTargetSpectatorNucleons = aNucleus.GetA_asInt();
  NumberOfNNcollisions = 0;

  // reset/recalculate everything for the new interaction
  theParameters->InitForInteraction( theProjectile.GetDefinition(), aNucleus.GetA_asInt(),
                                     aNucleus.GetZ_asInt(), PlabPerParticle ); 

  if ( theAdditionalString.size() != 0 ) {
    std::for_each( theAdditionalString.begin(), theAdditionalString.end(), 
                   DeleteVSplitableHadron() );
  }
  theAdditionalString.clear();

  #ifdef debugFTFmodel
  G4cout << "FTF end of Init" << G4endl << G4endl;
  #endif

  // In the case of Hydrogen target, for non-ion hadron projectiles,
  // do NOT simulate quasi-elastic (by forcing to 0 the probability of
  // elastic scatering in theParameters - which is used only by FTF).
  // This is necessary because in this case quasi-elastic on a target nucleus
  // with only one nucleon would be identical to the hadron elastic scattering,
  // and the latter is already included in the elastic process 
  // (i.e. G4HadronElasticProcess).
  if ( std::abs( theProjectile.GetDefinition()->GetBaryonNumber() ) <= 1  &&
       aNucleus.GetA_asInt() < 2 ) theParameters->SetProbabilityOfElasticScatt( 0.0 );

  if ( SampleBinInterval() ) theParticipants.SetBminBmax( GetBmin(), GetBmax() );
}


//============================================================================

G4ExcitedStringVector* G4FTFModel::GetStrings() { 

  #ifdef debugFTFmodel
  G4cout << "G4FTFModel::GetStrings() " << G4endl;
  #endif

  G4ExcitedStringVector* theStrings = new G4ExcitedStringVector;
  theParticipants.GetList( theProjectile, theParameters );

  SetImpactParameter( theParticipants.GetImpactParameter() );

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

    BuildStrings( theStrings );

    #ifdef debugFTFmodel
    G4cout << "FTF BuildStrings " << theStrings << " OK" << G4endl
           << "FTF GetResiduals of Nuclei " << G4endl;
    #endif

    GetResiduals();

    /*
    if ( theParameters != 0 ) {
      delete theParameters;
      theParameters = 0;
    }
    */
  } else if ( ! GetProjectileNucleus() ) {
    // Erase the hadron projectile
    std::vector< G4VSplitableHadron* > primaries;
    theParticipants.StartLoop();
    while ( theParticipants.Next() ) {  /* Loop checking, 10.08.2015, A.Ribon */
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
  for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfProjectile; ++i ) {
    aNucleon = TheInvolvedNucleonsOfProjectile[i]->GetSplitableHadron();
    if ( aNucleon ) delete aNucleon;
  } 
  NumberOfInvolvedNucleonsOfProjectile = 0;

  // Erase the target nucleons
  for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; ++i ) {
    aNucleon = TheInvolvedNucleonsOfTarget[i]->GetSplitableHadron();
    if ( aNucleon ) delete aNucleon;
  } 
  NumberOfInvolvedNucleonsOfTarget = 0;

  #ifdef debugFTFmodel
  G4cout << "End of FTF. Go to fragmentation" << G4endl
         << "To continue - enter 1, to stop - ^C" << G4endl;
  #endif

  theParticipants.Clean();

  return theStrings;
}


//============================================================================

void G4FTFModel::StoreInvolvedNucleon() {
  //To store nucleons involved in the interaction

  NumberOfInvolvedNucleonsOfTarget = 0;

  G4V3DNucleus* theTargetNucleus = GetTargetNucleus();
  theTargetNucleus->StartLoop();

  G4Nucleon* aNucleon;
  while ( ( aNucleon = theTargetNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
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
  while ( ( aProjectileNucleon = theProjectileNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
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

  #ifdef debugReggeonCascade
  G4cout << "G4FTFModel::ReggeonCascade -----------" << G4endl
         << "theProjectile.GetTotalMomentum() " << theProjectile.GetTotalMomentum() << G4endl
         << "theProjectile.GetTotalEnergy() " << theProjectile.GetTotalEnergy() << G4endl
         << "ExcitationE/WN " << theParameters->GetExcitationEnergyPerWoundedNucleon() << G4endl;
  #endif

  G4int InitNINt = NumberOfInvolvedNucleonsOfTarget;

  // Reggeon cascading in target nucleus
  for ( G4int InvTN = 0; InvTN < InitNINt; InvTN++ ) { 
    G4Nucleon* aTargetNucleon = TheInvolvedNucleonsOfTarget[ InvTN ];

    G4double CreationTime = aTargetNucleon->GetSplitableHadron()->GetTimeOfCreation();

    G4double XofWoundedNucleon = aTargetNucleon->GetPosition().x();
    G4double YofWoundedNucleon = aTargetNucleon->GetPosition().y();
           
    G4V3DNucleus* theTargetNucleus = GetTargetNucleus();
    theTargetNucleus->StartLoop();

    G4Nucleon* Neighbour(0);
    while ( ( Neighbour = theTargetNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
      if ( ! Neighbour->AreYouHit() ) {
        G4double impact2 = sqr( XofWoundedNucleon - Neighbour->GetPosition().x() ) +
                           sqr( YofWoundedNucleon - Neighbour->GetPosition().y() );

        if ( G4UniformRand() < theParameters->GetCofNuclearDestruction() *
                               G4Exp( -impact2 / theParameters->GetR2ofNuclearDestruction() )
           ) {  
          // The neighbour nucleon is involved in the reggeon cascade
          TheInvolvedNucleonsOfTarget[ NumberOfInvolvedNucleonsOfTarget ] = Neighbour;
          NumberOfInvolvedNucleonsOfTarget++;

          G4VSplitableHadron* targetSplitable; 
          targetSplitable = new G4DiffractiveSplitableHadron( *Neighbour ); 

          Neighbour->Hit( targetSplitable );
          targetSplitable->SetTimeOfCreation( CreationTime ); 
          targetSplitable->SetStatus( 3 );  // 2->3
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
  G4int InitNINp = NumberOfInvolvedNucleonsOfProjectile;

  //for ( G4int InvPN = 0; InvPN < NumberOfInvolvedNucleonsOfProjectile; InvPN++ ) { 
  for ( G4int InvPN = 0; InvPN < InitNINp; InvPN++ ) { 
    G4Nucleon* aProjectileNucleon = TheInvolvedNucleonsOfProjectile[ InvPN ];

    G4double CreationTime = aProjectileNucleon->GetSplitableHadron()->GetTimeOfCreation();

    G4double XofWoundedNucleon = aProjectileNucleon->GetPosition().x();
    G4double YofWoundedNucleon = aProjectileNucleon->GetPosition().y();
           
    G4V3DNucleus* theProjectileNucleus = GetProjectileNucleus();
    theProjectileNucleus->StartLoop();

    G4Nucleon* Neighbour( 0 );
    while ( ( Neighbour = theProjectileNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
      if ( ! Neighbour->AreYouHit() ) {
        G4double impact2= sqr( XofWoundedNucleon - Neighbour->GetPosition().x() ) +
                          sqr( YofWoundedNucleon - Neighbour->GetPosition().y() );

        if ( G4UniformRand() < theParameters->GetCofNuclearDestructionPr() * 
                               G4Exp( -impact2 / theParameters->GetR2ofNuclearDestruction() )
           ) {
          // The neighbour nucleon is involved in the reggeon cascade
          TheInvolvedNucleonsOfProjectile[ NumberOfInvolvedNucleonsOfProjectile ] = Neighbour;
          NumberOfInvolvedNucleonsOfProjectile++;

          G4VSplitableHadron* projectileSplitable; 
          projectileSplitable = new G4DiffractiveSplitableHadron( *Neighbour ); 

          Neighbour->Hit( projectileSplitable );
          projectileSplitable->SetTimeOfCreation( CreationTime );
          projectileSplitable->SetStatus( 3 );
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

  G4bool isProjectileNucleus = false;
  if ( GetProjectileNucleus() ) isProjectileNucleus = true;

  #ifdef debugPutOnMassShell
  G4cout << "PutOnMassShell start " << G4endl;
  if ( isProjectileNucleus ) {
    G4cout << "PutOnMassShell for Nucleus_Nucleus " << G4endl;
  }
  #endif

  G4LorentzVector Pprojectile( theProjectile.GetMomentum(), theProjectile.GetTotalEnergy() );
  if ( Pprojectile.z() < 0.0 ) return false;

  G4bool isOk = true;
  
  G4LorentzVector Ptarget( 0.0, 0.0, 0.0, 0.0 );
  G4LorentzVector PtargetResidual( 0.0, 0.0, 0.0, 0.0 );
  G4double SumMasses = 0.0;
  G4V3DNucleus* theTargetNucleus = GetTargetNucleus();
  G4double TargetResidualMass = 0.0; 

  #ifdef debugPutOnMassShell
  G4cout << "Target : ";
  #endif
  isOk = ComputeNucleusProperties( theTargetNucleus, Ptarget, PtargetResidual, SumMasses,
                                   TargetResidualExcitationEnergy, TargetResidualMass,
                                   TargetResidualMassNumber, TargetResidualCharge );
  if ( ! isOk ) return false;

  G4double Mprojectile  = 0.0;
  G4double M2projectile = 0.0;
  G4LorentzVector Pproj( 0.0, 0.0, 0.0, 0.0 );
  G4LorentzVector PprojResidual( 0.0, 0.0, 0.0, 0.0 );
  G4V3DNucleus* thePrNucleus = GetProjectileNucleus();
  G4double PrResidualMass = 0.0;

  if ( ! isProjectileNucleus ) {  // hadron-nucleus collision
    Mprojectile  = Pprojectile.mag();
    M2projectile = Pprojectile.mag2();
    SumMasses += Mprojectile + 20.0*MeV;
  } else {  // nucleus-nucleus or antinucleus-nucleus collision
    #ifdef debugPutOnMassShell
    G4cout << "Projectile : ";
    #endif
    isOk = ComputeNucleusProperties( thePrNucleus, Pproj, PprojResidual, SumMasses,
                                     ProjectileResidualExcitationEnergy, PrResidualMass,
                                     ProjectileResidualMassNumber, ProjectileResidualCharge );
    if ( ! isOk ) return false;
  }

  G4LorentzVector Psum = Pprojectile + Ptarget;   
  G4double SqrtS = Psum.mag();
  G4double     S = Psum.mag2();

  #ifdef debugPutOnMassShell
  G4cout << "Psum " << Psum/GeV << " GeV" << G4endl << "SqrtS " << SqrtS/GeV << " GeV" << G4endl
         << "SumMasses, PrResidualMass and TargetResidualMass " << SumMasses/GeV << " " 
         << PrResidualMass/GeV << " " << TargetResidualMass/GeV << " GeV" << G4endl;
  #endif

  if ( SqrtS < SumMasses ) return false;  // It is impossible to simulate after putting nuclear nucleons on mass-shell

  // Try to consider also the excitation energy of the residual nucleus, if this is
  // possible, with the available energy; otherwise, set the excitation energy to zero.
  G4double savedSumMasses = SumMasses;
  if ( isProjectileNucleus ) {
    SumMasses -= std::sqrt( sqr( PrResidualMass ) + PprojResidual.perp2() );
    SumMasses += std::sqrt( sqr( PrResidualMass + ProjectileResidualExcitationEnergy ) 
                            + PprojResidual.perp2() ); 
  }
  SumMasses -= std::sqrt( sqr( TargetResidualMass ) + PtargetResidual.perp2() );
  SumMasses += std::sqrt( sqr( TargetResidualMass + TargetResidualExcitationEnergy )
                          + PtargetResidual.perp2() );

  if ( SqrtS < SumMasses ) {
    SumMasses = savedSumMasses;
    if ( isProjectileNucleus ) ProjectileResidualExcitationEnergy = 0.0;
    TargetResidualExcitationEnergy = 0.0;
  }

  TargetResidualMass += TargetResidualExcitationEnergy;
  if ( isProjectileNucleus ) PrResidualMass += ProjectileResidualExcitationEnergy;

  #ifdef debugPutOnMassShell
  if ( isProjectileNucleus ) {
    G4cout << "PrResidualMass ProjResidualExcitationEnergy " << PrResidualMass/GeV << " "
	   << ProjectileResidualExcitationEnergy << " MeV" << G4endl;
  }
  G4cout << "TargetResidualMass TargetResidualExcitationEnergy " << TargetResidualMass/GeV << " " 
         << TargetResidualExcitationEnergy << " MeV" << G4endl
         << "Sum masses " << SumMasses/GeV << G4endl;
  #endif

  // Sampling of nucleons what can transfer to delta-isobars
  if ( isProjectileNucleus  &&  thePrNucleus->GetMassNumber() != 1 ) {
      isOk = GenerateDeltaIsobar( SqrtS, NumberOfInvolvedNucleonsOfProjectile,
                                  TheInvolvedNucleonsOfProjectile, SumMasses );       
  }
  if ( theTargetNucleus->GetMassNumber() != 1 ) {
    isOk = isOk  &&  GenerateDeltaIsobar( SqrtS, NumberOfInvolvedNucleonsOfTarget,
                                          TheInvolvedNucleonsOfTarget, SumMasses );
  }
  if ( ! isOk ) return false;

  // Now we know that it is kinematically possible to produce a final state made
  // of the involved nucleons (or corresponding delta-isobars) and a residual nucleus.
  // We have to sample the kinematical variables which will allow to define the 4-momenta
  // of the final state. The sampled kinematical variables refer to the center-of-mass frame.
  // Notice that the sampling of the transverse momentum corresponds to take into account
  // Fermi motion.

  G4LorentzRotation toCms( -1*Psum.boostVector() );
  G4LorentzVector Ptmp = toCms*Pprojectile;
  if ( Ptmp.pz() <= 0.0 ) return false;  // "String" moving backwards in c.m.s., abort collision!

  G4LorentzRotation toLab( toCms.inverse() );
  
  G4double YprojectileNucleus = 0.0;
  if ( isProjectileNucleus ) {
    Ptmp = toCms*Pproj;                      
    YprojectileNucleus = Ptmp.rapidity();
  }
  Ptmp = toCms*Ptarget;                      
  G4double YtargetNucleus = Ptmp.rapidity();

  // Ascribing of the involved nucleons Pt and Xminus
  G4double DcorP = 0.0;
  if ( isProjectileNucleus ) DcorP = theParameters->GetDofNuclearDestruction() / thePrNucleus->GetMassNumber();
  G4double DcorT       = theParameters->GetDofNuclearDestruction() / theTargetNucleus->GetMassNumber();
  G4double AveragePt2  = theParameters->GetPt2ofNuclearDestruction();
  G4double maxPtSquare = theParameters->GetMaxPt2ofNuclearDestruction();

  #ifdef debugPutOnMassShell
  if ( isProjectileNucleus ) {
    G4cout << "Y projectileNucleus " << YprojectileNucleus << G4endl;
  }
  G4cout << "Y targetNucleus     " << YtargetNucleus << G4endl 
         << "Dcor " << theParameters->GetDofNuclearDestruction()
         << " DcorP DcorT " << DcorP << " " << DcorT << " AveragePt2 " << AveragePt2 << G4endl;
  #endif

  G4double M2proj = M2projectile;  // Initialization needed only for hadron-nucleus collisions
  G4double WplusProjectile = 0.0;
  G4double M2target = 0.0;
  G4double WminusTarget = 0.0;
  G4int NumberOfTries = 0;
  G4double ScaleFactor = 2.0;
  G4bool OuterSuccess = true;

  const G4int maxNumberOfLoops = 1000;
  G4int loopCounter = 0;
  do {  // while ( ! OuterSuccess )
    OuterSuccess = true;
    const G4int maxNumberOfInnerLoops = 10000;
    do {  // while ( SqrtS < Mprojectile + std::sqrt( M2target ) )
      NumberOfTries++;
      if ( NumberOfTries == 100*(NumberOfTries/100) ) {
        // After many tries, it is convenient to reduce the values of DcorP, DcorT and
        // AveragePt2, so that the sampled momenta (respectively, pz, and pt) of the
	// involved nucleons (or corresponding delta-isomers) are smaller, and therefore
        // it is more likely to satisfy the momentum conservation.
        ScaleFactor /= 2.0;
        DcorP       *= ScaleFactor;
        DcorT       *= ScaleFactor;
        AveragePt2  *= ScaleFactor;
      }
      if ( isProjectileNucleus ) {
        // Sampling of kinematical properties of projectile nucleons
        isOk = SamplingNucleonKinematics( AveragePt2, maxPtSquare, DcorP, 
                                          thePrNucleus, PprojResidual, 
                                          PrResidualMass, ProjectileResidualMassNumber,
                                          NumberOfInvolvedNucleonsOfProjectile, 
                                          TheInvolvedNucleonsOfProjectile, M2proj );
      }
      // Sampling of kinematical properties of target nucleons
      isOk = isOk  &&  SamplingNucleonKinematics( AveragePt2, maxPtSquare, DcorT, 
                                                  theTargetNucleus, PtargetResidual, 
                                                  TargetResidualMass, TargetResidualMassNumber,
                                                  NumberOfInvolvedNucleonsOfTarget, 
                                                  TheInvolvedNucleonsOfTarget, M2target );
      #ifdef debugPutOnMassShell
      G4cout << "SqrtS, Mp+Mt, Mp, Mt " << SqrtS/GeV << " " 
             << ( std::sqrt( M2proj ) + std::sqrt( M2target) )/GeV << " "
             << std::sqrt( M2proj )/GeV << " " << std::sqrt( M2target )/GeV << G4endl;
      #endif
      if ( ! isOk ) return false;
    } while ( ( SqrtS < std::sqrt( M2proj ) + std::sqrt( M2target ) ) &&
              NumberOfTries < maxNumberOfInnerLoops );  /* Loop checking, 10.08.2015, A.Ribon */
    if ( NumberOfTries >= maxNumberOfInnerLoops ) {
      #ifdef debugPutOnMassShell
      G4cout << "BAD situation: forced exit of the inner while loop!" << G4endl;
      #endif
      return false;
    }
    if ( isProjectileNucleus ) {
      isOk = CheckKinematics( S, SqrtS, M2proj, M2target, YprojectileNucleus, true, 
                              NumberOfInvolvedNucleonsOfProjectile, 
                              TheInvolvedNucleonsOfProjectile,
                              WminusTarget, WplusProjectile, OuterSuccess );
    }
    isOk = isOk  &&  CheckKinematics( S, SqrtS, M2proj, M2target, YtargetNucleus, false, 
                                      NumberOfInvolvedNucleonsOfTarget, TheInvolvedNucleonsOfTarget,
                                      WminusTarget, WplusProjectile, OuterSuccess );
    if ( ! isOk ) return false;
  } while ( ( ! OuterSuccess ) &&
            ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */
  if ( loopCounter >= maxNumberOfLoops ) {
    #ifdef debugPutOnMassShell
    G4cout << "BAD situation: forced exit of the while loop!" << G4endl;
    #endif
    return false;
  }

  // Now the sampling is completed, and we can determine the kinematics of the
  // whole system. This is done first in the center-of-mass frame, and then it is boosted
  // to the lab frame. The transverse momentum of the residual nucleus is determined as
  // the recoil of each hadron (nucleon or delta) which is emitted, i.e. in such a way
  // to conserve (by construction) the transverse momentum.

  if ( ! isProjectileNucleus ) {  // hadron-nucleus collision

    G4double Pzprojectile = WplusProjectile/2.0 - M2projectile/2.0/WplusProjectile;
    G4double Eprojectile  = WplusProjectile/2.0 + M2projectile/2.0/WplusProjectile;
    Pprojectile.setPz( Pzprojectile ); 
    Pprojectile.setE( Eprojectile );

    #ifdef debugPutOnMassShell
    G4cout << "Proj after in CMS " << Pprojectile << G4endl;
    #endif

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

  } else {  // nucleus-nucleus or antinucleus-nucleus collision

    isOk = FinalizeKinematics( WplusProjectile, true, toLab, PrResidualMass, 
                               ProjectileResidualMassNumber, NumberOfInvolvedNucleonsOfProjectile,
                               TheInvolvedNucleonsOfProjectile, ProjectileResidual4Momentum );

    #ifdef debugPutOnMassShell
    G4cout << "Projectile Residual4Momentum in CMS " << ProjectileResidual4Momentum << G4endl;
    #endif

    if ( ! isOk ) return false;

    ProjectileResidual4Momentum.transform( toLab );

    #ifdef debugPutOnMassShell
    G4cout << "Projectile Residual4Momentum in Lab " << ProjectileResidual4Momentum << G4endl;
    #endif

  }

  isOk = FinalizeKinematics( WminusTarget, false, toLab, TargetResidualMass, 
                             TargetResidualMassNumber, NumberOfInvolvedNucleonsOfTarget,
                             TheInvolvedNucleonsOfTarget, TargetResidual4Momentum );

  #ifdef debugPutOnMassShell
  G4cout << "Target Residual4Momentum in CMS " << TargetResidual4Momentum << G4endl;
  #endif

  if ( ! isOk ) return false;

  TargetResidual4Momentum.transform( toLab );

  #ifdef debugPutOnMassShell
  G4cout << "Target Residual4Momentum in Lab " << TargetResidual4Momentum << G4endl;
  #endif

  return true;

}


//============================================================================

G4bool G4FTFModel::ExciteParticipants() {

  #ifdef debugBuildString
  G4cout << "G4FTFModel::ExciteParticipants() " << G4endl;
  #endif

  G4bool Success( false );
  G4int MaxNumOfInelCollisions = G4int( theParameters->GetMaxNumberOfCollisions() );
  if ( MaxNumOfInelCollisions > 0 ) {  //  Plab > Pbound, normal application of FTF is possible
    G4double ProbMaxNumber = theParameters->GetMaxNumberOfCollisions() - MaxNumOfInelCollisions;
    if ( G4UniformRand() < ProbMaxNumber ) MaxNumOfInelCollisions++;
  } else { 
    // Plab < Pbound, normal application of FTF is impossible,low energy corrections applied
    MaxNumOfInelCollisions = 1;
  }

  #ifdef debugBuildString
  G4cout << "MaxNumOfInelCollisions per hadron/nucleon " << MaxNumOfInelCollisions << G4endl;
  #endif

  G4int CurrentInteraction( 0 );
  theParticipants.StartLoop();

  G4bool InnerSuccess( true );
  while ( theParticipants.Next() ) {  /* Loop checking, 10.08.2015, A.Ribon */
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
         InnerSuccess = theElastic->ElasticScattering( projectile, target, theParameters );
      } else if ( G4UniformRand() > theParameters->GetProbabilityOfAnnihilation() ) { 
        // Inelastic scattering

        #ifdef debugBuildString
        G4cout << "Inelastic interaction" << G4endl
               << "MaxNumOfInelCollisions per hadron/nucleon " << MaxNumOfInelCollisions << G4endl;
        #endif

        if ( ! HighEnergyInter ) {
          G4bool Annihilation = false;
          G4bool Result = AdjustNucleons( projectile, ProjectileNucleon, target,
                                          TargetNucleon, Annihilation );
          if ( ! Result ) continue;
        } 
        if ( G4UniformRand() < 
             ( 1.0 - target->GetSoftCollisionCount()     / MaxNumOfInelCollisions )  *
             ( 1.0 - projectile->GetSoftCollisionCount() / MaxNumOfInelCollisions ) ) {
          //if ( ! HighEnergyInter ) {
          //  G4bool Annihilation = false;
          //  G4bool Result = AdjustNucleons( projectile, ProjectileNucleon, target,
          //                                  TargetNucleon, Annihilation );
          //  if ( ! Result ) continue;
          //} 
          if ( theExcitation->ExciteParticipants( projectile, target, theParameters, theElastic ) ) {
            InnerSuccess = true;
            NumberOfNNcollisions++;
            #ifdef debugBuildString
            G4cout << "FTF excitation Successfull " << G4endl;
            // G4cout << "After  pro " << projectile->Get4Momentum() << " " 
            //        << projectile->Get4Momentum().mag() << G4endl
            //        << "After  tar " << target->Get4Momentum() << " "
            //        << target->Get4Momentum().mag() << G4endl;
            #endif
          } else {
            InnerSuccess = theElastic->ElasticScattering( projectile, target, theParameters );
            #ifdef debugBuildString
            G4cout << "FTF excitation Non InnerSuccess of Elastic scattering " 
                   << InnerSuccess << G4endl;
            #endif
          }
        } else {  // The inelastic interactition was rejected -> elastic scattering
          #ifdef debugBuildString
          G4cout << "Elastic scat. at rejection inelastic scattering" << G4endl;
          #endif
          //if ( ! HighEnergyInter ) {
          //  G4bool Annihilation = false;
          //  G4bool Result = AdjustNucleons( projectile, ProjectileNucleon, target,
          //                                  TargetNucleon, Annihilation );
          //  if ( ! Result) continue;
          //} 
          InnerSuccess = theElastic->ElasticScattering( projectile, target, theParameters );
        } 
      } else {  // Annihilation

        #ifdef debugBuildString
        G4cout << "Annihilation" << G4endl;
        #endif

        // At last, annihilation
        if ( ! HighEnergyInter ) {
          G4bool Annihilation = true;
          G4bool Result = AdjustNucleons( projectile, ProjectileNucleon, target,
                                          TargetNucleon, Annihilation );
          if ( ! Result ) continue;
        }

        G4VSplitableHadron* AdditionalString = 0;
        if ( theAnnihilation->Annihilate( projectile, target, AdditionalString, theParameters ) ) {
          InnerSuccess = true;
          #ifdef debugBuildString
          G4cout << "Annihilation successfull. " << "*AdditionalString " 
                 << AdditionalString << G4endl;
          //G4cout << "After  pro " << projectile->Get4Momentum() << G4endl;
          //G4cout << "After  tar " << target->Get4Momentum() << G4endl;
          #endif

          if ( AdditionalString != 0 ) theAdditionalString.push_back( AdditionalString );

          NumberOfNNcollisions++;

          // Skipping possible interactions of the annihilated nucleons 
          while ( theParticipants.Next() ) {   /* Loop checking, 10.08.2015, A.Ribon */  
            G4InteractionContent& acollision = theParticipants.GetInteraction();
            G4VSplitableHadron* NextProjectileNucleon = acollision.GetProjectile();
            G4VSplitableHadron* NextTargetNucleon = acollision.GetTarget();
            if ( projectile == NextProjectileNucleon  ||  target == NextTargetNucleon ) {
              acollision.SetStatus( 0 );
            }
          }

          // Continue the interactions
          theParticipants.StartLoop(); 
          for ( G4int i = 0; i < CurrentInteraction; ++i ) theParticipants.Next();
	  
          /*
          if ( target->GetStatus() == 4 ) {
            // Skipping possible interactions of the annihilated nucleons 
            while ( theParticipants.Next() ) {    
              G4InteractionContent& acollision = theParticipants.GetInteraction();
              G4VSplitableHadron* NextProjectileNucleon = acollision.GetProjectile();
              G4VSplitableHadron* NextTargetNucleon = acollision.GetTarget();
              if ( target == NextTargetNucleon ) { acollision.SetStatus( 0 ); }
            }
          }
          theParticipants.StartLoop(); 
          for ( G4int I = 0; I < CurrentInteraction; ++I ) theParticipants.Next();
          */

        }
      } 
    }

    if( InnerSuccess ) Success = true;

    #ifdef debugBuildString
    G4cout << "----------------------------- Final properties " << G4endl
           << "projectile->GetStatus target->GetStatus " << projectile->GetStatus() 
           << " " << target->GetStatus() << G4endl << "projectile->GetSoftC  target->GetSoftC  "
           << projectile->GetSoftCollisionCount() << " " << target->GetSoftCollisionCount()
           << G4endl << "ExciteParticipants() Success? " << Success << G4endl;
    #endif

  }  // end of while ( theParticipants.Next() )

  return Success;
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
         << "Collis. pr tr " << SelectedAntiBaryon->GetSoftCollisionCount() << " "
         << SelectedTargetNucleon->GetSoftCollisionCount() << G4endl;
  #endif

  if ( SelectedAntiBaryon->GetSoftCollisionCount() != 0  && 
       SelectedTargetNucleon->GetSoftCollisionCount() != 0 ) {
    return true; // Selected hadrons were adjusted before.   
  }

  G4int interactionCase = 0;
  if (    ( ! GetProjectileNucleus()  &&
            SelectedAntiBaryon->GetSoftCollisionCount() == 0  &&
            SelectedTargetNucleon->GetSoftCollisionCount() == 0 )
       ||
          ( SelectedAntiBaryon->GetSoftCollisionCount() != 0  && 
            SelectedTargetNucleon->GetSoftCollisionCount() == 0 ) ) {
    // The case of hadron-nucleus interactions, or
    // the case when projectile nuclear nucleon participated in
    // a collision, but target nucleon did not participate.
    interactionCase = 1;
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
  } else if ( SelectedAntiBaryon->GetSoftCollisionCount() == 0  && 
              SelectedTargetNucleon->GetSoftCollisionCount() != 0 ) {
    // It is assumed that in the case there is ProjectileResidualNucleus
    interactionCase = 2;
    #ifdef debugAdjust
    G4cout << "case 2,  prcol=0 trcol#0" << G4endl;
    #endif
    if ( ProjectileResidualMassNumber < 1 ) {
      return false;
    }
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
  } else {  // It has to be a nucleus-nucleus interaction
    interactionCase = 3;
    #ifdef debugAdjust
    G4cout << "case 3,  prcol=0 trcol=0" << G4endl;
    #endif
    if ( ! GetProjectileNucleus() ) {
      return false;
    }
    #ifdef debugAdjust
    G4cout << "Proj res Init " << ProjectileResidual4Momentum << G4endl
           << "Targ res Init " << TargetResidual4Momentum << G4endl
           << "ProjectileResidualMassNumber ProjectileResidualCharge (ProjectileResidualLambdaNumber)" 
           << ProjectileResidualMassNumber << " " << ProjectileResidualCharge
           << " (" << ProjectileResidualLambdaNumber << ") " << G4endl
           << "TargetResidualMassNumber TargetResidualCharge " << TargetResidualMassNumber
           << " " << TargetResidualCharge << G4endl;
    #endif
  }

  CommonVariables common;
  G4int returnCode = AdjustNucleonsAlgorithm_beforeSampling( interactionCase, SelectedAntiBaryon,
                                                             ProjectileNucleon, SelectedTargetNucleon,
                                                             TargetNucleon, Annihilation, common );
  G4bool returnResult = false;
  if ( returnCode == 0 ) {
    returnResult = true;  // Successfully ended: no need of extra work
  } else if ( returnCode == 1 ) {
    // The part before sampling has been successfully completed: now try the sampling
    returnResult = AdjustNucleonsAlgorithm_Sampling( interactionCase, common );
    if ( returnResult ) {  // The sampling has completed successfully: do the last part
      AdjustNucleonsAlgorithm_afterSampling( interactionCase, SelectedAntiBaryon, 
                                             SelectedTargetNucleon, common ); 
    }
  }

  return returnResult;
}

//-------------------------------------------------------------------

G4int G4FTFModel::AdjustNucleonsAlgorithm_beforeSampling( G4int interactionCase,
                                                          G4VSplitableHadron* SelectedAntiBaryon,
                                                          G4Nucleon* ProjectileNucleon,
                                                          G4VSplitableHadron* SelectedTargetNucleon,
                                                          G4Nucleon* TargetNucleon,
                                                          G4bool Annihilation,
                                                          G4FTFModel::CommonVariables& common ) {
  // First of the three utility methods used only by AdjustNucleons: prepare for sampling.
  // This method returns an integer code - instead of a boolean, with the following meaning:
  //   "0" : successfully ended and nothing else needs to be done (i.e. no sampling);
  //   "1" : successfully completed, but the work needs to be continued, i.e. try to sample;
  //  "99" : unsuccessfully ended, nothing else can be done.
  G4int returnCode = 99;

  G4double ExcitationEnergyPerWoundedNucleon = theParameters->GetExcitationEnergyPerWoundedNucleon();
 
  // some checks and initializations
  if ( interactionCase == 1 ) {
    common.Psum = SelectedAntiBaryon->Get4Momentum() + TargetResidual4Momentum;
    #ifdef debugAdjust
    G4cout << "Targ res Init " << TargetResidual4Momentum << G4endl;
    #endif
    common.Pprojectile = SelectedAntiBaryon->Get4Momentum();
  } else if ( interactionCase == 2 ) {
    common.Psum = ProjectileResidual4Momentum + SelectedTargetNucleon->Get4Momentum();
    common.Pprojectile = ProjectileResidual4Momentum;
  } else if ( interactionCase == 3 ) {
    common.Psum = ProjectileResidual4Momentum + TargetResidual4Momentum; 
    common.Pprojectile = ProjectileResidual4Momentum;
  }

  // transform momenta to cms and then rotate parallel to z axis
  common.toCms = G4LorentzRotation( -1*common.Psum.boostVector() ); 
  common.Ptmp = common.toCms * common.Pprojectile;
  common.toCms.rotateZ( -1*common.Ptmp.phi() );
  common.toCms.rotateY( -1*common.Ptmp.theta() );
  common.Pprojectile.transform( common.toCms );
  common.toLab = common.toCms.inverse();
  common.SqrtS = common.Psum.mag();
  common.S = sqr( common.SqrtS );

  // get properties of the target residual and/or projectile residual
  G4bool Stopping = false;
  if ( interactionCase == 1 ) {
    common.TResidualMassNumber = TargetResidualMassNumber - 1;
    common.TResidualCharge =   TargetResidualCharge 
                             - G4int( TargetNucleon->GetDefinition()->GetPDGCharge() );
    common.TResidualExcitationEnergy =   TargetResidualExcitationEnergy
                                       - ExcitationEnergyPerWoundedNucleon*G4Log( G4UniformRand() );
    if ( common.TResidualMassNumber <= 1 ) {
      common.TResidualExcitationEnergy = 0.0;
    }
    if ( common.TResidualMassNumber != 0 ) {
      common.TResidualMass = G4ParticleTable::GetParticleTable()->GetIonTable()
                             ->GetIonMass( common.TResidualCharge, common.TResidualMassNumber );
    }
    common.TNucleonMass = TargetNucleon->GetDefinition()->GetPDGMass();
    common.SumMasses =   SelectedAntiBaryon->Get4Momentum().mag() + common.TNucleonMass 
                       + common.TResidualMass;
    #ifdef debugAdjust
    G4cout << "Annihilation " << Annihilation << G4endl;
    #endif
  } else if ( interactionCase == 2 ) {
    common.Ptarget = common.toCms * SelectedTargetNucleon->Get4Momentum();
    common.TResidualMassNumber = ProjectileResidualMassNumber - 1;
    common.TResidualCharge =   ProjectileResidualCharge 
                             - std::abs( G4int(ProjectileNucleon->GetDefinition()->GetPDGCharge()) );
    common.TResidualExcitationEnergy =   ProjectileResidualExcitationEnergy 
                                       - ExcitationEnergyPerWoundedNucleon*G4Log( G4UniformRand() );
    if ( common.TResidualMassNumber <= 1 ) {
      common.TResidualExcitationEnergy = 0.0;
    }
    if ( common.TResidualMassNumber != 0 ) {
      common.TResidualMass = G4ParticleTable::GetParticleTable()->GetIonTable()
                             ->GetIonMass( common.TResidualCharge, common.TResidualMassNumber );
    }
    common.TNucleonMass = ProjectileNucleon->GetDefinition()->GetPDGMass();
    common.SumMasses =   SelectedTargetNucleon->Get4Momentum().mag() + common.TNucleonMass 
                       + common.TResidualMass;
    #ifdef debugAdjust
    G4cout << "SelectedTN.mag() PNMass + PResidualMass " 
           << SelectedTargetNucleon->Get4Momentum().mag() << " " 
           << common.TNucleonMass << " " << common.TResidualMass << G4endl;
    #endif
  } else if ( interactionCase == 3 ) {
    common.Ptarget = common.toCms * TargetResidual4Momentum;
    common.PResidualMassNumber = ProjectileResidualMassNumber - 1;
    common.PResidualCharge =   ProjectileResidualCharge 
                             - std::abs( G4int(ProjectileNucleon->GetDefinition()->GetPDGCharge()) );
    common.PResidualLambdaNumber = ProjectileResidualLambdaNumber;
    if ( ProjectileNucleon->GetDefinition() == G4Lambda::Definition()  ||
	 ProjectileNucleon->GetDefinition() == G4AntiLambda::Definition() ) {
      --common.PResidualLambdaNumber;  
    }
    common.PResidualExcitationEnergy =   ProjectileResidualExcitationEnergy 
                                       - ExcitationEnergyPerWoundedNucleon*G4Log( G4UniformRand() );
    if ( common.PResidualMassNumber <= 1 ) {
      common.PResidualExcitationEnergy = 0.0;
    }
    if ( common.PResidualMassNumber != 0 ) {
      if ( common.PResidualMassNumber == 1 ) {
        if ( std::abs( common.PResidualCharge ) == 1 ) {
          common.PResidualMass = G4Proton::Definition()->GetPDGMass();       
        } else if ( common.PResidualLambdaNumber == 1 ) {
          common.PResidualMass = G4Lambda::Definition()->GetPDGMass();
        } else {
          common.PResidualMass = G4Neutron::Definition()->GetPDGMass();
        }
      } else {
        if ( common.PResidualLambdaNumber > 0 ) {
          if ( common.PResidualMassNumber == 2 ) {
            common.PResidualMass = G4Lambda::Definition()->GetPDGMass();
	    if ( std::abs( common.PResidualCharge ) == 1 ) {   // lambda + proton
	      common.PResidualMass += G4Proton::Definition()->GetPDGMass();
	    } else if ( common.PResidualLambdaNumber == 1 ) {  // lambda + neutron
	      common.PResidualMass += G4Neutron::Definition()->GetPDGMass();
	    } else {                                           // lambda + lambda
	      common.PResidualMass += G4Lambda::Definition()->GetPDGMass();
	    }
	  } else {
	    common.PResidualMass = G4HyperNucleiProperties::GetNuclearMass( common.PResidualMassNumber,
	  								    std::abs( common.PResidualCharge ),
									    common.PResidualLambdaNumber );
	  }
        } else {
	  common.PResidualMass = G4ParticleTable::GetParticleTable()->GetIonTable()->
	                         GetIonMass( std::abs( common.PResidualCharge ), common.PResidualMassNumber );
        }
      }
    }
    common.PNucleonMass = ProjectileNucleon->GetDefinition()->GetPDGMass();  // On-shell (anti-)nucleon mass 
    common.TResidualMassNumber = TargetResidualMassNumber - 1;
    common.TResidualCharge =   TargetResidualCharge
                             - G4int( TargetNucleon->GetDefinition()->GetPDGCharge() );
    common.TResidualExcitationEnergy =   TargetResidualExcitationEnergy 
                                       - ExcitationEnergyPerWoundedNucleon*G4Log( G4UniformRand() );
    if ( common.TResidualMassNumber <= 1 ) {
      common.TResidualExcitationEnergy = 0.0;
    }
    if ( common.TResidualMassNumber != 0 ) {
      common.TResidualMass = G4ParticleTable::GetParticleTable()->GetIonTable()->
                             GetIonMass( common.TResidualCharge, common.TResidualMassNumber );
    }
    common.TNucleonMass = TargetNucleon->GetDefinition()->GetPDGMass();  // On-shell nucleon mass
    common.SumMasses =   common.PNucleonMass + common.PResidualMass + common.TNucleonMass 
                       + common.TResidualMass;
    #ifdef debugAdjust
    G4cout << "PNucleonMass PResidualMass TNucleonMass TResidualMass " << common.PNucleonMass 
           << " " << common.PResidualMass << " " << common.TNucleonMass << " " 
           << common.TResidualMass << G4endl
           << "PResidualExcitationEnergy " << common.PResidualExcitationEnergy << G4endl
           << "TResidualExcitationEnergy " << common.TResidualExcitationEnergy << G4endl;
    #endif
  }  // End-if on interactionCase

  if ( ! Annihilation ) {
    if ( common.SqrtS < common.SumMasses ) {
      #ifdef debugAdjust
      G4cout << "SqrtS < SumMasses " << common.SqrtS << " " << common.SumMasses << G4endl;
      #endif
      return returnCode;  // Unsuccessfully ended, nothing else can be done
    } 
    if ( interactionCase == 1  ||  interactionCase == 2 ) {
      if ( common.SqrtS < common.SumMasses + common.TResidualExcitationEnergy ) {
        #ifdef debugAdjust
        G4cout << "TResidualExcitationEnergy : before " << common.TResidualExcitationEnergy << G4endl;
        #endif
        common.TResidualExcitationEnergy = common.SqrtS - common.SumMasses;
        #ifdef debugAdjust
        G4cout << "TResidualExcitationEnergy : after " << common.TResidualExcitationEnergy << G4endl;
        #endif
        Stopping = true;
        return returnCode;  // unsuccessfully ended, nothing else can be done
      }
    } else if ( interactionCase == 3 ) {
      #ifdef debugAdjust
      G4cout << "SqrtS < SumMasses + PResidualExcitationEnergy + TResidualExcitationEnergy "
             << common.SqrtS << " " << common.SumMasses + common.PResidualExcitationEnergy + common.TResidualExcitationEnergy
             << G4endl;
      #endif
      if ( common.SqrtS <   common.SumMasses + common.PResidualExcitationEnergy
                          + common.TResidualExcitationEnergy ) { 
        Stopping = true;
        if ( common.PResidualExcitationEnergy <= 0.0 ) {
          common.TResidualExcitationEnergy = common.SqrtS - common.SumMasses;
        } else if ( common.TResidualExcitationEnergy <= 0.0 ) {
          common.PResidualExcitationEnergy = common.SqrtS - common.SumMasses;
        } else {
          G4double Fraction = ( common.SqrtS - common.SumMasses )
            /  ( common.PResidualExcitationEnergy + common.TResidualExcitationEnergy );
          common.PResidualExcitationEnergy *= Fraction;
          common.TResidualExcitationEnergy *= Fraction;
        }
      }
    }
  }  // End-if on ! Annihilation

  if ( Annihilation ) {
    #ifdef debugAdjust
    G4cout << "SqrtS < SumMasses - TNucleonMass " << common.SqrtS << " " 
           << common.SumMasses - common.TNucleonMass << G4endl;
    #endif
    if ( common.SqrtS < common.SumMasses - common.TNucleonMass ) {
      return returnCode;  // unsuccessfully ended, nothing else can be done
    } 
    #ifdef debugAdjust
    G4cout << "SqrtS < SumMasses " << common.SqrtS << " " << common.SumMasses << G4endl;
    #endif
    if ( common.SqrtS < common.SumMasses ) {
      if ( interactionCase == 2  ||  interactionCase == 3 ) {
        common.TResidualExcitationEnergy = 0.0;
      }  
      common.TNucleonMass =   common.SqrtS - ( common.SumMasses - common.TNucleonMass )
                            - common.TResidualExcitationEnergy;  // Off-shell nucleon mass 
      #ifdef debugAdjust
      G4cout << "TNucleonMass " << common.TNucleonMass << G4endl;
      #endif
      common.SumMasses = common.SqrtS - common.TResidualExcitationEnergy;
      Stopping = true;
      #ifdef debugAdjust
      G4cout << "SqrtS < SumMasses " << common.SqrtS << " " << common.SumMasses << G4endl;
      #endif
    }
    if ( interactionCase == 1  ||  interactionCase == 2 ) {
      if ( common.SqrtS < common.SumMasses + common.TResidualExcitationEnergy ) {
        common.TResidualExcitationEnergy = common.SqrtS - common.SumMasses; 
        Stopping = true;
      }
    } else if ( interactionCase == 3 ) {
      if ( common.SqrtS <   common.SumMasses + common.PResidualExcitationEnergy
                          + common.TResidualExcitationEnergy ) {
        Stopping = true;
        if ( common.PResidualExcitationEnergy <= 0.0 ) {
          common.TResidualExcitationEnergy = common.SqrtS - common.SumMasses;
        } else if ( common.TResidualExcitationEnergy <= 0.0 ) {
          common.PResidualExcitationEnergy = common.SqrtS - common.SumMasses;
        } else {
          G4double Fraction = ( common.SqrtS - common.SumMasses ) / 
            ( common.PResidualExcitationEnergy + common.TResidualExcitationEnergy );
          common.PResidualExcitationEnergy *= Fraction;
          common.TResidualExcitationEnergy *= Fraction;
        }
      }
    }
  }  // End-if on Annihilation

  #ifdef debugAdjust
  G4cout << "Stopping " << Stopping << G4endl;
  #endif

  if ( Stopping ) {
    // All 3-momenta of particles = 0
    common.Ptmp.setPx( 0.0 ); common.Ptmp.setPy( 0.0 ); common.Ptmp.setPz( 0.0 );
    // New projectile
    if ( interactionCase == 1 ) {
      common.Ptmp.setE( SelectedAntiBaryon->Get4Momentum().mag() );
    } else if ( interactionCase == 2 ) {
      common.Ptmp.setE( common.TNucleonMass );
    } else if ( interactionCase == 3 ) {
      common.Ptmp.setE( common.PNucleonMass );
    }
    #ifdef debugAdjust
    G4cout << "Proj stop " << common.Ptmp << G4endl;
    #endif
    common.Pprojectile = common.Ptmp; 
    common.Pprojectile.transform( common.toLab );  // From center-of-mass to Lab frame
    //---AR-Jul2019 : To avoid unphysical projectile (anti-)fragments at rest, save the
    //                original momentum of the anti-baryon in the center-of-mass frame.
    G4LorentzVector saveSelectedAntiBaryon4Momentum = SelectedAntiBaryon->Get4Momentum();
    saveSelectedAntiBaryon4Momentum.transform( common.toCms );  // From Lab to center-of-mass frame
    //---
    SelectedAntiBaryon->Set4Momentum( common.Pprojectile );
    // New target nucleon
    if ( interactionCase == 1  ||  interactionCase == 3 ) {
      common.Ptmp.setE( common.TNucleonMass );
    } else if ( interactionCase == 2 ) {
      common.Ptmp.setE( SelectedTargetNucleon->Get4Momentum().mag() );
    }
    #ifdef debugAdjust
    G4cout << "Targ stop " << common.Ptmp << G4endl;
    #endif
    common.Ptarget = common.Ptmp; 
    common.Ptarget.transform( common.toLab );  // From center-of-mass to Lab frame
    //---AR-Jul2019 : To avoid unphysical target fragments at rest, save the original
    //                momentum of the target nucleon in the center-of-mass frame.
    G4LorentzVector saveSelectedTargetNucleon4Momentum = SelectedTargetNucleon->Get4Momentum();
    saveSelectedTargetNucleon4Momentum.transform( common.toCms );  // From Lab to center-of-mass frame
    //---
    SelectedTargetNucleon->Set4Momentum( common.Ptarget );
    // New target residual
    if ( interactionCase == 1  ||  interactionCase == 3 ) {
      common.Ptmp.setPx( 0.0 ); common.Ptmp.setPy( 0.0 ); common.Ptmp.setPz( 0.0 );
      TargetResidualMassNumber       = common.TResidualMassNumber;
      TargetResidualCharge           = common.TResidualCharge;
      TargetResidualExcitationEnergy = common.TResidualExcitationEnergy;
      //---AR-Jul2019 : To avoid unphysical target fragments at rest, use the saved
      //                original momentum of the target nucleon (instead of setting 0).
      //                This is a rough and simple approach!
      //common.Ptmp.setE( common.TResidualMass + TargetResidualExcitationEnergy );
      common.Ptmp.setPx( -saveSelectedTargetNucleon4Momentum.x() ); 
      common.Ptmp.setPy( -saveSelectedTargetNucleon4Momentum.y() ); 
      common.Ptmp.setPz( -saveSelectedTargetNucleon4Momentum.z() );
      common.Ptmp.setE( std::sqrt( sqr( common.TResidualMass + TargetResidualExcitationEnergy ) + common.Ptmp.vect().mag2() ) ); 
      //---
      #ifdef debugAdjust
      G4cout << "Targ Resi stop " << common.Ptmp << G4endl;
      #endif
      common.Ptmp.transform( common.toLab );  // From center-of-mass to Lab frame
      TargetResidual4Momentum = common.Ptmp;
    }
    // New projectile residual
    if ( interactionCase == 2  ||  interactionCase == 3 ) {
      common.Ptmp.setPx( 0.0 ); common.Ptmp.setPy( 0.0 ); common.Ptmp.setPz( 0.0 );
      if ( interactionCase == 2 ) {
        ProjectileResidualMassNumber       = common.TResidualMassNumber;
        ProjectileResidualCharge           = common.TResidualCharge;
	ProjectileResidualLambdaNumber     = 0;  // The target nucleus and its residual are never hypernuclei
        ProjectileResidualExcitationEnergy = common.TResidualExcitationEnergy;
        common.Ptmp.setE( common.TResidualMass + ProjectileResidualExcitationEnergy ); 
      } else {
        ProjectileResidualMassNumber       = common.PResidualMassNumber;
        ProjectileResidualCharge           = common.PResidualCharge;
	ProjectileResidualLambdaNumber     = common.PResidualLambdaNumber;
        ProjectileResidualExcitationEnergy = common.PResidualExcitationEnergy;
        //---AR-Jul2019 : To avoid unphysical projectile (anti-)fragments at rest, use the
        //                saved original momentum of the anti-baryon (instead of setting 0).
        //                This is a rough and simple approach!
        //common.Ptmp.setE( common.PResidualMass + ProjectileResidualExcitationEnergy );
        common.Ptmp.setPx( -saveSelectedAntiBaryon4Momentum.x() ); 
        common.Ptmp.setPy( -saveSelectedAntiBaryon4Momentum.y() ); 
        common.Ptmp.setPz( -saveSelectedAntiBaryon4Momentum.z() );
        common.Ptmp.setE( std::sqrt( sqr( common.PResidualMass + ProjectileResidualExcitationEnergy ) + common.Ptmp.vect().mag2() ) ); 
        //---
      }
      #ifdef debugAdjust
      G4cout << "Proj Resi stop " << common.Ptmp << G4endl;
      #endif
      common.Ptmp.transform( common.toLab );  // From center-of-mass to Lab frame
      ProjectileResidual4Momentum = common.Ptmp;
    }
    return returnCode = 0;  // successfully ended and nothing else needs to be done (i.e. no sampling)
  }  // End-if on Stopping

  // Initializations before sampling
  if ( interactionCase == 1 ) {
    common.Mprojectile  = common.Pprojectile.mag();
    common.M2projectile = common.Pprojectile.mag2();
    common.TResidual4Momentum = common.toCms * TargetResidual4Momentum;
    common.YtargetNucleus = common.TResidual4Momentum.rapidity();
    common.TResidualMass += common.TResidualExcitationEnergy;
  } else if ( interactionCase == 2 ) {
    common.Mtarget  = common.Ptarget.mag();
    common.M2target = common.Ptarget.mag2();
    common.TResidual4Momentum = common.toCms * ProjectileResidual4Momentum;
    common.YprojectileNucleus = common.TResidual4Momentum.rapidity();
    common.TResidualMass += common.TResidualExcitationEnergy;
  } else if ( interactionCase == 3 ) {
    common.PResidual4Momentum = common.toCms * ProjectileResidual4Momentum;
    common.YprojectileNucleus = common.PResidual4Momentum.rapidity();
    common.TResidual4Momentum = common.toCms*TargetResidual4Momentum;
    common.YtargetNucleus = common.TResidual4Momentum.rapidity();
    common.PResidualMass += common.PResidualExcitationEnergy;
    common.TResidualMass += common.TResidualExcitationEnergy;
  }
  #ifdef debugAdjust
  G4cout << "YprojectileNucleus " << common.YprojectileNucleus << G4endl;
  #endif

  return returnCode = 1;  // successfully completed, but the work needs to be continued, i.e. try to sample
}


//-------------------------------------------------------------------

G4bool G4FTFModel::AdjustNucleonsAlgorithm_Sampling( G4int interactionCase,
                                                     G4FTFModel::CommonVariables& common ) {
  // Second of the three utility methods used only by AdjustNucleons: do the sampling.
  // This method returns "false" if it fails to sample properly, else it returns "true".

  // Ascribing of the involved nucleons Pt and X 
  G4double Dcor = theParameters->GetDofNuclearDestruction();
  G4double DcorP = 0.0, DcorT = 0.0;
  if ( ProjectileResidualMassNumber != 0 ) DcorP = Dcor / G4double(ProjectileResidualMassNumber);
  if ( TargetResidualMassNumber != 0 )     DcorT = Dcor / G4double(TargetResidualMassNumber);
  G4double AveragePt2 = theParameters->GetPt2ofNuclearDestruction();
  G4double maxPtSquare = theParameters->GetMaxPt2ofNuclearDestruction();

  G4double ScaleFactor = 1.0;
  G4bool OuterSuccess = true;
  const G4int maxNumberOfLoops = 1000;
  const G4int maxNumberOfTries = 10000;
  G4int loopCounter = 0;
  G4int NumberOfTries = 0;
  do {  // Outmost do while loop
    OuterSuccess = true;
    G4bool loopCondition = false;
    do {  // Intermediate do while loop
      if ( NumberOfTries == 100*(NumberOfTries/100) ) { 
        // At large number of tries it would be better to reduce the values
        ScaleFactor /= 2.0;
        DcorP      *= ScaleFactor;
        DcorT      *= ScaleFactor;
        AveragePt2 *= ScaleFactor;
        #ifdef debugAdjust
        //G4cout << "NumberOfTries ScaleFactor " << NumberOfTries << " " << ScaleFactor << G4endl;
        #endif
      }

      // Some kinematics
      if ( interactionCase == 1 ) {
      } else if ( interactionCase == 2 ) {
        #ifdef debugAdjust
        G4cout << "ProjectileResidualMassNumber " << ProjectileResidualMassNumber << G4endl;
        #endif
        if ( ProjectileResidualMassNumber > 1 ) {
          common.PtNucleon = GaussianPt( AveragePt2, maxPtSquare );
        } else {
          common.PtNucleon = G4ThreeVector( 0.0, 0.0, 0.0 );
        }
        common.PtResidual = - common.PtNucleon;
        common.Mprojectile =   std::sqrt( sqr( common.TNucleonMass ) + common.PtNucleon.mag2() )
                             + std::sqrt( sqr( common.TResidualMass ) + common.PtResidual.mag2() );
        #ifdef debugAdjust
        G4cout << "SqrtS < Mtarget + Mprojectile " << common.SqrtS << "  " << common.Mtarget
               << " " << common.Mprojectile << " "  << common.Mtarget + common.Mprojectile << G4endl;
        #endif
        common.M2projectile = sqr( common.Mprojectile );
        if ( common.SqrtS < common.Mtarget + common.Mprojectile ) {
          OuterSuccess = false; 
          loopCondition = true;
          continue;
        }
      } else if ( interactionCase == 3 ) {
        if ( ProjectileResidualMassNumber > 1 ) {            
          common.PtNucleonP = GaussianPt( AveragePt2, maxPtSquare );
        } else {
          common.PtNucleonP = G4ThreeVector( 0.0, 0.0, 0.0 );
        }
        common.PtResidualP = - common.PtNucleonP;
        if ( TargetResidualMassNumber > 1 ) { 
          common.PtNucleonT = GaussianPt( AveragePt2, maxPtSquare );
        } else {
          common.PtNucleonT = G4ThreeVector( 0.0, 0.0, 0.0 );
        }
        common.PtResidualT = - common.PtNucleonT;
        common.Mprojectile =   std::sqrt( sqr( common.PNucleonMass )  + common.PtNucleonP.mag2() ) 
                             + std::sqrt( sqr( common.PResidualMass ) + common.PtResidualP.mag2() );
        common.M2projectile = sqr( common.Mprojectile );
        common.Mtarget =   std::sqrt( sqr( common.TNucleonMass )  + common.PtNucleonT.mag2() )
                         + std::sqrt( sqr( common.TResidualMass ) + common.PtResidualT.mag2() );
        common.M2target = sqr( common.Mtarget );
        if ( common.SqrtS < common.Mprojectile + common.Mtarget ) {
          OuterSuccess = false; 
          loopCondition = true;
          continue;
        }
      }  // End-if on interactionCase

      G4int numberOfTimesExecuteInnerLoop = 1;
      if ( interactionCase == 3 ) numberOfTimesExecuteInnerLoop = 2;
      for ( G4int iExecute = 0; iExecute < numberOfTimesExecuteInnerLoop; iExecute++ ) {

        G4bool InnerSuccess = true;
        G4bool isTargetToBeHandled = ( interactionCase == 1 || 
                                       ( interactionCase == 3 && iExecute == 1 ) );
        G4bool condition = false;
        if ( isTargetToBeHandled ) {
          condition = ( TargetResidualMassNumber > 1 );
	} else {  // Projectile to be handled
          condition = ( ProjectileResidualMassNumber > 1 );
        }
        if ( condition ) { 
          const G4int maxNumberOfInnerLoops = 1000;
          G4int innerLoopCounter = 0;
          do {  // Inner do while loop
            InnerSuccess = true;
            if ( isTargetToBeHandled ) {
              G4double Xcenter = 0.0;
              if ( interactionCase == 1 ) {
                common.PtNucleon = GaussianPt( AveragePt2, maxPtSquare );
                common.PtResidual = - common.PtNucleon;
                common.Mtarget =   std::sqrt( sqr( common.TNucleonMass )  + common.PtNucleon.mag2() ) 
                                 + std::sqrt( sqr( common.TResidualMass ) + common.PtResidual.mag2() );
                if ( common.SqrtS < common.Mprojectile + common.Mtarget ) {
                  InnerSuccess = false; 
                  continue;
                }
                Xcenter = std::sqrt( sqr( common.TNucleonMass ) + common.PtNucleon.mag2() ) 
                          / common.Mtarget;
              } else {
                Xcenter = std::sqrt( sqr( common.TNucleonMass ) + common.PtNucleonT.mag2() ) 
                          / common.Mtarget;
              }
              G4ThreeVector tmpX = GaussianPt( DcorT*DcorT, 1.0 );
              common.XminusNucleon = Xcenter + tmpX.x();
              if ( common.XminusNucleon <= 0.0  ||  common.XminusNucleon >= 1.0 ) {
                InnerSuccess = false; 
                continue;
              }
              common.XminusResidual = 1.0 - common.XminusNucleon;
            } else {  // Projectile to be handled
              G4ThreeVector tmpX = GaussianPt( DcorP*DcorP, 1.0 );
              G4double Xcenter = 0.0;
              if ( interactionCase == 2 ) {
                Xcenter = std::sqrt( sqr( common.TNucleonMass ) + common.PtNucleon.mag2() ) 
                          / common.Mprojectile;
              } else {
                Xcenter = std::sqrt( sqr( common.PNucleonMass ) + common.PtNucleonP.mag2() ) 
                          / common.Mprojectile;
              }
              common.XplusNucleon = Xcenter + tmpX.x();
              if ( common.XplusNucleon <= 0.0  ||  common.XplusNucleon >= 1.0 ) {
                InnerSuccess = false; 
                continue;
              }
              common.XplusResidual = 1.0 - common.XplusNucleon;
            }  // End-if on isTargetToBeHandled
          } while ( ( ! InnerSuccess ) &&                            // Inner do while loop
                      ++innerLoopCounter < maxNumberOfInnerLoops );  /* Loop checking, 10.08.2015, A.Ribon */ 
          if ( innerLoopCounter >= maxNumberOfInnerLoops ) {
            #ifdef debugAdjust
            G4cout << "BAD situation: forced exit of the inner while loop!" << G4endl;
            #endif
            return false;
          }
        } else {  // condition is false
          if ( isTargetToBeHandled ) {
            common.XminusNucleon  = 1.0;
            common.XminusResidual = 1.0;  // It must be 0, but in the calculation of Pz, E is problematic
          } else {  // Projectile to be handled
            common.XplusNucleon  = 1.0;
            common.XplusResidual = 1.0;   // It must be 0, but in the calculation of Pz, E is problematic
          } 
        }  // End-if on condition

      }  // End of for loop on iExecute

      if ( interactionCase == 1 ) {
        common.M2target =    ( sqr( common.TNucleonMass )  + common.PtNucleon.mag2() ) 
                             / common.XminusNucleon 
                          +  ( sqr( common.TResidualMass ) + common.PtResidual.mag2() ) 
                             / common.XminusResidual;
        loopCondition = ( common.SqrtS < common.Mprojectile + std::sqrt( common.M2target ) );
      } else if ( interactionCase == 2 ) {
        #ifdef debugAdjust
        G4cout << "TNucleonMass PtNucleon XplusNucleon " << common.TNucleonMass << " " 
               << common.PtNucleon << " " << common.XplusNucleon << G4endl
               << "TResidualMass PtResidual XplusResidual " << common.TResidualMass << " " 
               << common.PtResidual << "  " << common.XplusResidual << G4endl;
        #endif
        common.M2projectile =    ( sqr( common.TNucleonMass )  + common.PtNucleon.mag2() ) 
                                 / common.XplusNucleon 
                              +  ( sqr( common.TResidualMass ) + common.PtResidual.mag2() ) 
                                 / common.XplusResidual;
        #ifdef debugAdjust
        G4cout << "SqrtS < Mtarget + std::sqrt(M2projectile) " << common.SqrtS << "  " 
               << common.Mtarget << " " << std::sqrt( common.M2projectile ) << " "
               << common.Mtarget + std::sqrt( common.M2projectile ) << G4endl;
        #endif
        loopCondition = ( common.SqrtS < common.Mtarget + std::sqrt( common.M2projectile ) );
      } else if ( interactionCase == 3 ) {
        #ifdef debugAdjust
        G4cout << "PtNucleonP " << common.PtNucleonP << " " << common.PtResidualP << G4endl
               << "XplusNucleon XplusResidual " << common.XplusNucleon 
               << " " << common.XplusResidual << G4endl
               << "PtNucleonT " << common.PtNucleonT << " " << common.PtResidualT << G4endl
               << "XminusNucleon XminusResidual " << common.XminusNucleon 
               << " " << common.XminusResidual << G4endl;
        #endif
        common.M2projectile =   ( sqr( common.PNucleonMass ) + common.PtNucleonP.mag2() ) 
                                / common.XplusNucleon 
                              + ( sqr( common.PResidualMass) + common.PtResidualP.mag2() ) 
                                / common.XplusResidual;
        common.M2target =    ( sqr( common.TNucleonMass )  + common.PtNucleonT.mag2() ) 
                             / common.XminusNucleon 
                          +  ( sqr( common.TResidualMass ) + common.PtResidualT.mag2() ) 
                             / common.XminusResidual;
        loopCondition = ( common.SqrtS < (   std::sqrt( common.M2projectile ) 
					   + std::sqrt( common.M2target ) ) );
      }  // End-if on interactionCase

    } while ( loopCondition &&                       // Intermediate do while loop
              ++NumberOfTries < maxNumberOfTries );  /* Loop checking, 10.08.2015, A.Ribon */
    if ( NumberOfTries >= maxNumberOfTries ) { 
      #ifdef debugAdjust
      G4cout << "BAD situation: forced exit of the intermediate while loop!" << G4endl;
      #endif
      return false;
    }

    // kinematics
    G4double Yprojectile = 0.0, YprojectileNucleon = 0.0, Ytarget = 0.0, YtargetNucleon = 0.0;
    G4double DecayMomentum2 = sqr( common.S ) + sqr( common.M2projectile ) + sqr( common.M2target )
                              - 2.0 * ( common.S * ( common.M2projectile + common.M2target )
                                        + common.M2projectile * common.M2target );
    if ( interactionCase == 1 ) {
      common.WminusTarget = (   common.S - common.M2projectile + common.M2target 
                              + std::sqrt( DecayMomentum2 ) ) / 2.0 / common.SqrtS;
      common.WplusProjectile = common.SqrtS - common.M2target / common.WminusTarget;
      common.Pzprojectile =   common.WplusProjectile / 2.0 
                            - common.M2projectile / 2.0 / common.WplusProjectile;
      common.Eprojectile =    common.WplusProjectile / 2.0 
                            + common.M2projectile / 2.0 / common.WplusProjectile;
      Yprojectile  = 0.5 * G4Log(   ( common.Eprojectile + common.Pzprojectile )
                                  / ( common.Eprojectile - common.Pzprojectile ) );
      #ifdef debugAdjust
      G4cout << "DecayMomentum2 " << DecayMomentum2 << G4endl
             << "WminusTarget WplusProjectile " << common.WminusTarget 
             << " " << common.WplusProjectile << G4endl 
             << "Yprojectile " << Yprojectile << G4endl;
      #endif
      common.Mt2targetNucleon = sqr( common.TNucleonMass ) + common.PtNucleon.mag2();
      common.PztargetNucleon = - common.WminusTarget * common.XminusNucleon / 2.0
                               + common.Mt2targetNucleon 
                                 / ( 2.0 * common.WminusTarget * common.XminusNucleon );
      common.EtargetNucleon =   common.WminusTarget * common.XminusNucleon / 2.0
                              + common.Mt2targetNucleon
                                / ( 2.0 * common.WminusTarget * common.XminusNucleon );
      YtargetNucleon = 0.5 * G4Log(   ( common.EtargetNucleon + common.PztargetNucleon )
                                    / ( common.EtargetNucleon - common.PztargetNucleon ) );
      #ifdef debugAdjust
      G4cout << "YtN Ytr YtN-Ytr " << " " << YtargetNucleon << " " << common.YtargetNucleus 
             << " " << YtargetNucleon - common.YtargetNucleus << G4endl
             << "YtN Ypr YtN-Ypr " << " " << YtargetNucleon << " " << Yprojectile
             << " " << YtargetNucleon - Yprojectile << G4endl;
      #endif
      if ( std::abs( YtargetNucleon - common.YtargetNucleus ) > 2  ||
           Yprojectile < YtargetNucleon ) {
        OuterSuccess = false; 
        continue;
      }
    } else if ( interactionCase == 2 ) {
      common.WplusProjectile = (   common.S + common.M2projectile - common.M2target 
                                 + std::sqrt( DecayMomentum2 ) ) / 2.0 / common.SqrtS;
      common.WminusTarget = common.SqrtS - common.M2projectile / common.WplusProjectile;
      common.Pztarget = - common.WminusTarget / 2.0 + common.M2target / 2.0 / common.WminusTarget;
      common.Etarget =    common.WminusTarget / 2.0 + common.M2target / 2.0 / common.WminusTarget;
      Ytarget = 0.5 * G4Log(   ( common.Etarget + common.Pztarget )
                             / ( common.Etarget - common.Pztarget ) );
      #ifdef debugAdjust
      G4cout << "DecayMomentum2 " << DecayMomentum2 << G4endl
             << "WminusTarget WplusProjectile " << common.WminusTarget 
             << " " << common.WplusProjectile << G4endl 
             << "Ytarget " << Ytarget << G4endl;
      #endif
      common.Mt2projectileNucleon = sqr( common.TNucleonMass ) + common.PtNucleon.mag2();
      common.PzprojectileNucleon =   common.WplusProjectile * common.XplusNucleon / 2.0
                                   - common.Mt2projectileNucleon
                                     / ( 2.0 * common.WplusProjectile * common.XplusNucleon );
      common.EprojectileNucleon =    common.WplusProjectile * common.XplusNucleon / 2.0 
                                   + common.Mt2projectileNucleon
                                     / ( 2.0 * common.WplusProjectile * common.XplusNucleon );
      YprojectileNucleon = 0.5 * G4Log(   ( common.EprojectileNucleon + common.PzprojectileNucleon )
                                        / ( common.EprojectileNucleon - common.PzprojectileNucleon) );
      #ifdef debugAdjust
      G4cout << "YpN Ypr YpN-Ypr " << " " << YprojectileNucleon << " " << common.YprojectileNucleus
             << " " << YprojectileNucleon - common.YprojectileNucleus << G4endl
             << "YpN Ytr YpN-Ytr " << " " << YprojectileNucleon << " " << Ytarget
             << " " << YprojectileNucleon - Ytarget << G4endl;
      #endif
      if ( std::abs( YprojectileNucleon - common.YprojectileNucleus ) > 2  ||
           Ytarget > YprojectileNucleon ) {
        OuterSuccess = false; 
        continue;
      }
    } else if ( interactionCase == 3 ) {
      common.WplusProjectile = (   common.S + common.M2projectile - common.M2target 
                                 + std::sqrt( DecayMomentum2 ) ) / 2.0 / common.SqrtS;
      common.WminusTarget = common.SqrtS - common.M2projectile / common.WplusProjectile;
      common.Mt2projectileNucleon = sqr( common.PNucleonMass ) + common.PtNucleonP.mag2();
      common.PzprojectileNucleon =   common.WplusProjectile * common.XplusNucleon / 2.0
                                   - common.Mt2projectileNucleon
                                     / ( 2.0 * common.WplusProjectile * common.XplusNucleon );
      common.EprojectileNucleon =    common.WplusProjectile * common.XplusNucleon / 2.0
                                   + common.Mt2projectileNucleon
                                     / ( 2.0 * common.WplusProjectile * common.XplusNucleon );
      YprojectileNucleon = 0.5 * G4Log(   ( common.EprojectileNucleon + common.PzprojectileNucleon )
                                        / ( common.EprojectileNucleon - common.PzprojectileNucleon ) );
      common.Mt2targetNucleon = sqr( common.TNucleonMass ) + common.PtNucleonT.mag2();
      common.PztargetNucleon = - common.WminusTarget * common.XminusNucleon / 2.0
                               + common.Mt2targetNucleon
                                 / ( 2.0 * common.WminusTarget * common.XminusNucleon );
      common.EtargetNucleon =    common.WminusTarget * common.XminusNucleon / 2.0
                               + common.Mt2targetNucleon
                                 / ( 2.0 * common.WminusTarget * common.XminusNucleon );
      YtargetNucleon = 0.5 * G4Log(   ( common.EtargetNucleon + common.PztargetNucleon )
                                    / ( common.EtargetNucleon - common.PztargetNucleon ) ); 
      if ( std::abs( YtargetNucleon - common.YtargetNucleus ) > 2          ||
           std::abs( YprojectileNucleon - common.YprojectileNucleus ) > 2  ||
           YprojectileNucleon < YtargetNucleon ) {
        OuterSuccess = false;
        continue;
      }
    }  // End-if on interactionCase

  } while ( ( ! OuterSuccess ) &&                // Outmost do while loop
            ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */
  if ( loopCounter >= maxNumberOfLoops ) {
    #ifdef debugAdjust
    G4cout << "BAD situation: forced exit of the while loop!" << G4endl;
    #endif
    return false;
  }

  return true;
}

//-------------------------------------------------------------------

void G4FTFModel::AdjustNucleonsAlgorithm_afterSampling( G4int interactionCase,
                                                        G4VSplitableHadron* SelectedAntiBaryon,
                                                        G4VSplitableHadron* SelectedTargetNucleon,
                                                        G4FTFModel::CommonVariables& common ) {
  // Third of the three utility methods used only by AdjustNucleons: do the final kinematics
  // and transform back.

  // New projectile
  if ( interactionCase == 1 ) {
    common.Pprojectile.setPz( common.Pzprojectile );  
    common.Pprojectile.setE( common.Eprojectile );
  } else if ( interactionCase == 2 ) {
    common.Pprojectile.setPx( common.PtNucleon.x() ); 
    common.Pprojectile.setPy( common.PtNucleon.y() );
    common.Pprojectile.setPz( common.PzprojectileNucleon );
    common.Pprojectile.setE( common.EprojectileNucleon ); 
  } else if ( interactionCase == 3 ) {
    common.Pprojectile.setPx( common.PtNucleonP.x() );
    common.Pprojectile.setPy( common.PtNucleonP.y() );
    common.Pprojectile.setPz( common.PzprojectileNucleon );
    common.Pprojectile.setE( common.EprojectileNucleon );
  }
  #ifdef debugAdjust
  G4cout << "Proj after in CMS " << common.Pprojectile << G4endl;
  #endif
  common.Pprojectile.transform( common.toLab );
  SelectedAntiBaryon->Set4Momentum( common.Pprojectile );
  #ifdef debugAdjust
  G4cout << "Proj after in Lab " << common.Pprojectile << G4endl;
  #endif

  // New target nucleon
  if ( interactionCase == 1 ) {
    common.Ptarget.setPx( common.PtNucleon.x() );
    common.Ptarget.setPy( common.PtNucleon.y() );
    common.Ptarget.setPz( common.PztargetNucleon );
    common.Ptarget.setE( common.EtargetNucleon ); 
  } else if ( interactionCase == 2 ) {
    common.Ptarget.setPz( common.Pztarget ); 
    common.Ptarget.setE( common.Etarget );
  } else if ( interactionCase == 3 ) {
    common.Ptarget.setPx( common.PtNucleonT.x() );
    common.Ptarget.setPy( common.PtNucleonT.y() );
    common.Ptarget.setPz( common.PztargetNucleon );
    common.Ptarget.setE( common.EtargetNucleon );
  }
  #ifdef debugAdjust
  G4cout << "Targ after in CMS " << common.Ptarget << G4endl;
  #endif
  common.Ptarget.transform( common.toLab );
  SelectedTargetNucleon->Set4Momentum( common.Ptarget );
  #ifdef debugAdjust
  G4cout << "Targ after in Lab " << common.Ptarget << G4endl;
  #endif

  // New target residual
  if ( interactionCase == 1  ||  interactionCase == 3 ) {
    TargetResidualMassNumber       = common.TResidualMassNumber;
    TargetResidualCharge           = common.TResidualCharge;
    TargetResidualExcitationEnergy = common.TResidualExcitationEnergy;
    #ifdef debugAdjust
    G4cout << "TargetResidualMassNumber TargetResidualCharge TargetResidualExcitationEnergy "
           << TargetResidualMassNumber << " " << TargetResidualCharge << " "
           << TargetResidualExcitationEnergy << G4endl;
    #endif
    if ( TargetResidualMassNumber != 0 ) {
      G4double Mt2 = 0.0;
      if ( interactionCase == 1 ) {
        Mt2 = sqr( common.TResidualMass ) + common.PtResidual.mag2();
        TargetResidual4Momentum.setPx( common.PtResidual.x() );
        TargetResidual4Momentum.setPy( common.PtResidual.y() );
      } else {  // interactionCase == 3
        Mt2 = sqr( common.TResidualMass ) + common.PtResidualT.mag2();
        TargetResidual4Momentum.setPx( common.PtResidualT.x() );
        TargetResidual4Momentum.setPy( common.PtResidualT.y() );
      }
      G4double Pz = - common.WminusTarget * common.XminusResidual / 2.0
                    + Mt2 / ( 2.0 * common.WminusTarget * common.XminusResidual );
      G4double E =    common.WminusTarget * common.XminusResidual / 2.0 
                    + Mt2 / ( 2.0 * common.WminusTarget * common.XminusResidual );
      TargetResidual4Momentum.setPz( Pz );
      TargetResidual4Momentum.setE( E ) ;
      TargetResidual4Momentum.transform( common.toLab );
    } else {
      TargetResidual4Momentum = G4LorentzVector( 0.0, 0.0, 0.0, 0.0 );
    }
    #ifdef debugAdjust
    G4cout << "Tr N R " << common.Ptarget << G4endl << "       " << TargetResidual4Momentum << G4endl;
    #endif
  }

  // New projectile residual
  if ( interactionCase == 2  ||  interactionCase == 3 ) {
    if ( interactionCase == 2 ) {
      ProjectileResidualMassNumber       = common.TResidualMassNumber;
      ProjectileResidualCharge           = common.TResidualCharge;
      ProjectileResidualExcitationEnergy = common.TResidualExcitationEnergy;
      ProjectileResidualLambdaNumber     = common.PResidualLambdaNumber;
    } else {  // interactionCase == 3
      ProjectileResidualMassNumber       = common.PResidualMassNumber;
      ProjectileResidualCharge           = common.PResidualCharge;
      ProjectileResidualExcitationEnergy = common.PResidualExcitationEnergy;
      ProjectileResidualLambdaNumber     = common.PResidualLambdaNumber;
    }
    #ifdef debugAdjust
    G4cout << "ProjectileResidualMassNumber ProjectileResidualCharge  Lambdas ProjectileResidualExcitationEnergy "
           << ProjectileResidualMassNumber << " " << ProjectileResidualCharge << " "
           << ProjectileResidualLambdaNumber << " "
           << ProjectileResidualExcitationEnergy << G4endl;
    #endif
    if ( ProjectileResidualMassNumber != 0 ) {
      G4double Mt2 = 0.0;
      if ( interactionCase == 2 ) {
        Mt2 = sqr( common.TResidualMass ) + common.PtResidual.mag2();
        ProjectileResidual4Momentum.setPx( common.PtResidual.x() );
        ProjectileResidual4Momentum.setPy( common.PtResidual.y() );
      } else {  // interactionCase == 3
        Mt2 = sqr( common.PResidualMass ) + common.PtResidualP.mag2();
        ProjectileResidual4Momentum.setPx( common.PtResidualP.x() );
        ProjectileResidual4Momentum.setPy( common.PtResidualP.y() );
      }
      G4double Pz =   common.WplusProjectile * common.XplusResidual / 2.0
                    - Mt2 / ( 2.0 * common.WplusProjectile * common.XplusResidual );
      G4double E  =   common.WplusProjectile * common.XplusResidual / 2.0 
                    + Mt2 / ( 2.0 * common.WplusProjectile * common.XplusResidual );
      ProjectileResidual4Momentum.setPz( Pz );
      ProjectileResidual4Momentum.setE( E );
      ProjectileResidual4Momentum.transform( common.toLab );
    } else {
      ProjectileResidual4Momentum = G4LorentzVector( 0.0, 0.0, 0.0, 0.0 );
    }
    #ifdef debugAdjust
    G4cout << "Pr N R " << common.Pprojectile << G4endl 
           << "       " << ProjectileResidual4Momentum << G4endl;
    #endif
  }

}


//============================================================================

void G4FTFModel::BuildStrings( G4ExcitedStringVector* strings ) {
  // Loop over all collisions; find all primaries, and all targets 
  // (targets may be duplicate in the List (to unique G4VSplitableHadrons) ).

  G4ExcitedString* FirstString( 0 );     // If there will be a kink,
  G4ExcitedString* SecondString( 0 );    // two strings will be produced.

  if ( ! GetProjectileNucleus() ) {

    std::vector< G4VSplitableHadron* > primaries;
    theParticipants.StartLoop();
    while ( theParticipants.Next() ) {  /* Loop checking, 10.08.2015, A.Ribon */
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
      //G4cout << "primaries[ ahadron ] " << primaries[ ahadron ] << G4endl;
      //if ( primaries[ ahadron ]->GetStatus() <= 1 ) isProjectile = true;
      FirstString = 0; SecondString = 0;
      if ( primaries[ahadron]->GetStatus() == 0 ) {
       theExcitation->CreateStrings( primaries[ ahadron ], isProjectile, 
                                     FirstString, SecondString, theParameters );
       NumberOfProjectileSpectatorNucleons--;
      } else if ( primaries[ahadron]->GetStatus() == 1  
               && primaries[ahadron]->GetSoftCollisionCount() != 0 ) {
       theExcitation->CreateStrings( primaries[ ahadron ], isProjectile, 
                                     FirstString, SecondString, theParameters );
       NumberOfProjectileSpectatorNucleons--;
      } else if ( primaries[ahadron]->GetStatus() == 1  
               && primaries[ahadron]->GetSoftCollisionCount() == 0 ) {
       G4LorentzVector ParticleMomentum=primaries[ahadron]->Get4Momentum();
       G4KineticTrack* aTrack = new G4KineticTrack( primaries[ahadron]->GetDefinition(),
                                                    primaries[ahadron]->GetTimeOfCreation(),
                                                    primaries[ahadron]->GetPosition(),
                                                    ParticleMomentum );
       FirstString = new G4ExcitedString( aTrack );
      } else if (primaries[ahadron]->GetStatus() == 2) {
       G4LorentzVector ParticleMomentum=primaries[ahadron]->Get4Momentum();
       G4KineticTrack* aTrack = new G4KineticTrack( primaries[ahadron]->GetDefinition(),
                                                    primaries[ahadron]->GetTimeOfCreation(),
                                                    primaries[ahadron]->GetPosition(),
                                                    ParticleMomentum );
       FirstString = new G4ExcitedString( aTrack );
       NumberOfProjectileSpectatorNucleons--;
      } else {
        G4cout << "Something wrong in FTF Model Build String" << G4endl;
      }

      if ( FirstString  != 0 ) strings->push_back( FirstString );
      if ( SecondString != 0 ) strings->push_back( SecondString );

      #ifdef debugBuildString
      G4cout << "FirstString & SecondString? " << FirstString << " " << SecondString << G4endl;
      if ( FirstString->IsExcited() ) {
        G4cout << "Quarks on the FirstString ends " << FirstString->GetRightParton()->GetPDGcode()
               << " " << FirstString->GetLeftParton()->GetPDGcode() << G4endl;
      } else {
        G4cout << "Kinetic track is stored" << G4endl;
      }
      #endif

    }

    #ifdef debugBuildString
    if ( FirstString->IsExcited() ) {
      G4cout << "Check 1 string " << strings->operator[](0)->GetRightParton()->GetPDGcode() 
             << " " << strings->operator[](0)->GetLeftParton()->GetPDGcode() << G4endl << G4endl;
    }
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
                       ->GetSoftCollisionCount()<<G4endl;
      #endif

      G4VSplitableHadron* aProjectile = 
          TheInvolvedNucleonsOfProjectile[ ahadron ]->GetSplitableHadron();

      #ifdef debugBuildString
      G4cout << G4endl << "ahadron aProjectile Status " << ahadron << " " << aProjectile
             << " " << aProjectile->GetStatus() << G4endl;
      #endif

      FirstString = 0; SecondString = 0;
      if ( aProjectile->GetStatus() == 0 ) {  // A nucleon took part in non-diffractive interaction

        #ifdef debugBuildString
        G4cout << "Case1 aProjectile->GetStatus() == 0 " << G4endl;
        #endif

        theExcitation->CreateStrings( 
                           TheInvolvedNucleonsOfProjectile[ ahadron ]->GetSplitableHadron(),
                           isProjectile, FirstString, SecondString, theParameters );
        NumberOfProjectileSpectatorNucleons--;
      } else if ( aProjectile->GetStatus() == 1 && aProjectile->GetSoftCollisionCount() != 0 ) { 
        // Nucleon took part in diffractive interaction

        #ifdef debugBuildString
        G4cout << "Case2 aProjectile->GetStatus() !=0 St==1 SoftCol!=0" << G4endl;
        #endif

        theExcitation->CreateStrings( 
                           TheInvolvedNucleonsOfProjectile[ ahadron ]->GetSplitableHadron(),
                           isProjectile, FirstString, SecondString, theParameters );
        NumberOfProjectileSpectatorNucleons--;
      } else if ( aProjectile->GetStatus() == 1  &&  aProjectile->GetSoftCollisionCount() == 0  &&
                  HighEnergyInter ) {
        // Nucleon was considered as a paricipant of an interaction,
        // but the interaction was skipped due to annihilation.
        // It is now considered as an involved nucleon at high energies.

        #ifdef debugBuildString
        G4cout << "Case3 aProjectile->GetStatus() !=0 St==1 SoftCol==0" << G4endl;
        #endif

        G4LorentzVector ParticleMomentum = aProjectile->Get4Momentum();
        G4KineticTrack* aTrack = new G4KineticTrack( aProjectile->GetDefinition(),
                                                     aProjectile->GetTimeOfCreation(),
                                                     aProjectile->GetPosition(),
                                                     ParticleMomentum );
        FirstString = new G4ExcitedString( aTrack );

        #ifdef debugBuildString
        G4cout << " Strings are built for nucleon marked for an interaction, but"
               << " the interaction was skipped." << G4endl;
        #endif

      } else if ( aProjectile->GetStatus() == 2  ||  aProjectile->GetStatus() == 3 ) {
        // Nucleon which was involved in the Reggeon cascading

        #ifdef debugBuildString
        G4cout << "Case4 aProjectile->GetStatus() !=0 St==2 " << G4endl;
        #endif

        G4LorentzVector ParticleMomentum = aProjectile->Get4Momentum();
        G4KineticTrack* aTrack = new G4KineticTrack( aProjectile->GetDefinition(),
                                                     aProjectile->GetTimeOfCreation(),
                                                     aProjectile->GetPosition(),
                                                     ParticleMomentum );
        FirstString = new G4ExcitedString( aTrack );

        #ifdef debugBuildString
        G4cout << " Strings are build for involved nucleon." << G4endl;
        #endif

        if ( aProjectile->GetStatus() == 2 ) NumberOfProjectileSpectatorNucleons--;
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

      if ( FirstString  != 0 ) strings->push_back( FirstString );
      if ( SecondString != 0 ) strings->push_back( SecondString );
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
           << aNucleon->GetStatus() << " " << aNucleon->GetSoftCollisionCount()<<G4endl;;
    #endif

    FirstString = 0 ; SecondString = 0;

    if ( aNucleon->GetStatus() == 0 ) { // A nucleon took part in non-diffractive interaction
      theExcitation->CreateStrings( aNucleon, isProjectile, 
                                    FirstString, SecondString, theParameters );
      NumberOfTargetSpectatorNucleons--;

      #ifdef debugBuildString
      G4cout << " 1 case A string is build" << G4endl;
      #endif

    } else if ( aNucleon->GetStatus() == 1  &&  aNucleon->GetSoftCollisionCount() != 0 ) {
      // A nucleon took part in diffractive interaction
      theExcitation->CreateStrings( aNucleon, isProjectile,
                                    FirstString, SecondString, theParameters );

      #ifdef debugBuildString
      G4cout << " 2 case A string is build, nucleon was excited." << G4endl;
      #endif

      NumberOfTargetSpectatorNucleons--;

    } else if ( aNucleon->GetStatus() == 1  &&  aNucleon->GetSoftCollisionCount() == 0  &&
                HighEnergyInter ) {
      // A nucleon was considered as a participant but due to annihilation
      // its interactions were skipped. It will be considered as involved one
      // at high energies.

      G4LorentzVector ParticleMomentum = aNucleon->Get4Momentum();
      G4KineticTrack* aTrack = new G4KineticTrack( aNucleon->GetDefinition(),
                                                   aNucleon->GetTimeOfCreation(),
                                                   aNucleon->GetPosition(),
                                                   ParticleMomentum );

      FirstString = new G4ExcitedString( aTrack );

      #ifdef debugBuildString
      G4cout << "3 case A string is build" << G4endl;
      #endif

    } else if ( aNucleon->GetStatus() == 1  &&  aNucleon->GetSoftCollisionCount() == 0  &&
                ! HighEnergyInter ) {
      // A nucleon was considered as a participant but due to annihilation
      // its interactions were skipped. It will be returned to nucleus
      // at low energies energies.
      aNucleon->SetStatus( 5 );  // 4->5
      // ??? delete aNucleon;

      #ifdef debugBuildString
      G4cout << "4 case A string is not build" << G4endl;
      #endif

    } else if ( aNucleon->GetStatus() == 2  ||   // A nucleon took part in quark exchange
                aNucleon->GetStatus() == 3  ) {  // A nucleon was involved in Reggeon cascading
      G4LorentzVector ParticleMomentum = aNucleon->Get4Momentum();
      G4KineticTrack* aTrack = new G4KineticTrack( aNucleon->GetDefinition(),
                                                   aNucleon->GetTimeOfCreation(),
                                                   aNucleon->GetPosition(),
                                                   ParticleMomentum );
      FirstString = new G4ExcitedString( aTrack );

      #ifdef debugBuildString
      G4cout << "5 case A string is build" << G4endl;
      #endif

      if ( aNucleon->GetStatus() == 2 ) NumberOfTargetSpectatorNucleons--;

    } else {

      #ifdef debugBuildString
      G4cout << "6 case No string" << G4endl;
      #endif

    }

    if ( FirstString  != 0 ) strings->push_back( FirstString );
    if ( SecondString != 0 ) strings->push_back( SecondString );

  }

  #ifdef debugBuildString
  G4cout << G4endl << "theAdditionalString.size() " << theAdditionalString.size() 
         << G4endl << G4endl;
  #endif

  isProjectile = true;
  if ( theAdditionalString.size() != 0 ) {
    for ( unsigned int  ahadron = 0; ahadron < theAdditionalString.size(); ahadron++ ) {
      //if ( theAdditionalString[ ahadron ]->GetStatus() <= 1 ) isProjectile = true; 
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

  return;
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

    for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; ++i ) {
      G4Nucleon* aNucleon = TheInvolvedNucleonsOfTarget[i];

      #ifdef debugFTFmodel
      G4VSplitableHadron* targetSplitable = aNucleon->GetSplitableHadron();
      G4cout << i << " Hit? " << aNucleon->AreYouHit() << " pointer " << targetSplitable << G4endl;
      if ( targetSplitable ) G4cout << i << "Status " << targetSplitable->GetStatus() << G4endl;
      #endif

      G4LorentzVector tmp = -DeltaPResidualNucleus;
      aNucleon->SetMomentum( tmp );
      aNucleon->SetBindingEnergy( DeltaExcitationE );
    }

    if ( TargetResidualMassNumber != 0 ) {
      G4ThreeVector bstToCM = TargetResidual4Momentum.findBoostToCM();

      G4V3DNucleus* theTargetNucleus = GetTargetNucleus();
      G4LorentzVector residualMomentum( 0.0, 0.0, 0.0, 0.0 );
      G4Nucleon* aNucleon = 0;
      theTargetNucleus->StartLoop();
      while ( ( aNucleon = theTargetNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
        if ( ! aNucleon->AreYouHit() ) { 
          G4LorentzVector tmp = aNucleon->Get4Momentum(); tmp.boost( bstToCM );
          aNucleon->SetMomentum( tmp );
          residualMomentum += tmp;
        }
      }

      residualMomentum /= TargetResidualMassNumber;

      G4double Mass = TargetResidual4Momentum.mag();
      G4double SumMasses = 0.0;
  
      aNucleon = 0;
      theTargetNucleus->StartLoop();
      while ( ( aNucleon = theTargetNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
        if ( ! aNucleon->AreYouHit() ) { 
          G4LorentzVector tmp = aNucleon->Get4Momentum() - residualMomentum;
          G4double E = std::sqrt( tmp.vect().mag2() +
                                  sqr( aNucleon->GetDefinition()->GetPDGMass() - aNucleon->GetBindingEnergy() ) );
          tmp.setE( E );  aNucleon->SetMomentum( tmp );
          SumMasses += E;
        }
      }

      G4double Chigh = Mass / SumMasses; G4double Clow = 0.0; G4double C;
      const G4int maxNumberOfLoops = 1000;
      G4int loopCounter = 0;
      do {
        C = ( Chigh + Clow ) / 2.0;
        SumMasses = 0.0;
        aNucleon = 0;
        theTargetNucleus->StartLoop();
        while ( ( aNucleon = theTargetNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
          if ( ! aNucleon->AreYouHit() ) { 
            G4LorentzVector tmp = aNucleon->Get4Momentum();
            G4double E = std::sqrt( tmp.vect().mag2()*sqr(C) +
                                    sqr( aNucleon->GetDefinition()->GetPDGMass() - aNucleon->GetBindingEnergy() ) );
            SumMasses += E;
          }
        }

        if ( SumMasses > Mass ) Chigh = C;
        else                    Clow  = C;

      } while ( Chigh - Clow > 0.01  &&
                ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */
      if ( loopCounter >= maxNumberOfLoops ) {
        #ifdef debugFTFmodel
        G4cout << "BAD situation: forced exit of the first while loop in G4FTFModel::GetResidual" << G4endl
               << "\t return immediately from the method!" << G4endl;
        #endif
        return;
      }

      aNucleon = 0;
      theTargetNucleus->StartLoop();
      while ( ( aNucleon = theTargetNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
        if ( !aNucleon->AreYouHit() ) { 
          G4LorentzVector tmp = aNucleon->Get4Momentum()*C;
          G4double E = std::sqrt( tmp.vect().mag2()+
                                  sqr( aNucleon->GetDefinition()->GetPDGMass()-aNucleon->GetBindingEnergy() ) );
          tmp.setE( E ); tmp.boost( -bstToCM );  
          aNucleon->SetMomentum( tmp );     
        }
      }
    }

    if ( ! GetProjectileNucleus() ) return;  // The projectile is a hadron

    #ifdef debugFTFmodel
    G4cout << "NumberOfInvolvedNucleonsOfProjectile " << NumberOfInvolvedNucleonsOfProjectile
           << G4endl << "ProjectileResidualExcitationEnergy ProjectileResidual4Momentum "
           << ProjectileResidualExcitationEnergy << "  " << ProjectileResidual4Momentum << G4endl;
    #endif

    DeltaExcitationE = ProjectileResidualExcitationEnergy /
                       G4double( NumberOfInvolvedNucleonsOfProjectile );
    DeltaPResidualNucleus = ProjectileResidual4Momentum /
                            G4double( NumberOfInvolvedNucleonsOfProjectile );

    for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfProjectile; ++i ) {
      G4Nucleon* aNucleon = TheInvolvedNucleonsOfProjectile[i];

      #ifdef debugFTFmodel
      G4VSplitableHadron* projSplitable = aNucleon->GetSplitableHadron();
      G4cout << i << " Hit? " << aNucleon->AreYouHit() << " pointer " << projSplitable << G4endl;
      if ( projSplitable ) G4cout << i << "Status " << projSplitable->GetStatus() << G4endl;
      #endif

      G4LorentzVector tmp = -DeltaPResidualNucleus;
      aNucleon->SetMomentum( tmp );
      aNucleon->SetBindingEnergy( DeltaExcitationE );
    }

    if ( ProjectileResidualMassNumber != 0 ) {
      G4ThreeVector bstToCM = ProjectileResidual4Momentum.findBoostToCM();

      G4V3DNucleus* theProjectileNucleus = GetProjectileNucleus();
      G4LorentzVector residualMomentum( 0.0, 0.0, 0.0, 0.0);
      G4Nucleon* aNucleon = 0;
      theProjectileNucleus->StartLoop();
      while ( ( aNucleon = theProjectileNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
        if ( ! aNucleon->AreYouHit() ) { 
          G4LorentzVector tmp = aNucleon->Get4Momentum(); tmp.boost( bstToCM );
          aNucleon->SetMomentum( tmp );
          residualMomentum += tmp;
        }
      }

      residualMomentum /= ProjectileResidualMassNumber;

      G4double Mass = ProjectileResidual4Momentum.mag();
      G4double SumMasses= 0.0;
  
      aNucleon = 0;
      theProjectileNucleus->StartLoop();
      while ( ( aNucleon = theProjectileNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
        if ( ! aNucleon->AreYouHit() ) { 
          G4LorentzVector tmp = aNucleon->Get4Momentum() - residualMomentum;
          G4double E=std::sqrt( tmp.vect().mag2() +
                                sqr(aNucleon->GetDefinition()->GetPDGMass()-aNucleon->GetBindingEnergy() ) );
          tmp.setE( E );  aNucleon->SetMomentum( tmp );
          SumMasses += E;
        }
      }

      G4double Chigh = Mass / SumMasses; G4double Clow = 0.0; G4double C;
      const G4int maxNumberOfLoops = 1000;
      G4int loopCounter = 0;
      do {
        C = ( Chigh + Clow ) / 2.0;

        SumMasses = 0.0;
        aNucleon = 0;
        theProjectileNucleus->StartLoop();
        while ( ( aNucleon = theProjectileNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
          if ( ! aNucleon->AreYouHit() ) { 
            G4LorentzVector tmp = aNucleon->Get4Momentum();
            G4double E = std::sqrt( tmp.vect().mag2()*sqr(C) +
                                    sqr( aNucleon->GetDefinition()->GetPDGMass() - aNucleon->GetBindingEnergy() ) );
            SumMasses += E;
          }
        }

        if ( SumMasses > Mass) Chigh = C;
        else                   Clow  = C;

      } while ( Chigh - Clow > 0.01  &&
                ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */
      if ( loopCounter >= maxNumberOfLoops ) {
        #ifdef debugFTFmodel
        G4cout << "BAD situation: forced exit of the second while loop in G4FTFModel::GetResidual" << G4endl
               << "\t return immediately from the method!" << G4endl;
        #endif
        return;
      }

      aNucleon = 0;
      theProjectileNucleus->StartLoop();
      while ( ( aNucleon = theProjectileNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
        if ( ! aNucleon->AreYouHit() ) { 
          G4LorentzVector tmp = aNucleon->Get4Momentum()*C;
          G4double E = std::sqrt( tmp.vect().mag2() +
                                  sqr( aNucleon->GetDefinition()->GetPDGMass() - aNucleon->GetBindingEnergy() ) );
          tmp.setE( E ); tmp.boost( -bstToCM );  
          aNucleon->SetMomentum( tmp ); 
        }
      }
    }   // End of if ( ProjectileResidualMassNumber != 0 )
  
    #ifdef debugFTFmodel
    G4cout << "End projectile" << G4endl;
    #endif
   
  } else {  // Related to the condition: if ( HighEnergyInter )

    #ifdef debugFTFmodel
    G4cout << "Low energy interaction: Target nucleus --------------" << G4endl
           << "Tr ResidualMassNumber Tr ResidualCharge Tr ResidualExcitationEnergy  "
           << TargetResidualMassNumber << " " << TargetResidualCharge << " "
           << TargetResidualExcitationEnergy << G4endl;
    #endif

    G4int NumberOfTargetParticipant( 0 );
    for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; ++i ) {
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

    for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfTarget; ++i ) {
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

    if ( ! GetProjectileNucleus() ) return;  // The projectile is a hadron

    #ifdef debugFTFmodel
    G4cout << "Low energy interaction: Projectile nucleus --------------" << G4endl
           << "Pr ResidualMassNumber Pr ResidualCharge Pr ResidualExcitationEnergy "
           << ProjectileResidualMassNumber << " " << ProjectileResidualCharge << " "
           << ProjectileResidualExcitationEnergy << G4endl;
    #endif

    G4int NumberOfProjectileParticipant( 0 );
    for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfProjectile; ++i ) {
      G4Nucleon* aNucleon = TheInvolvedNucleonsOfProjectile[i];
      G4VSplitableHadron* projectileSplitable = aNucleon->GetSplitableHadron();
      if ( projectileSplitable->GetSoftCollisionCount() != 0 ) NumberOfProjectileParticipant++;
    }
 
    #ifdef debugFTFmodel
    G4cout << "NumberOfProjectileParticipant" << G4endl;
    #endif

    DeltaExcitationE = 0.0;
    DeltaPResidualNucleus = G4LorentzVector( 0.0, 0.0, 0.0, 0.0 );

    if ( NumberOfProjectileParticipant != 0 ) {
      DeltaExcitationE = ProjectileResidualExcitationEnergy / G4double( NumberOfProjectileParticipant );
      DeltaPResidualNucleus = ProjectileResidual4Momentum / G4double( NumberOfProjectileParticipant );
    }
    //G4cout << "DeltaExcitationE DeltaPResidualNucleus " << DeltaExcitationE
    //       << " " << DeltaPResidualNucleus << G4endl;
    for ( G4int i = 0; i < NumberOfInvolvedNucleonsOfProjectile; ++i ) {
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

  }  // End of the condition: if ( HighEnergyInter )

  #ifdef debugFTFmodel
  G4cout << "End GetResiduals -----------------" << G4endl;
  #endif

}


//============================================================================

G4ThreeVector G4FTFModel::GaussianPt( G4double AveragePt2, G4double maxPtSquare ) const {

  G4double Pt2( 0.0 ), Pt( 0.0 );

  if (AveragePt2 > 0.0) {
    const G4double ymax = maxPtSquare/AveragePt2;
    if ( ymax < 200. ) {
      Pt2 = -AveragePt2 * G4Log( 1.0 + G4UniformRand() * ( G4Exp( -ymax ) -1.0 ) );
    } else {
      Pt2 = -AveragePt2 * G4Log( 1.0 - G4UniformRand() ); 
    }
    Pt = std::sqrt( Pt2 );
  }

  G4double phi = G4UniformRand() * twopi;
 
  return G4ThreeVector( Pt*std::cos(phi), Pt*std::sin(phi), 0.0 );    
}

//============================================================================

G4bool G4FTFModel::
ComputeNucleusProperties( G4V3DNucleus* nucleus,               // input parameter 
                          G4LorentzVector& nucleusMomentum,    // input & output parameter
                          G4LorentzVector& residualMomentum,   // input & output parameter
                          G4double& sumMasses,                 // input & output parameter
                          G4double& residualExcitationEnergy,  // input & output parameter
                          G4double& residualMass,              // input & output parameter
                          G4int& residualMassNumber,           // input & output parameter
                          G4int& residualCharge ) {            // input & output parameter

  // This method, which is called only by PutOnMassShell, computes some nucleus properties for:
  // -  either the target nucleus (which is never an antinucleus): this for any kind
  //    of hadronic interaction (hadron-nucleus, nucleus-nucleus, antinucleus-nucleus);
  // -  or the projectile nucleus or antinucleus: this only in the case of nucleus-nucleus
  //    or antinucleus-nucleus interaction.
  // This method assumes that the all the parameters have been initialized by the caller;
  // the action of this method consists in modifying all these parameters, except the
  // first one. The return value is "false" only in the case the pointer to the nucleus
  // is null.

  if ( ! nucleus ) return false;

  G4double ExcitationEnergyPerWoundedNucleon = 
    theParameters->GetExcitationEnergyPerWoundedNucleon();

  // Loop over the nucleons of the nucleus. 
  // The nucleons that have been involved in the interaction (either from Glauber or
  // Reggeon Cascading) will be candidate to be emitted.
  // All the remaining nucleons will be the nucleons of the candidate residual nucleus.
  // The variable sumMasses is the amount of energy corresponding to:
  //     1. transverse mass of each involved nucleon
  //     2. 20.0*MeV separation energy for each involved nucleon
  //     3. transverse mass of the residual nucleus
  // In this first evaluation of sumMasses, the excitation energy of the residual nucleus
  // (residualExcitationEnergy, estimated by adding a constant value to each involved
  // nucleon) is not taken into account.
  G4int residualNumberOfLambdas = 0;  // Projectile nucleus and its residual can be a hypernucleus
  G4Nucleon* aNucleon = 0;
  nucleus->StartLoop();
  while ( ( aNucleon = nucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
    nucleusMomentum += aNucleon->Get4Momentum();
    if ( aNucleon->AreYouHit() ) {  // Involved nucleons
      // Consider in sumMasses the nominal, i.e. on-shell, masses of the nucleons
      // (not the current masses, which could be different because the nucleons are off-shell).
      sumMasses += std::sqrt( sqr( aNucleon->GetDefinition()->GetPDGMass() ) 
                              +  aNucleon->Get4Momentum().perp2() );                     
      sumMasses += 20.0*MeV;  // Separation energy for a nucleon

      //residualExcitationEnergy += ExcitationEnergyPerWoundedNucleon;  // In G4 10.1
      residualExcitationEnergy += -ExcitationEnergyPerWoundedNucleon*G4Log( G4UniformRand() );

      residualMassNumber--;
      // The absolute value below is needed only in the case of anti-nucleus.
      residualCharge -= std::abs( G4int( aNucleon->GetDefinition()->GetPDGCharge() ) );
    } else {   // Spectator nucleons
      residualMomentum += aNucleon->Get4Momentum();
      if ( aNucleon->GetDefinition() == G4Lambda::Definition()  ||
	   aNucleon->GetDefinition() == G4AntiLambda::Definition() ) {
	++residualNumberOfLambdas;
      }
    }
  }
  #ifdef debugPutOnMassShell
  G4cout << "ExcitationEnergyPerWoundedNucleon " << ExcitationEnergyPerWoundedNucleon << G4endl
         << "\t Residual Charge, MassNumber (Number of Lambdas)" << residualCharge << " "
	 << residualMassNumber << " (" << residualNumberOfLambdas << ") "
         << G4endl << "\t Initial Momentum " << nucleusMomentum
         << G4endl << "\t Residual Momentum   " << residualMomentum << G4endl;
  #endif
  residualMomentum.setPz( 0.0 ); 
  residualMomentum.setE( 0.0 );
  if ( residualMassNumber == 0 ) {
    residualMass = 0.0;
    residualExcitationEnergy = 0.0;
  } else {
    if ( residualMassNumber == 1 ) {
      if ( std::abs( residualCharge ) == 1 ) {
        residualMass = G4Proton::Definition()->GetPDGMass();
      } else if ( residualNumberOfLambdas == 1 ) {
        residualMass = G4Lambda::Definition()->GetPDGMass();
      } else {
        residualMass = G4Neutron::Definition()->GetPDGMass();
      } 
      residualExcitationEnergy = 0.0;
    } else {
      if ( residualNumberOfLambdas > 0 ) {
        if ( residualMassNumber == 2 ) {
	  residualMass = G4Lambda::Definition()->GetPDGMass();
          if ( std::abs( residualCharge ) == 1 ) {      // lambda + proton
            residualMass += G4Proton::Definition()->GetPDGMass();
	  } else if ( residualNumberOfLambdas == 1 ) {  // lambda + neutron
	    residualMass += G4Neutron::Definition()->GetPDGMass();
          } else {                                      // lambda + lambda
	    residualMass += G4Lambda::Definition()->GetPDGMass();
          }
        } else {
	  residualMass = G4HyperNucleiProperties::GetNuclearMass( residualMassNumber, std::abs( residualCharge ),
	  							  residualNumberOfLambdas );
        }
      } else {
        residualMass = G4ParticleTable::GetParticleTable()->GetIonTable()->
	               GetIonMass( std::abs( residualCharge ), residualMassNumber );
      }
    }
    residualMass += residualExcitationEnergy;
  }
  sumMasses += std::sqrt( sqr( residualMass ) + residualMomentum.perp2() );
  return true;
}


//============================================================================

G4bool G4FTFModel::
GenerateDeltaIsobar( const G4double sqrtS,                  // input parameter
                     const G4int numberOfInvolvedNucleons,  // input parameter
                     G4Nucleon* involvedNucleons[],         // input & output parameter
                     G4double& sumMasses ) {                // input & output parameter

  // This method, which is called only by PutOnMassShell, check whether is possible to
  // re-interpret some of the involved nucleons as delta-isobars:
  // - either by replacing a proton (2212) with a Delta+ (2214),
  // - or by replacing a neutron (2112) with a Delta0 (2114).
  // The on-shell mass of these delta-isobars is ~1232 MeV, so  ~292-294 MeV  heavier than
  // the corresponding nucleon on-shell mass. However  400.0*MeV  is considered to estimate
  // the max number of deltas compatible with the available energy.
  // The delta-isobars are considered with the same transverse momentum as their
  // corresponding nucleons.
  // This method assumes that all the parameters have been initialized by the caller;
  // the action of this method consists in modifying (eventually) involveNucleons and
  // sumMasses. The return value is "false" only in the case that the input parameters
  // have unphysical values.

  if ( sqrtS < 0.0  ||  numberOfInvolvedNucleons <= 0  ||  sumMasses < 0.0 ) return false;

  const G4double probDeltaIsobar = 0.05;

  G4int maxNumberOfDeltas = G4int( (sqrtS - sumMasses)/(400.0*MeV) );
  G4int numberOfDeltas = 0;

  for ( G4int i = 0; i < numberOfInvolvedNucleons; ++i ) {

    if ( G4UniformRand() < probDeltaIsobar  &&  numberOfDeltas < maxNumberOfDeltas ) {
      numberOfDeltas++;
      if ( ! involvedNucleons[i] ) continue;
      // Skip any eventual lambda (that can be present in a projectile hypernucleus)
      if ( involvedNucleons[i]->GetDefinition() == G4Lambda::Definition()  ||
	   involvedNucleons[i]->GetDefinition() == G4AntiLambda::Definition() ) continue;
      G4VSplitableHadron* splitableHadron = involvedNucleons[i]->GetSplitableHadron();
      G4double massNuc = std::sqrt( sqr( splitableHadron->GetDefinition()->GetPDGMass() )
                                    + splitableHadron->Get4Momentum().perp2() );
      // The absolute value below is needed in the case of an antinucleus. 
      G4int pdgCode = std::abs( splitableHadron->GetDefinition()->GetPDGEncoding() );
      const G4ParticleDefinition* old_def = splitableHadron->GetDefinition();
      G4int newPdgCode = pdgCode/10; newPdgCode = newPdgCode*10 + 4; // Delta
      if ( splitableHadron->GetDefinition()->GetPDGEncoding() < 0 ) newPdgCode *= -1;
      const G4ParticleDefinition* ptr = 
        G4ParticleTable::GetParticleTable()->FindParticle( newPdgCode );
      splitableHadron->SetDefinition( ptr );
      G4double massDelta = std::sqrt( sqr( splitableHadron->GetDefinition()->GetPDGMass() )
                                      + splitableHadron->Get4Momentum().perp2() );
      //G4cout << i << " " << sqrtS/GeV << " " << sumMasses/GeV << " " << massDelta/GeV
      //       << " " << massNuc << G4endl;
      if ( sqrtS < sumMasses + massDelta - massNuc ) {  // Change cannot be accepted!
        splitableHadron->SetDefinition( old_def );
        break;
      } else {  // Change is accepted
        sumMasses += ( massDelta - massNuc );
      }
    } 
  }

  return true;
}


//============================================================================

G4bool G4FTFModel::
SamplingNucleonKinematics( G4double averagePt2,                   // input parameter
                           const G4double maxPt2,                 // input parameter
                           G4double dCor,                         // input parameter
                           G4V3DNucleus* nucleus,                 // input parameter
                           const G4LorentzVector& pResidual,      // input parameter
                           const G4double residualMass,           // input parameter
                           const G4int residualMassNumber,        // input parameter
                           const G4int numberOfInvolvedNucleons,  // input parameter 
                           G4Nucleon* involvedNucleons[],         // input & output parameter
                           G4double& mass2 ) {                    // output parameter

  // This method, which is called only by PutOnMassShell, does the sampling of:
  // -  either the target nucleons: this for any kind of hadronic interactions
  //    (hadron-nucleus, nucleus-nucleus, antinucleus-nucleus);
  // -  or the projectile nucleons or antinucleons: this only in the case of
  //    nucleus-nucleus or antinucleus-nucleus interactions, respectively.
  // This method assumes that all the parameters have been initialized by the caller;
  // the action of this method consists in changing the properties of the nucleons
  // whose pointers are in the vector involvedNucleons, as well as changing the
  // variable mass2.
#ifdef debugPutOnMassShell
  G4cout << "G4FTFModel::SamplingNucleonKinematics:" << G4endl;
  G4cout << " averagePt2= " << averagePt2 << " maxPt2= " << maxPt2
	 << " dCor= " << dCor << " resMass(GeV)= " << residualMass/GeV
	 << " resMassN= " << residualMassNumber 
         << " nNuc= " << numberOfInvolvedNucleons
	 << " lv= " << pResidual << G4endl;   
#endif

  if ( ! nucleus  ||  numberOfInvolvedNucleons < 1 ) return false;

  if ( residualMassNumber == 0  &&  numberOfInvolvedNucleons == 1 ) {
    dCor = 0.0; 
    averagePt2 = 0.0;
  } 

  G4bool success = true;                            

  G4double SumMasses = residualMass; 
  G4double invN = 1.0 / (G4double)numberOfInvolvedNucleons;

  // to avoid problems due to precision lost a tolerance is added
  const G4double eps = 1.e-10;  
  const G4int maxNumberOfLoops = 1000;
  G4int loopCounter = 0;
  do {

    success = true;

    // Sampling of nucleon Pt 
    G4ThreeVector ptSum( 0.0, 0.0, 0.0 );
    if( averagePt2 > 0.0 ) {
      for ( G4int i = 0; i < numberOfInvolvedNucleons; ++i ) {
	G4Nucleon* aNucleon = involvedNucleons[i];
	if ( ! aNucleon ) continue;
	G4ThreeVector tmpPt = GaussianPt( averagePt2, maxPt2 );
	ptSum += tmpPt;
	G4LorentzVector tmp( tmpPt.x(), tmpPt.y(), 0.0, 0.0 );
	aNucleon->SetMomentum( tmp );
      }
    }

    G4double deltaPx = ( ptSum.x() - pResidual.x() )*invN;
    G4double deltaPy = ( ptSum.y() - pResidual.y() )*invN;

    SumMasses = residualMass;
    for ( G4int i = 0; i < numberOfInvolvedNucleons; ++i ) {
      G4Nucleon* aNucleon = involvedNucleons[i];
      if ( ! aNucleon ) continue;
      G4double px = aNucleon->Get4Momentum().px() - deltaPx;
      G4double py = aNucleon->Get4Momentum().py() - deltaPy;
      G4double MtN = std::sqrt( sqr( aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass() )
                                + sqr( px ) + sqr( py ) );
      SumMasses += MtN;
      G4LorentzVector tmp( px, py, 0.0, MtN);
      aNucleon->SetMomentum( tmp );
    }

    // Sampling X of nucleon
    G4double xSum = 0.0;

    for ( G4int i = 0; i < numberOfInvolvedNucleons; ++i ) {
      G4Nucleon* aNucleon = involvedNucleons[i];
      if ( ! aNucleon ) continue;
     
      G4double x = 0.0;
      if( 0.0 != dCor ) {
	G4ThreeVector tmpX = GaussianPt( dCor*dCor, 1.0 );
        x = tmpX.x();
      } 
      x += aNucleon->Get4Momentum().e()/SumMasses;
      if ( x < -eps  ||  x > 1.0 + eps ) { 
        success = false; 
        break;
      }
      x = std::min(1.0, std::max(x, 0.0));
      xSum += x;
      // The energy is in the lab (instead of cms) frame but it will not be used

      G4LorentzVector tmp( aNucleon->Get4Momentum().x(), 
                           aNucleon->Get4Momentum().y(), 
                           x, aNucleon->Get4Momentum().e() );
      aNucleon->SetMomentum( tmp );
    }

    if ( xSum < -eps || xSum > 1.0 + eps ) success = false;
    if ( ! success ) continue;

    G4double delta = ( residualMassNumber == 0 ) ? std::min( xSum - 1.0, 0.0 )*invN : 0.0;

    xSum = 1.0;
    mass2 = 0.0;
    for ( G4int i = 0; i < numberOfInvolvedNucleons; ++i ) {
      G4Nucleon* aNucleon = involvedNucleons[i];
      if ( ! aNucleon ) continue;
      G4double x = aNucleon->Get4Momentum().pz() - delta;
      xSum -= x;

      if ( residualMassNumber == 0 ) {
        if ( x <= -eps || x > 1.0 + eps ) {
          success = false; 
          break;
        }
      } else {
        if ( x <= -eps || x > 1.0 + eps || xSum <= -eps || xSum > 1.0 + eps ) {
          success = false; 
          break;
        }
      } 
      x = std::min( 1.0, std::max(x, eps) );

      mass2 += sqr( aNucleon->Get4Momentum().e() ) / x;

      G4LorentzVector tmp( aNucleon->Get4Momentum().px(), aNucleon->Get4Momentum().py(), 
                           x, aNucleon->Get4Momentum().e() );
      aNucleon->SetMomentum( tmp );
    }
    if ( ! success ) continue;
    xSum = std::min( 1.0, std::max(xSum, eps) );

    if ( residualMassNumber > 0 ) mass2 += ( sqr( residualMass ) + pResidual.perp2() ) / xSum;
    
    #ifdef debugPutOnMassShell
    G4cout << "success: " << success << " Mt(GeV)= " 
	   << std::sqrt( mass2 )/GeV << G4endl;
    #endif

  } while ( ( ! success ) &&
            ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */
  return ( loopCounter < maxNumberOfLoops );
}


//============================================================================

G4bool G4FTFModel::
CheckKinematics( const G4double sValue,                 // input parameter
                 const G4double sqrtS,                  // input parameter
                 const G4double projectileMass2,        // input parameter
                 const G4double targetMass2,            // input parameter
                 const G4double nucleusY,               // input parameter
                 const G4bool isProjectileNucleus,      // input parameter
                 const G4int numberOfInvolvedNucleons,  // input parameter 
                 G4Nucleon* involvedNucleons[],         // input parameter
                 G4double& targetWminus,                // output parameter
                 G4double& projectileWplus,             // output parameter
                 G4bool& success ) {                    // input & output parameter

  // This method, which is called only by PutOnMassShell, checks whether the
  // kinematics is acceptable or not.
  // This method assumes that all the parameters have been initialized by the caller;
  // notice that the input boolean parameter isProjectileNucleus is meant to be true
  // only in the case of nucleus or antinucleus projectile.
  // The action of this method consists in computing targetWminus and projectileWplus
  // and setting the parameter success to false in the case that the kinematics should
  // be rejeted.

  G4double decayMomentum2 = sqr( sValue ) + sqr( projectileMass2 ) + sqr( targetMass2 )
                            - 2.0*( sValue*( projectileMass2 + targetMass2 ) 
                                    + projectileMass2*targetMass2 );
  targetWminus = ( sValue - projectileMass2 + targetMass2 + std::sqrt( decayMomentum2 ) )
                 / 2.0 / sqrtS;
  projectileWplus = sqrtS - targetMass2/targetWminus;
  G4double projectilePz = projectileWplus/2.0 - projectileMass2/2.0/projectileWplus;
  G4double projectileE  = projectileWplus/2.0 + projectileMass2/2.0/projectileWplus;
  G4double projectileY  = 0.5 * G4Log( (projectileE + projectilePz)/
                                       (projectileE - projectilePz) );
  G4double targetPz = -targetWminus/2.0 + targetMass2/2.0/targetWminus;
  G4double targetE  =  targetWminus/2.0 + targetMass2/2.0/targetWminus;
  G4double targetY  = 0.5 * G4Log( (targetE + targetPz)/(targetE - targetPz) );

  #ifdef debugPutOnMassShell
  G4cout << "decayMomentum2 " << decayMomentum2 << G4endl 
         << "\t targetWminus projectileWplus " << targetWminus << " " << projectileWplus << G4endl
         << "\t projectileY targetY " << projectileY << " " << targetY << G4endl;
  if ( isProjectileNucleus ) {
    G4cout << "Order# of Wounded nucleon i, nucleon Y proj Y nuclY - proj Y " << G4endl;
  } else {
    G4cout << "Order# of Wounded nucleon i, nucleon Y targ Y nuclY - targ Y " << G4endl;
  }
  G4cout << G4endl;
  #endif
  
  for ( G4int i = 0; i < numberOfInvolvedNucleons; ++i ) {
    G4Nucleon* aNucleon = involvedNucleons[i];
    if ( ! aNucleon ) continue;
    G4LorentzVector tmp = aNucleon->Get4Momentum();
    G4double mt2 = sqr( tmp.x() ) + sqr( tmp.y() ) +
                   sqr( aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass() );
    G4double x = tmp.z();
    G4double pz = -targetWminus*x/2.0 + mt2/(2.0*targetWminus*x);
    G4double e =   targetWminus*x/2.0 + mt2/(2.0*targetWminus*x);
    if ( isProjectileNucleus ) {
      pz = projectileWplus*x/2.0 - mt2/(2.0*projectileWplus*x);
      e =  projectileWplus*x/2.0 + mt2/(2.0*projectileWplus*x);
    }
    G4double nucleonY = 0.5 * G4Log( (e + pz)/(e - pz) ); 

    #ifdef debugPutOnMassShell
    if( isProjectileNucleus ) {
      G4cout << " " << i << " " << nucleonY << " " << projectileY << " " <<nucleonY - projectileY << G4endl;
    } else {
      G4cout << " " << i << " " << nucleonY << " " << targetY     << " " <<nucleonY - targetY << G4endl;
    }
    G4cout << G4endl;
    #endif

    if ( std::abs( nucleonY - nucleusY ) > 2  ||  
         ( isProjectileNucleus  &&  targetY > nucleonY )  ||
         ( ! isProjectileNucleus  &&  projectileY < nucleonY ) ) {
      success = false; 
      break;
    } 
  }
  return true;
}  

  
//============================================================================

G4bool G4FTFModel::
FinalizeKinematics( const G4double w,                            // input parameter
                    const G4bool isProjectileNucleus,            // input parameter
                    const G4LorentzRotation& boostFromCmsToLab,  // input parameter
                    const G4double residualMass,                 // input parameter
                    const G4int residualMassNumber,              // input parameter
                    const G4int numberOfInvolvedNucleons,        // input parameter 
                    G4Nucleon* involvedNucleons[],               // input & output parameter
	            G4LorentzVector& residual4Momentum ) {       // output parameter

  // This method, which is called only by PutOnMassShell, finalizes the kinematics:
  // this method is called when we are sure that the sampling of the kinematics is
  // acceptable.
  // This method assumes that all the parameters have been initialized by the caller;
  // notice that the input boolean parameter isProjectileNucleus is meant to be true
  // only in the case of nucleus or antinucleus projectile: this information is needed
  // because the sign of pz (in the center-of-mass frame) in this case is opposite
  // with respect to the case of a normal hadron projectile.
  // The action of this method consists in modifying the momenta of the nucleons
  // (in the lab frame) and computing the residual 4-momentum (in the center-of-mass
  // frame).

  G4ThreeVector residual3Momentum( 0.0, 0.0, 1.0 );

  for ( G4int i = 0; i < numberOfInvolvedNucleons; ++i ) {
    G4Nucleon* aNucleon = involvedNucleons[i];
    if ( ! aNucleon ) continue;
    G4LorentzVector tmp = aNucleon->Get4Momentum();
    residual3Momentum -= tmp.vect();
    G4double mt2 = sqr( tmp.x() ) + sqr( tmp.y() ) +
                   sqr( aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass() );
    G4double x = tmp.z();
    G4double pz = -w * x / 2.0  +  mt2 / ( 2.0 * w * x );
    G4double e  =  w * x / 2.0  +  mt2 / ( 2.0 * w * x );
    // Reverse the sign of pz in the case of nucleus or antinucleus projectile
    if ( isProjectileNucleus ) pz *= -1.0;
    tmp.setPz( pz ); 
    tmp.setE( e );
    tmp.transform( boostFromCmsToLab );
    aNucleon->SetMomentum( tmp );
    G4VSplitableHadron* splitableHadron = aNucleon->GetSplitableHadron();
    splitableHadron->Set4Momentum( tmp );
  }

  G4double residualMt2 = sqr( residualMass ) + sqr( residual3Momentum.x() )
                       + sqr( residual3Momentum.y() );

  #ifdef debugPutOnMassShell
  if ( isProjectileNucleus ) {
    G4cout << "Wminus Proj and residual3Momentum.z() " << w << " " << residual3Momentum.z() << G4endl;
  } else {
    G4cout << "Wplus  Targ and residual3Momentum.z() " << w << " " << residual3Momentum.z() << G4endl;
  }
  #endif

  G4double residualPz = 0.0;
  G4double residualE  = 0.0;
  if ( residualMassNumber != 0 ) {
    residualPz = -w * residual3Momentum.z() / 2.0 + 
                  residualMt2 / ( 2.0 * w * residual3Momentum.z() );
    residualE  =  w * residual3Momentum.z() / 2.0 + 
                  residualMt2 / ( 2.0 * w * residual3Momentum.z() );
    // Reverse the sign of residualPz in the case of nucleus or antinucleus projectile
    if ( isProjectileNucleus ) residualPz *= -1.0;
  }

  residual4Momentum.setPx( residual3Momentum.x() );
  residual4Momentum.setPy( residual3Momentum.y() );
  residual4Momentum.setPz( residualPz ); 
  residual4Momentum.setE( residualE );

  return true;
}

  
//============================================================================

void G4FTFModel::ModelDescription( std::ostream& desc ) const {
  desc << "                 FTF (Fritiof) Model               \n" 
       << "The FTF model is based on the well-known FRITIOF   \n"
       << "model (B. Andersson et al., Nucl. Phys. B281, 289  \n"
       << "(1987)). Its first program implementation was given\n"
       << "by B. Nilsson-Almquist and E. Stenlund (Comp. Phys.\n"
       << "Comm. 43, 387 (1987)). The Fritiof model assumes   \n"
       << "that all hadron-hadron interactions are binary     \n"
       << "reactions, h_1+h_2->h_1'+h_2' where h_1' and h_2'  \n"
       << "are excited states of the hadrons with continuous  \n"
       << "mass spectra. The excited hadrons are considered as\n"
       << "QCD-strings, and the corresponding LUND-string     \n"
       << "fragmentation model is applied for a simulation of \n"
       << "their decays.                                      \n"
       << "   The Fritiof model assumes that in the course of \n"
       << "a hadron-nucleus interaction a string originated   \n"
       << "from the projectile can interact with various intra\n"
       << "nuclear nucleons and becomes into highly excited   \n"
       << "states. The probability of multiple interactions is\n"
       << "calculated in the Glauber approximation. A cascading\n"
       << "of secondary particles was neglected as a rule. Due\n"
       << "to these, the original Fritiof model fails to des- \n"
       << "cribe a nuclear destruction and slow particle spectra.\n"
       << "   In order to overcome the difficulties we enlarge\n"
       << "the model by the reggeon theory inspired model of  \n"
       << "nuclear desctruction (Kh. Abdel-Waged and V.V. Uzhi-\n"
       << "nsky, Phys. Atom. Nucl. 60, 828 (1997); Yad. Fiz. 60, 925\n"
       << "(1997)). Momenta of the nucleons ejected from a nuc-\n"
       << "leus in the reggeon cascading are sampled according\n"
       << "to a Fermi motion algorithm presented in (EMU-01   \n"
       << "Collaboration (M.I. Adamovich et al.) Zeit. fur Phys.\n"
       << "A358, 337 (1997)).                                 \n"
       << "   New features were also added to the Fritiof model\n"
       << "implemented in Geant4: a simulation of elastic had-\n"
       << "ron-nucleon scatterings, a simulation of binary \n"
       << "reactions like NN>NN* in hadron-nucleon interactions,\n"
       << "a separate simulation of single diffractive and non-\n"
       << " diffractive events. These allowed to describe after\n"
       << "model parameter tuning a wide set of experimental  \n"
       << "data.                                              \n";
}

