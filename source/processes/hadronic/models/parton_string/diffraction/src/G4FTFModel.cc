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
// $Id: G4FTFModel.cc 102029 2016-12-16 14:53:08Z gcosmo $
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
#include "G4KineticTrack.hh"

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
  G4VPartonStringModel::SetThisPointer( this );
  theParameters = 0;
  NumberOfInvolvedNucleonsOfTarget = 0;
  NumberOfInvolvedNucleonsOfProjectile= 0;
  for ( G4int i = 0; i < 250; i++ ) {
    TheInvolvedNucleonsOfTarget[i] = 0;
    TheInvolvedNucleonsOfProjectile[i] = 0;
  }

//  LowEnergyLimit = 2000.0*MeV;   // Uzhi March 2015
  LowEnergyLimit = 1000.0*MeV;     // Uzhi May 2015

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

struct DeleteVSplitableHadron { void operator()( G4VSplitableHadron* aH ) { delete aH; } };


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

  theParticipants.Clean();

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
      while ( ( aNucleon = theParticipants.theProjectileNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
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
//theParticipants.Init( aNucleus.GetA_asInt(), 0 ); //For h+neutron // Uzhi March 2016

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

//  if ( (std::abs( theProjectile.GetDefinition()->GetBaryonNumber() ) <= 1 ) &&        // Uzhi 29.05.2015
//       (aNucleus.GetA_asInt() < 2) ) theParameters->SetProbabilityOfElasticScatt(0.);

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
  //G4int Uzhi; G4cin >> Uzhi;
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
          targetSplitable->SetStatus( 3 );                       // 2->3  Uzhi Oct 2014
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

//  for ( G4int InvPN = 0; InvPN < NumberOfInvolvedNucleonsOfProjectile; InvPN++ ) { 
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
  if ( GetProjectileNucleus() ) {
    isProjectileNucleus = true;
  }

  #ifdef debugPutOnMassShell
  G4cout << "PutOnMassShell start " << G4endl;
  if ( isProjectileNucleus ) {
    G4cout << "PutOnMassShell for Nucleus_Nucleus " << G4endl;
  }
  #endif

  G4LorentzVector Pprojectile( theProjectile.GetMomentum(), theProjectile.GetTotalEnergy() );
  if ( Pprojectile.z() < 0.0 ) {
    return false;
  }

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

  if ( SqrtS < SumMasses ) {
    return false;  // It is impossible to simulate after putting nuclear nucleons on mass-shell.
  }

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
    if ( isProjectileNucleus ) {
      ProjectileResidualExcitationEnergy = 0.0;
    }
    TargetResidualExcitationEnergy = 0.0;
  }

  TargetResidualMass += TargetResidualExcitationEnergy;
  if ( isProjectileNucleus ) {
    PrResidualMass += ProjectileResidualExcitationEnergy;
  }

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
    isOk = isOk  &&
           GenerateDeltaIsobar( SqrtS, NumberOfInvolvedNucleonsOfTarget,
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
  if ( Ptmp.pz() <= 0.0 ) {  // "String" moving backwards in c.m.s., abort collision!
    return false; 
  }

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
  if ( isProjectileNucleus ) {
    DcorP = theParameters->GetDofNuclearDestruction() / thePrNucleus->GetMassNumber();
  }
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
  G4double ScaleFactor = 1.0;
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
      isOk = isOk  &&
             SamplingNucleonKinematics( AveragePt2, maxPtSquare, DcorT, 
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
    isOk = isOk  &&
           CheckKinematics( S, SqrtS, M2proj, M2target, YtargetNucleus, false, 
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
             ( 1.0 - target->GetSoftCollisionCount()     / MaxNumOfInelCollisions )  *  // Uzhi March 2015
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

            Successfull = theElastic->ElasticScattering( projectile, target, theParameters )
                          &&  Successfull; 
//                          ||  Successfull; 
            #ifdef debugBuildString
            G4cout << "FTF excitation Non Successfull -> Elastic scattering " 
                   << Successfull << G4endl;
            #endif
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
// //  Uzhi March 2016       
        // Skipping possible interactions of the annihilated nucleons 
        while ( theParticipants.Next() ) {   /* Loop checking, 10.08.2015, A.Ribon */  
          G4InteractionContent& acollision = theParticipants.GetInteraction();
          G4VSplitableHadron* NextProjectileNucleon = acollision.GetProjectile();
          G4VSplitableHadron* NextTargetNucleon = acollision.GetTarget();
          if ( projectile == NextProjectileNucleon  ||  target == NextTargetNucleon ) {
            acollision.SetStatus( 0 );
          }
        }
// // Uzhi March 2016
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
/* Uzhi March 2016
if(target->GetStatus() == 4){
        // Skipping possible interactions of the annihilated nucleons 
        while ( theParticipants.Next() ) {    
          G4InteractionContent& acollision = theParticipants.GetInteraction();
          G4VSplitableHadron* NextProjectileNucleon = acollision.GetProjectile();
          G4VSplitableHadron* NextTargetNucleon = acollision.GetTarget();
          if ( target == NextTargetNucleon ) {acollision.SetStatus( 0 );}
        }
}
        theParticipants.StartLoop(); 
        for ( G4int I = 0; I < CurrentInteraction; I++ ) theParticipants.Next();
*/ //Uzhi March 2016
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
//Uzhi    G4double TResidualExcitationEnergy = TargetResidualExcitationEnergy + 
//                                         ExcitationEnergyPerWoundedNucleon;
    G4double TResidualExcitationEnergy = TargetResidualExcitationEnergy -        // Uzhi April 2015
                                         ExcitationEnergyPerWoundedNucleon*G4Log( G4UniformRand());
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

        SumMasses = SqrtS - TResidualExcitationEnergy;
        //TResidualExcitationEnergy =0.0;             
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

    const G4int maxNumberOfLoops = 1000;
    G4int loopCounter = 0;
    do { // while ( ! OuterSuccess )
      OuterSuccess = true;

      const G4int maxNumberOfTries = 10000;
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
          const G4int maxNumberOfInnerLoops = 1000;
          G4int innerLoopCounter = 0;
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
          } while ( ( ! InerSuccess ) &&
                    ++innerLoopCounter < maxNumberOfInnerLoops );  /* Loop checking, 10.08.2015, A.Ribon */
          if ( innerLoopCounter >= maxNumberOfInnerLoops ) {
            #ifdef debugAdjust
            G4cout << "BAD situation: forced exit of the inner while loop!" << G4endl;
            #endif
            return false;
          }

        } else {
          XminusNucleon  = 1.0;
          XminusResidual = 1.0;  // It must be 0, but in the case calculation of Pz,
                                 // E is problematic.
        }

        M2target = ( sqr( TNucleonMass ) + PtNucleon.mag2() ) / XminusNucleon + 
                   ( sqr( TResidualMass ) + PtResidual.mag2() ) / XminusResidual;

      } while ( ( SqrtS < Mprojectile + std::sqrt( M2target) ) &&
                ++NumberOfTries < maxNumberOfTries );  /* Loop checking, 10.08.2015, A.Ribon */
      if ( NumberOfTries >= maxNumberOfTries ) {
        #ifdef debugAdjust
        G4cout << "BAD situation: forced exit of the intermediate while loop!" << G4endl;
        #endif
        return false;
      }

      G4double DecayMomentum2 = sqr( S ) + sqr( M2projectile ) + sqr( M2target )
                                - 2.0*S*M2projectile - 2.0*S*M2target - 2.0*M2projectile*M2target;

      WminusTarget = ( S - M2projectile + M2target + std::sqrt( DecayMomentum2 ) ) / 2.0 / SqrtS;
      WplusProjectile = SqrtS - M2target / WminusTarget;

      G4double Pzprojectile = WplusProjectile/2.0 - M2projectile/2.0/WplusProjectile;
      G4double Eprojectile  = WplusProjectile/2.0 + M2projectile/2.0/WplusProjectile;
      G4double Yprojectile  = 0.5 * G4Log( (Eprojectile + Pzprojectile) /
                                           (Eprojectile - Pzprojectile) );

      #ifdef debugAdjust
      G4cout << "DecayMomentum2 " << DecayMomentum2 << G4endl
             << "WminusTarget WplusProjectile " << WminusTarget << " " << WplusProjectile
             << G4endl << "Yprojectile " << Yprojectile << G4endl;
      #endif

      G4double Mt2 = sqr( TNucleonMass ) + PtNucleon.mag2();
      G4double Pz = -WminusTarget*XminusNucleon/2.0 + Mt2/(2.0*WminusTarget*XminusNucleon);
      G4double E  =  WminusTarget*XminusNucleon/2.0 + Mt2/(2.0*WminusTarget*XminusNucleon);
      G4double YtargetNucleon = 0.5 * G4Log( (E + Pz)/(E - Pz) ); 

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

    } while ( ( ! OuterSuccess ) &&
              ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */
    if ( loopCounter >= maxNumberOfLoops ) {
      #ifdef debugAdjust
      G4cout << "BAD situation: forced exit of the while loop!" << G4endl;
      #endif
      return false;
    }

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
//Uzhi    G4double TResidualExcitationEnergy = ProjectileResidualExcitationEnergy + 
//                                         ExcitationEnergyPerWoundedNucleon;
    G4double TResidualExcitationEnergy = ProjectileResidualExcitationEnergy -   // Uzhi April 2015 
                                         ExcitationEnergyPerWoundedNucleon*G4Log( G4UniformRand());
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

    const G4int maxNumberOfLoops = 1000;
    G4int loopCounter = 0;
    do { // while ( ! OuterSuccess )
  
      OuterSuccess = true;
      const G4int maxNumberOfTries = 10000;
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

        M2projectile = sqr( Mprojectile );
        if ( SqrtS < Mtarget + Mprojectile ) {
          OuterSuccess = false; 
          continue;
        }

        G4double Xcenter = std::sqrt( sqr( TNucleonMass ) + PtNucleon.mag2() ) / Mprojectile;

        G4bool InerSuccess = true;
        if ( ProjectileResidualMassNumber > 1 ) {
          const G4int maxNumberOfInnerLoops = 1000;
          G4int innerLoopCounter = 0;
          do {
            InerSuccess = true;
            G4ThreeVector tmpX = GaussianPt( DcorP*DcorP, 1.0 );
            XplusNucleon = Xcenter + tmpX.x();
            if ( XplusNucleon <= 0.0  ||  XplusNucleon >= 1.0 ) { 
              InerSuccess = false; 
              continue;
            }
            XplusResidual = 1.0 - XplusNucleon;
          } while ( ( ! InerSuccess ) &&
                    ++innerLoopCounter < maxNumberOfInnerLoops );  /* Loop checking, 10.08.2015, A.Ribon */
          if ( innerLoopCounter >= maxNumberOfInnerLoops ) {
            #ifdef debugAdjust
            G4cout << "BAD situation: forced exit of the inner while loop!" << G4endl;
            #endif
            return false;
          }

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

      } while ( ( SqrtS < Mtarget + std::sqrt( M2projectile ) ) &&
                ++NumberOfTries < maxNumberOfTries );  /* Loop checking, 10.08.2015, A.Ribon */
      if ( NumberOfTries >= maxNumberOfTries ) {
        #ifdef debugAdjust
        G4cout << "BAD situation: forced exit of the intermediate while loop!" << G4endl;
        #endif
        return false;
      }

      G4double DecayMomentum2 = sqr( S ) + sqr( M2projectile ) + sqr( M2target )
                                - 2.0*S*M2projectile - 2.0*S*M2target - 2.0*M2projectile*M2target;

      WplusProjectile = ( S + M2projectile - M2target + std::sqrt( DecayMomentum2 ) )/2.0/SqrtS;
      WminusTarget = SqrtS - M2projectile/WplusProjectile;

      G4double Pztarget = -WminusTarget/2.0 + M2target/2.0/WminusTarget;
      G4double Etarget =   WminusTarget/2.0 + M2target/2.0/WminusTarget;
      G4double Ytarget = 0.5 * G4Log( (Etarget + Pztarget)/(Etarget - Pztarget) );

      #ifdef debugAdjust
      G4cout << "DecayMomentum2 " << DecayMomentum2 << G4endl
             << "WminusTarget WplusProjectile " << WminusTarget << " " << WplusProjectile 
             << G4endl << "YtargetNucleon " << Ytarget << G4endl;
      #endif

      G4double Mt2 = sqr( TNucleonMass ) + PtNucleon.mag2();
      G4double Pz = WplusProjectile*XplusNucleon/2.0 - Mt2/(2.0*WplusProjectile*XplusNucleon);
      G4double E =  WplusProjectile*XplusNucleon/2.0 + Mt2/(2.0*WplusProjectile*XplusNucleon);
      G4double YprojectileNucleon = 0.5 * G4Log( (E + Pz)/(E - Pz) ); 

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
 
    } while ( ( ! OuterSuccess ) &&
              ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */
    if ( loopCounter >= maxNumberOfLoops ) {
      #ifdef debugAdjust
      G4cout << "BAD situation: forced exit of the while loop!" << G4endl;
      #endif
      return false;
    }

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
//Uzhi    G4double PResidualExcitationEnergy = ProjectileResidualExcitationEnergy +
//                                         ExcitationEnergyPerWoundedNucleon;
    G4double PResidualExcitationEnergy = ProjectileResidualExcitationEnergy - 
                                         ExcitationEnergyPerWoundedNucleon*G4Log( G4UniformRand());
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
//Uzhi    G4double TResidualExcitationEnergy = TargetResidualExcitationEnergy +
//                                         ExcitationEnergyPerWoundedNucleon;
    G4double TResidualExcitationEnergy = TargetResidualExcitationEnergy -   
                                         ExcitationEnergyPerWoundedNucleon*G4Log( G4UniformRand());
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
      Ptmp.setPx( 0.0 ); Ptmp.setPy( 0.0 ); Ptmp.setPz( 0.0 );
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

    const G4int maxNumberOfLoops = 1000;
    G4int loopCounter = 0;
    do { // while ( ! OuterSuccess )

      OuterSuccess = true;
      const G4int maxNumberOfTries = 10000;
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
        M2projectile = sqr( Mprojectile );

        G4double Mtarget = std::sqrt( sqr( TNucleonMass )  + PtNucleonT.mag2() ) + 
                           std::sqrt( sqr( TResidualMass ) + PtResidualT.mag2() );
        M2target = sqr( Mtarget );        

        if ( SqrtS < Mprojectile + Mtarget ) {
          OuterSuccess = false; 
          continue;
        }

        G4bool InerSuccess = true;

        if ( ProjectileResidualMassNumber > 1 ) { 
          const G4int maxNumberOfInnerLoops = 1000;
          G4int innerLoopCounter = 0;
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
          } while ( ( ! InerSuccess ) &&
                    ++innerLoopCounter < maxNumberOfInnerLoops );  /* Loop checking, 10.08.2015, A.Ribon */
          if ( innerLoopCounter >= maxNumberOfInnerLoops ) {
            #ifdef debugAdjust
            G4cout << "BAD situation: forced exit of the first inner while loop!" << G4endl;
            #endif
            return false;
          }

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

          const G4int maxNumberOfInnerLoops = 1000;
          G4int innerLoopCounter = 0;
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
          } while ( ( ! InerSuccess ) &&
                    ++innerLoopCounter < maxNumberOfInnerLoops );  /* Loop checking, 10.08.2015, A.Ribon */
          if ( innerLoopCounter >= maxNumberOfInnerLoops ) {
            #ifdef debugAdjust
            G4cout << "BAD situation: forced exit of the second inner while loop!" << G4endl;
            #endif
            return false;
          }
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

      } while ( ( SqrtS < std::sqrt( M2projectile ) + std::sqrt( M2target ) ) &&
                ++NumberOfTries < maxNumberOfTries );  /* Loop checking, 10.08.2015, A.Ribon */
      if ( NumberOfTries >= maxNumberOfTries ) {
        #ifdef debugAdjust
        G4cout << "BAD situation: forced exit of the intermediate while loop!" << G4endl;
        #endif
        return false;
      }

      G4double DecayMomentum2 = sqr( S ) + sqr( M2projectile ) + sqr( M2target )
                                - 2.0*S*M2projectile - 2.0*S*M2target - 2.0*M2projectile*M2target;

      WplusProjectile = ( S + M2projectile - M2target + std::sqrt( DecayMomentum2 ) )/2.0/SqrtS;
      WminusTarget = SqrtS - M2projectile/WplusProjectile;

      G4double Mt2 = sqr( PNucleonMass ) + PtNucleonP.mag2();
      G4double Pz = WplusProjectile*XplusNucleon/2.0 - Mt2/(2.0*WplusProjectile*XplusNucleon);
      G4double E =  WplusProjectile*XplusNucleon/2.0 + Mt2/(2.0*WplusProjectile*XplusNucleon);
      G4double YprojectileNucleon = 0.5 * G4Log( (E + Pz)/(E - Pz) );

      Mt2 = sqr( TNucleonMass ) + PtNucleonT.mag2();
      Pz = -WminusTarget*XminusNucleon/2.0 + Mt2/(2.0*WminusTarget*XminusNucleon);
      E =   WminusTarget*XminusNucleon/2.0 + Mt2/(2.0*WminusTarget*XminusNucleon);
      G4double YtargetNucleon = 0.5 * G4Log( (E + Pz)/(E - Pz) ); 

      if ( std::abs( YtargetNucleon - YtargetNucleus ) > 2         || 
           std::abs( YprojectileNucleon - YprojectileNucleus ) > 2 ||
           YprojectileNucleon < YtargetNucleon ) {        
        OuterSuccess = false;
        continue;
      } 

    } while ( ( ! OuterSuccess ) &&
              ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */
    if ( loopCounter >= maxNumberOfLoops ) {
      #ifdef debugAdjust
      G4cout << "BAD situation: forced exit of the while loop!" << G4endl;
      #endif
      return false;
    }

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
      //G4cout << "primaries[ahadron] " << primaries[ahadron] << G4endl;
      //if ( primaries[ahadron]->GetStatus() <= 1 ) isProjectile=true;
      FirstString = 0; SecondString = 0;
      if ( primaries[ahadron]->GetStatus() == 0 )                                    // Uzhi May 2016
      {
       theExcitation->CreateStrings( primaries[ ahadron ], isProjectile, 
                                     FirstString, SecondString, theParameters );
      }
      else if ( primaries[ahadron]->GetStatus() == 1  
             && primaries[ahadron]->GetSoftCollisionCount() != 0 )                   // Uzhi May 2016
      {
       theExcitation->CreateStrings( primaries[ ahadron ], isProjectile, 
                                     FirstString, SecondString, theParameters );
      }
      else if ( primaries[ahadron]->GetStatus() == 1  
             && primaries[ahadron]->GetSoftCollisionCount() == 0 )                   // Uzhi May 2016
      {
       G4LorentzVector ParticleMomentum=primaries[ahadron]->Get4Momentum();
       G4KineticTrack* aTrack=new G4KineticTrack(
                                  primaries[ahadron]->GetDefinition(),
                                  primaries[ahadron]->GetTimeOfCreation(),
                                  primaries[ahadron]->GetPosition(),
                                  ParticleMomentum);
       FirstString=new G4ExcitedString(aTrack);
      }
      else if(primaries[ahadron]->GetStatus() == 2)
      {
       G4LorentzVector ParticleMomentum=primaries[ahadron]->Get4Momentum();
       G4KineticTrack* aTrack=new G4KineticTrack(
                                  primaries[ahadron]->GetDefinition(),
                                  primaries[ahadron]->GetTimeOfCreation(),
                                  primaries[ahadron]->GetPosition(),
                                  ParticleMomentum);
       FirstString=new G4ExcitedString(aTrack);
      }
      else {G4cout<<"Something wrong in FTF Model Build String" << G4endl;}

      if ( FirstString  != 0 ) strings->push_back( FirstString );
      if ( SecondString != 0 ) strings->push_back( SecondString );

      #ifdef debugBuildString
      G4cout << "FirstString & SecondString? " << FirstString << " " << SecondString << G4endl;
      if(FirstString->IsExcited())
      {
       G4cout<< "Quarks on the FirstString ends " << FirstString->GetRightParton()->GetPDGcode()
             << " " << FirstString->GetLeftParton()->GetPDGcode() << G4endl;
      } else {G4cout<<"Kinetic track is stored"<<G4endl;}
      #endif

    }

    #ifdef debugBuildString
    if(FirstString->IsExcited())
    {
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
      if ( aProjectile->GetStatus() == 0 ) { // A nucleon took part in non-diffractive interaction

        #ifdef debugBuildString
        G4cout << "Case1 aProjectile->GetStatus() == 0 " << G4endl;
        #endif

        theExcitation->CreateStrings( 
                           TheInvolvedNucleonsOfProjectile[ ahadron ]->GetSplitableHadron(),
                           isProjectile, FirstString, SecondString, theParameters );

      } else if ( aProjectile->GetStatus() == 1 && aProjectile->GetSoftCollisionCount() != 0 ) { 
        // Nucleon took part in diffractive interaction

        #ifdef debugBuildString
        G4cout << "Case2 aProjectile->GetStatus() !=0 St==1 SoftCol!=0" << G4endl;
        #endif

        theExcitation->CreateStrings( 
                           TheInvolvedNucleonsOfProjectile[ ahadron ]->GetSplitableHadron(),
                           isProjectile, FirstString, SecondString, theParameters );

      } else if ( aProjectile->GetStatus() == 1  &&  aProjectile->GetSoftCollisionCount() == 0  &&
                  HighEnergyInter ) {
        // Nucleon was considered as a paricipant of an interaction,
        // but the interaction was skipped due to annihilation.
        // It is now considered as an involved nucleon at high energies.

        #ifdef debugBuildString
        G4cout << "Case3 aProjectile->GetStatus() !=0 St==1 SoftCol==0" << G4endl;
        #endif

        G4LorentzVector ParticleMomentum=aProjectile->Get4Momentum();
        G4KineticTrack* aTrack=new G4KineticTrack(
                                   aProjectile->GetDefinition(),
                                   aProjectile->GetTimeOfCreation(),
                                   aProjectile->GetPosition(),
                                   ParticleMomentum);
        FirstString=new G4ExcitedString(aTrack);

        #ifdef debugBuildString
        G4cout << " Strings are built for nucleon marked for an interaction, but"
               << " the interaction was skipped." << G4endl;
        #endif

      } else if ( (aProjectile->GetStatus() == 2) || (aProjectile->GetStatus() == 3) ) {    // Uzhi Nov. 2014
        // Nucleon which was involved in the Reggeon cascading

        #ifdef debugBuildString
        G4cout << "Case4 aProjectile->GetStatus() !=0 St==2 " << G4endl;
        #endif

        G4LorentzVector ParticleMomentum=aProjectile->Get4Momentum();
        G4KineticTrack* aTrack=new G4KineticTrack(
                                   aProjectile->GetDefinition(),
                                   aProjectile->GetTimeOfCreation(),
                                   aProjectile->GetPosition(),
                                   ParticleMomentum);
        FirstString=new G4ExcitedString(aTrack);

        #ifdef debugBuildString
        G4cout << " A track is build for involved nucleon." << G4endl;
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

      if ( FirstString  != 0 ) strings->push_back( FirstString );
      if ( SecondString != 0 ) strings->push_back( SecondString );
    } // end of for ( G4int ahadron = 0; ahadron < NumberOfInvolvedNucleonsOfProjectile
  }   // ens of if ( ! GetProjectileNucleus() )

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

    if ( aNucleon->GetStatus() == 0 ) { 
      // A nucleon took part in non-diffractive interaction
      theExcitation->CreateStrings( aNucleon, isProjectile, 
                                    FirstString, SecondString, theParameters );
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

    } else if ( aNucleon->GetStatus() == 1  &&  aNucleon->GetSoftCollisionCount() == 0  &&
                HighEnergyInter ) {
      // A nucleon was considered as a participant but due to annihilation
      // its interactions were skipped. It will be considered as involved one
      // at high energies.

      G4LorentzVector ParticleMomentum=aNucleon->Get4Momentum();
      G4KineticTrack* aTrack=new G4KineticTrack(
                                  aNucleon->GetDefinition(),
                                  aNucleon->GetTimeOfCreation(),
                                  aNucleon->GetPosition(),
                                  ParticleMomentum);

      FirstString=new G4ExcitedString(aTrack);

      #ifdef debugBuildString
      G4cout << "3 case A track is build" << G4endl;
      #endif

    } else if ( aNucleon->GetStatus() == 1  &&  aNucleon->GetSoftCollisionCount() == 0  &&
                ! HighEnergyInter ) {
      // A nucleon was considered as a participant but due to annihilation
      // its interactions were skipped. It will be returned to nucleus
      // at low energies energies.
      aNucleon->SetStatus( 5 );
      // ????????? delete aNucleon;

      #ifdef debugBuildString
      G4cout << "4 case A string is not build" << G4endl;
      #endif

    } else if(( aNucleon->GetStatus() == 2 )||   // A nucleon took part in quark exchange 
              ( aNucleon->GetStatus() == 3 )  ){ // A nucleon was involved in Reggeon cascading


      G4LorentzVector ParticleMomentum=aNucleon->Get4Momentum();
      G4KineticTrack* aTrack=new G4KineticTrack(
                                 aNucleon->GetDefinition(),
                                 aNucleon->GetTimeOfCreation(),
                                 aNucleon->GetPosition(), //FirstString->GetPosition(),
                                 ParticleMomentum);

      FirstString=new G4ExcitedString(aTrack);

      #ifdef debugBuildString
      G4cout << "5 case A track is build" << G4endl;
      #endif

    } else {

      #ifdef debugBuildString
      G4cout << "6 case No string" << G4endl;
      #endif

    }

    if ( FirstString  != 0 ) strings->push_back( FirstString );
    if ( SecondString != 0 ) strings->push_back( SecondString );

  }   // end of for ( G4int ahadron = 0; ahadron < NumberOfInvolvedNucleonsOfTarget

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

//-------------------------------------
    if( TargetResidualMassNumber != 0 )
    {
     G4ThreeVector bstToCM =TargetResidual4Momentum.findBoostToCM();

     G4V3DNucleus* theTargetNucleus = GetTargetNucleus();
     G4LorentzVector residualMomentum(0.,0.,0.,0.);
     G4Nucleon* aNucleon = 0;
     theTargetNucleus->StartLoop();
     while ( ( aNucleon = theTargetNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
       if ( !aNucleon->AreYouHit() ) { 
         G4LorentzVector tmp=aNucleon->Get4Momentum(); tmp.boost(bstToCM);
         aNucleon->SetMomentum(tmp);
         residualMomentum +=tmp;
       }
     }

     residualMomentum/=TargetResidualMassNumber;

     G4double Mass = TargetResidual4Momentum.mag();
     G4double SumMasses=0.;
  
     aNucleon = 0;
     theTargetNucleus->StartLoop();
     while ( ( aNucleon = theTargetNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
       if ( !aNucleon->AreYouHit() ) { 
         G4LorentzVector tmp=aNucleon->Get4Momentum() - residualMomentum;
         G4double E=std::sqrt(tmp.vect().mag2()+
                              sqr(aNucleon->GetDefinition()->GetPDGMass()-aNucleon->GetBindingEnergy()));
         tmp.setE(E);  aNucleon->SetMomentum(tmp);
         SumMasses+=E;
       }
     }

     G4double Chigh=Mass/SumMasses; G4double Clow=0; G4double C;
     const G4int maxNumberOfLoops = 1000;
     G4int loopCounter = 0;
     do
     {
      C=(Chigh+Clow)/2.;

      SumMasses=0.;
      aNucleon = 0;
      theTargetNucleus->StartLoop();
      while ( ( aNucleon = theTargetNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
        if ( !aNucleon->AreYouHit() ) { 
         G4LorentzVector tmp=aNucleon->Get4Momentum();
         G4double E=std::sqrt(tmp.vect().mag2()*sqr(C)+
                              sqr(aNucleon->GetDefinition()->GetPDGMass()-aNucleon->GetBindingEnergy()));
         SumMasses+=E;
        }
      }

      if(SumMasses > Mass) {Chigh=C;}
      else                 {Clow =C;}

     } while( (Chigh-Clow > 0.01) &&  // end do
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
        G4LorentzVector tmp=aNucleon->Get4Momentum()*C;
        G4double E=std::sqrt(tmp.vect().mag2()+
                             sqr(aNucleon->GetDefinition()->GetPDGMass()-aNucleon->GetBindingEnergy()));
        tmp.setE(E); tmp.boost(-bstToCM);  
        aNucleon->SetMomentum(tmp);     
       }
     }
    }   // End of if( TargetResidualMassNumber != 0 )
//-------------------------------------

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

//-------------------------------------
    if( ProjectileResidualMassNumber != 0 )
    {
     G4ThreeVector bstToCM =ProjectileResidual4Momentum.findBoostToCM();

     G4V3DNucleus* theProjectileNucleus = GetProjectileNucleus();
     G4LorentzVector residualMomentum(0.,0.,0.,0.);
     G4Nucleon* aNucleon = 0;
     theProjectileNucleus->StartLoop();
     while ( ( aNucleon = theProjectileNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
       if ( !aNucleon->AreYouHit() ) { 
         G4LorentzVector tmp=aNucleon->Get4Momentum(); tmp.boost(bstToCM);
         aNucleon->SetMomentum(tmp);
         residualMomentum +=tmp;
       }
     }

     residualMomentum/=ProjectileResidualMassNumber;

     G4double Mass = ProjectileResidual4Momentum.mag();
     G4double SumMasses=0.;
  
     aNucleon = 0;
     theProjectileNucleus->StartLoop();
     while ( ( aNucleon = theProjectileNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
       if ( !aNucleon->AreYouHit() ) { 
         G4LorentzVector tmp=aNucleon->Get4Momentum() - residualMomentum;
         G4double E=std::sqrt(tmp.vect().mag2()+
                              sqr(aNucleon->GetDefinition()->GetPDGMass()-aNucleon->GetBindingEnergy()));
         tmp.setE(E);  aNucleon->SetMomentum(tmp);
         SumMasses+=E;
       }
     }

     G4double Chigh=Mass/SumMasses; G4double Clow=0; G4double C;
     const G4int maxNumberOfLoops = 1000;
     G4int loopCounter = 0;
     do
     {
      C=(Chigh+Clow)/2.;

      SumMasses=0.;
      aNucleon = 0;
      theProjectileNucleus->StartLoop();
      while ( ( aNucleon = theProjectileNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
        if ( !aNucleon->AreYouHit() ) { 
         G4LorentzVector tmp=aNucleon->Get4Momentum();
         G4double E=std::sqrt(tmp.vect().mag2()*sqr(C)+
                              sqr(aNucleon->GetDefinition()->GetPDGMass()-aNucleon->GetBindingEnergy()));
         SumMasses+=E;
        }
      }

      if(SumMasses > Mass) {Chigh=C;}
      else                 {Clow =C;}

     } while( (Chigh-Clow > 0.01) &&  // end do
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
       if ( !aNucleon->AreYouHit() ) { 
        G4LorentzVector tmp=aNucleon->Get4Momentum()*C;
        G4double E=std::sqrt(tmp.vect().mag2()+
                             sqr(aNucleon->GetDefinition()->GetPDGMass()-aNucleon->GetBindingEnergy()));
        tmp.setE(E); tmp.boost(-bstToCM);  
        aNucleon->SetMomentum(tmp);     
       }
     }
    }   // End of if( ProjectileResidualMassNumber != 0 )
//-------------------------------------  
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
/*                                                                                Closed by Uzhi, May 2016
G4ThreeVector G4FTFModel::GaussianPt( G4double AveragePt2, G4double maxPtSquare ) const {
  //  @@ this method is used in FTFModel as well. Should go somewhere common!
G4cout<<"MaxPt^2 "<<maxPtSquare<<" "<<G4endl;
G4cout<<"Enter AveragePt2"<<G4endl; 
{do
 {
  G4cin>>AveragePt2;
  G4cout<<AveragePt2<<" "<<maxPtSquare<<" "<<G4endl;
  G4cout<<"Argument and G4Exp( -maxPtSquare/AveragePt2 ) "<<maxPtSquare/AveragePt2<<" "<<G4Exp( -maxPtSquare/AveragePt2)<<G4endl;
 } while(true);
}
  G4double Pt2( 0.0 );

  if ( AveragePt2 <= 0.0 ) {
    Pt2 = 0.0;
  } else {
    Pt2 = -AveragePt2 * G4Log( 1.0 + G4UniformRand() * 
                                        ( G4Exp( -maxPtSquare/AveragePt2 ) -1.0 ) );
  }

  G4double Pt = std::sqrt( Pt2 );
  G4double phi = G4UniformRand() * twopi;

  return G4ThreeVector( Pt*std::cos(phi), Pt*std::sin(phi), 0.0 );    
}
*/
//
//============================================================================ Uzhi 2016

G4ThreeVector G4FTFModel::GaussianPt( G4double AveragePt2, G4double maxPtSquare ) const {

  G4double Pt2( 0.0 );

  if(AveragePt2 > 0.0) {
    if(maxPtSquare/AveragePt2 < 1.0e+9) {
      Pt2 = -AveragePt2 * G4Log( 1.0 + G4UniformRand() * 
                                     ( G4Exp( -maxPtSquare/AveragePt2 ) -1.0 ) );
    } else {
      Pt2 = -AveragePt2 * G4Log( 1.0 - G4UniformRand() ); 
    }
  }

  G4double Pt = std::sqrt( Pt2 );
  G4double phi = G4UniformRand() * twopi;
 
  return G4ThreeVector( Pt*std::cos(phi), Pt*std::sin(phi), 0.0 );    
}
//

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

      residualExcitationEnergy += -ExcitationEnergyPerWoundedNucleon*
                                   G4Log( G4UniformRand());
      residualMassNumber--;
      // The absolute value below is needed only in the case of anti-nucleus.
      residualCharge -= std::abs( G4int( aNucleon->GetDefinition()->GetPDGCharge() ) );
    } else {   // Spectator nucleons
      residualMomentum += aNucleon->Get4Momentum();
    }
  }
  #ifdef debugPutOnMassShell
  G4cout << "ExcitationEnergyPerWoundedNucleon " << ExcitationEnergyPerWoundedNucleon << G4endl
         << "\t Residual Charge, MassNumber " << residualCharge << " " << residualMassNumber
         << G4endl << "\t Initial Momentum " << nucleusMomentum
         << G4endl << "\t Residual Momentum   " << residualMomentum << G4endl;
  #endif
  residualMomentum.setPz( 0.0 ); 
  residualMomentum.setE( 0.0 );
  if ( residualMassNumber == 0 ) {
    residualMass = 0.0;
    residualExcitationEnergy = 0.0;
  } else {
    residualMass = G4ParticleTable::GetParticleTable()->GetIonTable()->
                     GetIonMass( residualCharge, residualMassNumber );
    if ( residualMassNumber == 1 ) {
      residualExcitationEnergy = 0.0;
    }
    residualMass += residualExcitationEnergy;       // Uzhi March 2016 ????
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

  //const G4double ProbDeltaIsobar = 0.05;
  //const G4double ProbDeltaIsobar = 0.25;
  const G4double probDeltaIsobar = 0.05;  // A.R. 07.08.2013 0.10 -> 0.05 Uzhi March 2016

  G4int maxNumberOfDeltas = G4int( (sqrtS - sumMasses)/(400.0*MeV) );
  G4int numberOfDeltas = 0;

  for ( G4int i = 0; i < numberOfInvolvedNucleons; i++ ) {
    //G4cout << "i maxNumberOfDeltas probDeltaIsobar " << i << " " << maxNumberOfDeltas
    //       << " " << probDeltaIsobar << G4endl;
    if ( G4UniformRand() < probDeltaIsobar  &&  numberOfDeltas < maxNumberOfDeltas ) {
      numberOfDeltas++;
      if ( ! involvedNucleons[i] ) continue;
      G4VSplitableHadron* splitableHadron = involvedNucleons[i]->GetSplitableHadron();
      G4double massNuc = std::sqrt( sqr( splitableHadron->GetDefinition()->GetPDGMass() )
                                    + splitableHadron->Get4Momentum().perp2() );
      //AR The absolute value below is needed in the case of an antinucleus. 
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
        sumMasses += ( massDelta - massNuc );        // Uzhi March 2016 ???
      }
    } 
  }
  //G4cout << "maxNumberOfDeltas numberOfDeltas " << maxNumberOfDeltas << " " 
  //       << numberOfDeltas << G4endl;
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

  if ( ! nucleus ) return false;

  if ( residualMassNumber == 0  &&  numberOfInvolvedNucleons == 1 ) {
    dCor = 0.0; 
    averagePt2 = 0.0;
  } 

  G4bool success = true;                            
  G4double SumMasses = residualMass; 
//                                                         // Uzhi March 2016 ???
  for ( G4int i = 0; i < numberOfInvolvedNucleons; i++ ) {
    G4Nucleon* aNucleon = involvedNucleons[i];
    if ( ! aNucleon ) continue;
    SumMasses += aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass();
  }
//
  const G4int maxNumberOfLoops = 1000;
  G4int loopCounter = 0;
  do {  // while ( ! success )

    success = true;
//======================================= Sampling of nucleon Pt ===============
    G4ThreeVector ptSum( 0.0, 0.0, 0.0 );

    for ( G4int i = 0; i < numberOfInvolvedNucleons; i++ ) {
      G4Nucleon* aNucleon = involvedNucleons[i];
      if ( ! aNucleon ) continue;
      G4ThreeVector tmpPt = GaussianPt( averagePt2, maxPt2 );
      ptSum += tmpPt;

      G4LorentzVector tmp( tmpPt.x(), tmpPt.y(), 0., 0.);
      aNucleon->SetMomentum( tmp );
    }

    G4double deltaPx = ( ptSum.x() - pResidual.x() ) / numberOfInvolvedNucleons;
    G4double deltaPy = ( ptSum.y() - pResidual.y() ) / numberOfInvolvedNucleons;


    SumMasses = residualMass;
    for ( G4int i = 0; i < numberOfInvolvedNucleons; i++ ) {
      G4Nucleon* aNucleon = involvedNucleons[i];
      if ( ! aNucleon ) continue;
      G4double px = aNucleon->Get4Momentum().px() - deltaPx;
      G4double py = aNucleon->Get4Momentum().py() - deltaPy;
      G4double MtN = std::sqrt( sqr( aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass() )
                              + sqr( px ) + sqr( py ) );
      SumMasses += MtN;
      G4LorentzVector tmp( px, py, 0., MtN);
      aNucleon->SetMomentum( tmp );
    }
//======================================== Sampling X of nucleon ===============
    G4double xSum = 0.0;

    for ( G4int i = 0; i < numberOfInvolvedNucleons; i++ ) {
      G4Nucleon* aNucleon = involvedNucleons[i];
      if ( ! aNucleon ) continue;

      G4ThreeVector tmpX = GaussianPt( dCor*dCor, 1.0 );
//      G4double x = tmpX.x() +                                 // Uzhi 2016
//                   aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass()/SumMasses;
      G4double x = tmpX.x() + aNucleon->Get4Momentum().e()/SumMasses;
      if ( x < 0.0  ||  x > 1.0 ) { 
        success = false; 
        break;
      }
      xSum += x;
      //AR The energy is in the lab (instead of cms) frame but it will not be used.
      G4LorentzVector tmp( aNucleon->Get4Momentum().x(), aNucleon->Get4Momentum().y(), 
                                    x, aNucleon->Get4Momentum().e() );
      aNucleon->SetMomentum( tmp );
    }

    if ( xSum < 0.0  ||  xSum > 1.0 ) success = false;

    if ( ! success ) continue;

//    G4double deltaPx = ( ptSum.x() - pResidual.x() ) / numberOfInvolvedNucleons;  // Uzhi 2016
//    G4double deltaPy = ( ptSum.y() - pResidual.y() ) / numberOfInvolvedNucleons;
    G4double delta = 0.0;

    if ( residualMassNumber == 0 ) {
      delta = ( xSum - 1.0 ) / numberOfInvolvedNucleons;
    } else {
      delta = 0.0;
    }

    xSum = 1.0;
    mass2 = 0.0;

    for ( G4int i = 0; i < numberOfInvolvedNucleons; i++ ) {
      G4Nucleon* aNucleon = involvedNucleons[i];
      if ( ! aNucleon ) continue;
      G4double x = aNucleon->Get4Momentum().pz() - delta;
      xSum -= x;               
      if ( residualMassNumber == 0 ) {
        if ( x <= 0.0  ||  x > 1.0 ) {
          success = false; 
          break;
        }
      } else {
        if ( x <= 0.0  ||  x > 1.0  ||  xSum <= 0.0  ||  xSum > 1.0 ) {
          success = false; 
          break;
        }
      }                                          
/*                                                            // Uzhi 2016
      G4double px = aNucleon->Get4Momentum().px() - deltaPx;
      G4double py = aNucleon->Get4Momentum().py() - deltaPy;
      mass2 += ( sqr( aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass() )
                    + sqr( px ) + sqr( py ) ) / x;
      G4LorentzVector tmp( px, py, x, aNucleon->Get4Momentum().e() );
*/
      mass2 += sqr( aNucleon->Get4Momentum().e() ) / x;
      G4LorentzVector tmp( aNucleon->Get4Momentum().px(), aNucleon->Get4Momentum().py(), 
                                 x, aNucleon->Get4Momentum().e() );
      aNucleon->SetMomentum( tmp );
    }

    if ( ! success ) continue;
//=======================================================

    if ( success  &&  residualMassNumber != 0 ) {
      mass2 += ( sqr( residualMass ) + pResidual.perp2() ) / xSum;      // Uzhi 2016
//      mass2 += sqr( residualMass ) / xSum;
    }

    #ifdef debugPutOnMassShell
    G4cout << "success " << success << G4endl << " Mt " << std::sqrt( mass2 )/GeV << G4endl;
    #endif

  } while ( ( ! success ) &&
            ++loopCounter < maxNumberOfLoops );
  if ( loopCounter >= maxNumberOfLoops ) {
    return false;
  }

  return true;
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

//G4cout<<"sqrtS Mp Mt sum "<<std::sqrt(sValue)<<"  "<<std::sqrt(projectileMass2)<<" "<<std::sqrt(targetMass2)<<" "<<std::sqrt(projectileMass2)+std::sqrt(targetMass2)<<G4endl;

  G4double decayMomentum2 = sqr( sValue ) + sqr( projectileMass2 ) + sqr( targetMass2 )
                            - 2.0*sValue*projectileMass2 - 2.0*sValue*targetMass2 
                            - 2.0*projectileMass2*targetMass2;
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
  #endif

  for ( G4int i = 0; i < numberOfInvolvedNucleons; i++ ) {
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
    G4cout << "i nY pY nY-AY AY " << i << " " << nucleonY << " " << projectileY <<G4endl;
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

  for ( G4int i = 0; i < numberOfInvolvedNucleons; i++ ) {
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
  G4cout << "w residual3Momentum.z() " << w << " " << residual3Momentum.z() << G4endl;
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
