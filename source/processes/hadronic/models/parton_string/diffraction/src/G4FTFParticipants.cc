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
// $Id: G4FTFParticipants.cc 106965 2017-10-31 08:40:14Z gcosmo $
// GEANT4 tag $Name:  $
//

// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4FTFParticipants----------------
//             by Gunter Folger, June 1998.
//       class finding colliding particles in FTFPartonStringModel
//  Changed in a part by V. Uzhinsky in oder to put in correcpondence
//        with original FRITIOF mode. November - December 2006.
//  Ajusted for (anti) nucleus - nucleus interactions by V. Uzhinsky.
//                    (February 2011)
// ------------------------------------------------------------

#include <utility>
#include <vector>
#include <algorithm>

#include "G4FTFParticipants.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4FTFParameters.hh"                            
#include "G4DiffractiveSplitableHadron.hh"
#include "G4VSplitableHadron.hh"


//============================================================================

//#define debugFTFparticipant


//============================================================================

G4FTFParticipants::G4FTFParticipants() : currentInteraction( -1 ) {}


//============================================================================

G4FTFParticipants::G4FTFParticipants( const G4FTFParticipants& ) : 
  G4VParticipants(), currentInteraction( -1 ) 
{
  G4Exception( "G4FTFParticipants::G4FTFParticipants()", "HAD_FTF_001",
               FatalException, " Must not use copy ctor()" );
}


//============================================================================

G4FTFParticipants::~G4FTFParticipants() {}


//============================================================================

void G4FTFParticipants::GetList( const G4ReactionProduct& thePrimary, 
                                 G4FTFParameters* theParameters ) { 

  #ifdef debugFTFparticipant
  G4cout << "Participants::GetList" << G4endl 
         << "thePrimary " << thePrimary.GetMomentum() << G4endl << G4endl;
  #endif

  G4double betta_z = thePrimary.GetMomentum().z() / thePrimary.GetTotalEnergy();
  if ( betta_z < 1.0e-10 ) betta_z = 1.0e-10;

  StartLoop(); // reset Loop over Interactions

  for ( unsigned int i = 0; i < theInteractions.size(); i++ ) delete theInteractions[i];
  theInteractions.clear();

  G4double deltaxy = 2.0 * fermi;  // Extra nuclear radius

  if ( theProjectileNucleus == 0 ) {  // Hadron-nucleus or anti-baryon-nucleus interactions

    G4double impactX( 0.0 ), impactY( 0.0 );

    G4VSplitableHadron* primarySplitable = new G4DiffractiveSplitableHadron( thePrimary );

    #ifdef debugFTFparticipant
    G4cout << "Hadron-nucleus or anti-baryon-nucleus interactions" << G4endl;
    #endif

    G4double xyradius;                          
    xyradius = theNucleus->GetOuterRadius() + deltaxy; // Range of impact parameter sampling

    const G4int maxNumberOfLoops = 1000;
    G4int loopCounter = 0;                                          
    do {

      std::pair< G4double, G4double > theImpactParameter;
      theImpactParameter = theNucleus->ChooseImpactXandY( xyradius );
      impactX = theImpactParameter.first; 
      impactY = theImpactParameter.second;

      #ifdef debugFTFparticipant
      G4cout << "New interaction list," << " b= " 
             << std::sqrt( sqr(impactX ) + sqr( impactY ) )/fermi << G4endl;
      #endif

      G4ThreeVector thePosition( impactX, impactY, 0.0 );
      primarySplitable->SetPosition( thePosition );

      theNucleus->StartLoop();
      G4Nucleon* nucleon;

      #ifdef debugFTFparticipant
      G4int TrN( 0 );
      #endif

      while ( ( nucleon = theNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */

        G4double impact2 = sqr( impactX - nucleon->GetPosition().x() ) +
                           sqr( impactY - nucleon->GetPosition().y() );

        if ( theParameters->GetProbabilityOfInteraction( impact2/fermi/fermi ) > 
             G4UniformRand() ) {
          primarySplitable->SetStatus( 1 );  // It takes part in the interaction
          G4VSplitableHadron* targetSplitable = 0;
          if ( ! nucleon->AreYouHit() ) {
            targetSplitable = new G4DiffractiveSplitableHadron( *nucleon );
            nucleon->Hit( targetSplitable );
            targetSplitable->SetStatus( 1 ); // It takes part in the interaction

            #ifdef debugFTFparticipant
            G4cout << "Participated nucleons #, " << TrN << " " << "Splitable Pr* Tr* "
                   << primarySplitable << " " << targetSplitable << G4endl;
            #endif

          }
          G4InteractionContent* aInteraction = new G4InteractionContent( primarySplitable );
          G4Nucleon* PrNucleon = 0;
          aInteraction->SetProjectileNucleon( PrNucleon );
          aInteraction->SetTarget( targetSplitable );
          aInteraction->SetTargetNucleon( nucleon );     
          aInteraction->SetStatus( 1 );                  
          aInteraction->SetInteractionTime( ( primarySplitable->GetPosition().z() + 
                                              nucleon->GetPosition().z() ) / betta_z );
          theInteractions.push_back( aInteraction );
        }

        #ifdef debugFTFparticipant
        TrN++;
        #endif

      } 

    } while ( ( theInteractions.size() == 0 ) && 
              ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */
    if ( loopCounter >= maxNumberOfLoops ) {
      #ifdef debugFTFparticipant
      G4cout << "BAD situation: forced exit from the while loop!" << G4endl;
      #endif
      return;
    }

    #ifdef debugFTFparticipant
    G4cout << "Number of Hit nucleons " << theInteractions.size() << "\t Bx " << impactX/fermi
           << "\t By " << impactY/fermi << "\t B " 
           << std::sqrt( sqr( impactX ) + sqr( impactY ) )/fermi << G4endl << G4endl;
    #endif

    //SortInteractionsIncT(); // Not need because nucleons are sorted in increasing z-coordinates.
    ShiftInteractionTime();  // To put correct times and z-coordinates
    return;

  }  // end of if ( theProjectileNucleus == 0 )

  // Projectile and target are nuclei

  #ifdef debugFTFparticipant
  G4cout << "Projectile and target are nuclei" << G4endl;
  #endif

//G4cout<<theProjectileNucleus->GetOuterRadius()/fermi<<" "<<theNucleus->GetOuterRadius()/fermi<<" "<<deltaxy/fermi<<G4endl;

  G4double xyradius;                          
  xyradius = theProjectileNucleus->GetOuterRadius() +  // Range of impact parameter sampling
             theNucleus->GetOuterRadius() + deltaxy;

  G4double impactX( 0.0 ), impactY( 0.0 );

  const G4int maxNumberOfLoops = 1000;
  G4int loopCounter = 0;
  do {

    std::pair< G4double, G4double > theImpactParameter;
    theImpactParameter = theNucleus->ChooseImpactXandY( xyradius );
    impactX = theImpactParameter.first; 
    impactY = theImpactParameter.second;

    #ifdef debugFTFparticipant
    G4cout << "New interaction list, " << "b "
           << std::sqrt( sqr( impactX ) + sqr( impactY ) )/fermi << G4endl;
    #endif

    G4ThreeVector theBeamPosition( impactX, impactY, 0.0 );                    

    theProjectileNucleus->StartLoop();
    G4Nucleon* ProjectileNucleon;

    #ifdef debugFTFparticipant
    G4int PrNuclN( 0 );
    #endif

    while ( ( ProjectileNucleon = theProjectileNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */
 
      G4VSplitableHadron* ProjectileSplitable = 0;
      theNucleus->StartLoop();
      G4Nucleon* TargetNucleon = 0;

      #ifdef debugFTFparticipant
      G4int TrNuclN( 0 );
      #endif

      while ( ( TargetNucleon = theNucleus->GetNextNucleon() ) ) {  /* Loop checking, 10.08.2015, A.Ribon */

        G4double impact2 = sqr( impactX + ProjectileNucleon->GetPosition().x() - 
                                TargetNucleon->GetPosition().x() ) +
                           sqr( impactY + ProjectileNucleon->GetPosition().y() -
                                TargetNucleon->GetPosition().y() );
        G4VSplitableHadron* TargetSplitable = 0;
        if ( theParameters->GetProbabilityOfInteraction( impact2/fermi/fermi ) >
             G4UniformRand() ) {  // An interaction has happend!

          #ifdef debugFTFparticipant
          G4cout << G4endl << "An Interaction has happend" << G4endl << "Proj N mom " << PrNuclN
                 << " " << ProjectileNucleon->Get4Momentum() << "-------------" << G4endl
                 << "Targ N mom " << TrNuclN << " " << TargetNucleon->Get4Momentum() << G4endl
                 << "PrN TrN Z coords " << ProjectileNucleon->GetPosition().z()/fermi 
                 << " " << TargetNucleon->GetPosition().z()/fermi 
                 << " " << ProjectileNucleon->GetPosition().z()/fermi + 
                           TargetNucleon->GetPosition().z()/fermi << G4endl;
          #endif

          if ( ! ProjectileNucleon->AreYouHit() ) { 
            // Projectile nucleon was not involved until now.
            ProjectileSplitable = new G4DiffractiveSplitableHadron( *ProjectileNucleon );
            ProjectileNucleon->Hit( ProjectileSplitable );
            ProjectileSplitable->SetStatus( 1 );  // It takes part in the interaction
          } else {  // Projectile nucleon was involved before.
            ProjectileSplitable = ProjectileNucleon->GetSplitableHadron();
          }

          if ( ! TargetNucleon->AreYouHit() ) {  // Target nucleon was not involved until now
            TargetSplitable = new G4DiffractiveSplitableHadron( *TargetNucleon );
            TargetNucleon->Hit( TargetSplitable );
            TargetSplitable->SetStatus( 1 );   // It takes part in the interaction
          } else {  // Target nucleon was involved before.
            TargetSplitable = TargetNucleon->GetSplitableHadron();
          }

          G4InteractionContent* anInteraction = new G4InteractionContent( ProjectileSplitable );
          anInteraction->SetTarget( TargetSplitable );
          anInteraction->SetProjectileNucleon( ProjectileNucleon );
          anInteraction->SetTargetNucleon( TargetNucleon );
          anInteraction->SetInteractionTime( ( ProjectileNucleon->GetPosition().z() + 
                                               TargetNucleon->GetPosition().z() ) / betta_z );
          anInteraction->SetStatus( 1 );                     

          #ifdef debugFTFparticipant
          G4cout << "Part anInteraction->GetInteractionTime() " 
                 << anInteraction->GetInteractionTime()/fermi << G4endl
                 << "Splitable Pr* Tr* " << ProjectileSplitable << " " 
                 << TargetSplitable << G4endl;
          #endif

          theInteractions.push_back( anInteraction );

        } // End of an Interaction has happend!

        #ifdef debugFTFparticipant
        TrNuclN++;
        #endif

      }  // End of while ( ( TargetNucleon = theNucleus->GetNextNucleon() ) )

      #ifdef debugFTFparticipant
      PrNuclN++;
      #endif

    }  // End of while ( ( ProjectileNucleon = theProjectileNucleus->GetNextNucleon() ) )

    if ( theInteractions.size() != 0 ) theProjectileNucleus->DoTranslation( theBeamPosition );

  } while ( ( theInteractions.size() == 0 ) &&
            ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */
  if ( loopCounter >= maxNumberOfLoops ) {
    #ifdef debugFTFparticipant
    G4cout << "BAD situation: forced exit from the while loop!" << G4endl;
    #endif
    return;
  }

  SortInteractionsIncT();
  ShiftInteractionTime();

  #ifdef debugFTFparticipant
  G4cout << G4endl << "Number of primary collisions " << theInteractions.size() 
         << "\t Bx " << impactX/fermi << "\t By " << impactY/fermi
         << "\t B " << std::sqrt( sqr( impactX ) + sqr( impactY ) )/fermi << G4endl
         << "FTF participant End. #######################" << G4endl << G4endl;
  #endif
  return;
}


//============================================================================

bool G4FTFPartHelperForSortInT( const G4InteractionContent* Int1, 
                                const G4InteractionContent* Int2 ) {
  return Int1->GetInteractionTime() < Int2->GetInteractionTime();
}


//============================================================================

void G4FTFParticipants::SortInteractionsIncT() {  // on increased T 
  if ( theInteractions.size() < 2 ) return;  // Avoid unnecesary work
  std::sort( theInteractions.begin(), theInteractions.end(), G4FTFPartHelperForSortInT ); 
}


//============================================================================

void G4FTFParticipants::ShiftInteractionTime() {
  G4double InitialTime = theInteractions[0]->GetInteractionTime();
  for ( unsigned int i = 1; i < theInteractions.size(); i++ ) {
    G4double InterTime = theInteractions[i]->GetInteractionTime() - InitialTime;
    theInteractions[i]->SetInteractionTime( InterTime );
    G4InteractionContent* aCollision = theInteractions[i];
    G4VSplitableHadron* projectile = aCollision->GetProjectile();
    G4VSplitableHadron* target = aCollision->GetTarget();
    G4ThreeVector prPosition = projectile->GetPosition();
    prPosition.setZ( target->GetPosition().z() );
    projectile->SetPosition( prPosition );
    projectile->SetTimeOfCreation( InterTime );
    target->SetTimeOfCreation( InterTime );
  }
  return;
}


//============================================================================

void G4FTFParticipants::Clean() {
  for ( size_t i = 0; i < theInteractions.size(); i++ ) {
    if ( theInteractions[ i ] ) {
      delete theInteractions[ i ];
      theInteractions[ i ] = 0;
    }
  }
  theInteractions.clear();
  currentInteraction = -1;
}

