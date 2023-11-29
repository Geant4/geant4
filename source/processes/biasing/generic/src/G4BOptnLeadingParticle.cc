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
#include "G4BOptnLeadingParticle.hh"
#include "G4BiasingProcessInterface.hh"

#include <vector>
#include <map>


G4BOptnLeadingParticle::G4BOptnLeadingParticle(G4String name)
  : G4VBiasingOperation                ( name ),
    fRussianRouletteKillingProbability ( -1.0 )
{
}

G4BOptnLeadingParticle::~G4BOptnLeadingParticle()
{
}

G4VParticleChange* G4BOptnLeadingParticle::ApplyFinalStateBiasing( const G4BiasingProcessInterface* callingProcess,
								   const G4Track*                            track,
								   const G4Step*                               step,
								   G4bool&                                           )
{
  // -- collect wrapped process particle change:
  auto wrappedProcessParticleChange = callingProcess->GetWrappedProcess()->PostStepDoIt(*track,*step);

  // -- does nothing in case the primary stays alone or in weird situation where all are killed...
  if ( wrappedProcessParticleChange->GetNumberOfSecondaries() == 0 )                        return wrappedProcessParticleChange;
  if ( wrappedProcessParticleChange->GetTrackStatus()         == fKillTrackAndSecondaries ) return wrappedProcessParticleChange;
  // -- ... else, collect the secondaries in a same vector (the primary is pushed in this vector, if surviving, later see [**]):
  std::vector < G4Track* > secondariesAndPrimary;
  for ( auto i = 0 ; i < wrappedProcessParticleChange->GetNumberOfSecondaries() ; i++ ) secondariesAndPrimary.push_back(  wrappedProcessParticleChange->GetSecondary( i ) );
  

  // -- If case the primary survives, need to collect its new state. In the general case of the base class G4VParticleChange
  // -- this is tricky, as this class does not hold the primary changes (and we have to build a fake step and fake track
  // -- caring about the touchables, etc.) So we limit here to the G4ParticleChange case, checking the reality of this
  // -- class with a dynamic cast. If we don't have such an actual G4DynamicParticle object, we give up the biasing and return
  // -- the untrimmed process final state.
  // -- Note that this case does not happen often, as the technique is intended for inelastic processes. For case where several
  // -- particles can be produced without killing the primary, we have for example the electron-/positron-nuclear
  G4ParticleChange* castParticleChange ( nullptr );
  G4Track*          finalStatePrimary  ( nullptr );
  if ( ( wrappedProcessParticleChange->GetTrackStatus() != fStopAndKill ) )
    {
      // fFakePrimaryTrack->CopyTrackInfo( *track );
      // fFakeStep        ->InitializeStep( fFakePrimaryTrack );
      // wrappedProcessParticleChange->UpdateStepForPostStep( fFakeStep );
      // fFakeStep->UpdateTrack();
      castParticleChange = dynamic_cast< G4ParticleChange* >(  wrappedProcessParticleChange );
      if ( castParticleChange == nullptr )
	{
	  G4cout << " **** G4BOptnLeadingParticle::ApplyFinalStateBiasing(...) : can not bias for " << callingProcess->GetProcessName() << ", this is just a warning." << G4endl;
	  return wrappedProcessParticleChange;
	}
      finalStatePrimary = new G4Track( *track );
      finalStatePrimary->SetKineticEnergy    ( castParticleChange->GetEnergy()               );
      finalStatePrimary->SetWeight           ( castParticleChange->GetWeight()               );
      finalStatePrimary->SetMomentumDirection( *(castParticleChange->GetMomentumDirection()) );
      // -- [**] push the primary as the last track in the vector of tracks:
      secondariesAndPrimary.push_back( finalStatePrimary );
    }
  
  // -- Ensure the secondaries all have the primary weight:
  // ---- collect primary track weight, from updated by process if alive, or from original copy if died:
  G4double                 primaryWeight;
  if ( finalStatePrimary ) primaryWeight = finalStatePrimary->GetWeight();
  else                     primaryWeight = track            ->GetWeight();
  // ---- now set this same weight to all secondaries:
  for ( auto i = 0 ; i < wrappedProcessParticleChange->GetNumberOfSecondaries() ; i++ ) secondariesAndPrimary[ i ]->SetWeight( primaryWeight );
  

  // -- finds the leading particle, initialize a map of surviving tracks, tag as surviving the leading track:
  size_t   leadingIDX    =  0;
  G4double leadingEnergy = -1;
  std::map< G4Track*, G4bool > survivingMap;
  for ( size_t idx = 0; idx < secondariesAndPrimary.size(); idx++ )
    {
      survivingMap[ secondariesAndPrimary[idx] ] = false;
      if ( secondariesAndPrimary[idx]->GetKineticEnergy() > leadingEnergy )
	{
	  leadingEnergy = secondariesAndPrimary[idx]->GetKineticEnergy();
	  leadingIDX    = idx;
	}
    }
  survivingMap[ secondariesAndPrimary[leadingIDX] ] = true; // -- tag as surviving the leading particle
  
  
  // -- now make track vectors of given types ( choose type = abs(PDG) ), excluding the leading particle:
  std::map < G4int, std::vector< G4Track* > > typesAndTracks;
  for ( size_t idx = 0; idx < secondariesAndPrimary.size(); idx++ )
    {
      if ( idx == leadingIDX ) continue; // -- excludes the leading particle
      auto currentTrack = secondariesAndPrimary[idx];
      auto GROUP        = std::abs( currentTrack->GetDefinition()->GetPDGEncoding() ); // -- merge particles and anti-particles in the same category  -- §§ this might be proposed as an option in future
      if ( currentTrack->GetDefinition()->GetBaryonNumber() >= 2 ) GROUP = -1000; // -- merge all baryons above proton/neutron in one same group -- §§ might be proposed as an option too
      
      if ( typesAndTracks.find( GROUP ) == typesAndTracks.end() )
	{
	  std::vector< G4Track* > v;
	  v.push_back( currentTrack );
	  typesAndTracks[ GROUP ] = v;
	}
      else
	{
	  typesAndTracks[ GROUP ].push_back( currentTrack );
	}
    }
  // -- and on these vectors, randomly select the surviving particles:
  // ---- randomly select one surviving track per species
  // ---- for this surviving track, further apply a Russian roulette
  G4int nSecondaries = 0; // -- the number of secondaries to be returned
  for ( auto& typeAndTrack : typesAndTracks )
    {
      size_t nTracks = (typeAndTrack.second).size();
      G4Track* keptTrack;
      // -- select one track among ones in same species:
      if ( nTracks > 1 )
	{
	  auto keptTrackIDX = G4RandFlat::shootInt( nTracks );
	  keptTrack = (typeAndTrack.second)[keptTrackIDX];
	  keptTrack->SetWeight( keptTrack->GetWeight() * nTracks );
	}
      else
	{
	  keptTrack = (typeAndTrack.second)[0];
	}
      // -- further apply a Russian Roulette on it:
      G4bool keepTrack = false;
      if ( fRussianRouletteKillingProbability > 0.0 )
	{
	  if ( G4UniformRand() > fRussianRouletteKillingProbability )
	    {
	      keptTrack->SetWeight( keptTrack->GetWeight() / (1. - fRussianRouletteKillingProbability) );
	      keepTrack = true;
	    }
	}
      else keepTrack = true;
      if ( keepTrack )
	{
	  survivingMap[ keptTrack ] = true;
	  if ( keptTrack != finalStatePrimary ) nSecondaries++;
	}
    }
  // -- and if the leading is not the primary, we have to count it in nSecondaries:
  if ( secondariesAndPrimary[leadingIDX] != finalStatePrimary ) nSecondaries++;
  
  // -- verify if the primary is still alive or not after above selection:
  G4bool primarySurvived = false;
  if ( finalStatePrimary ) primarySurvived = survivingMap[ finalStatePrimary ];
    
  
  // -- fill the trimmed particle change:
  // ---- fill for the primary:
  fParticleChange.Initialize(*track);
  if ( primarySurvived )
    {
      fParticleChange.ProposeTrackStatus       ( wrappedProcessParticleChange->GetTrackStatus()   );
      fParticleChange.ProposeParentWeight      ( finalStatePrimary->GetWeight()                   ); // -- take weight from copy of primary, this one being updated in the random selection loop above
      fParticleChange.ProposeEnergy            ( finalStatePrimary->GetKineticEnergy()            );
      fParticleChange.ProposeMomentumDirection ( finalStatePrimary->GetMomentumDirection()        );
    }
  else
    {
      fParticleChange.ProposeTrackStatus ( fStopAndKill );
      fParticleChange.ProposeParentWeight( 0.0 );
      fParticleChange.ProposeEnergy      ( 0.0 );
    }
  // -- fill for surviving secondaries:
  fParticleChange.SetSecondaryWeightByProcess(true);
  fParticleChange.SetNumberOfSecondaries(nSecondaries);
  // ---- note we loop up to on the number of secondaries, which excludes the primary, last in secondariesAndPrimary vector:
  //////  G4cout << callingProcess->GetProcessName() << " :";
  for ( auto idx = 0 ; idx < wrappedProcessParticleChange->GetNumberOfSecondaries() ; idx++ )
    {
      G4Track* secondary = secondariesAndPrimary[idx];
      // ********************
      //// if ( !survivingMap[ secondary ] ) G4cout << " [";
      ///// else G4cout << " ";
      ///// G4cout << secondary->GetDefinition()->GetParticleName() << " " << secondary->GetKineticEnergy();
      ///// if ( !survivingMap[ secondary ] ) G4cout << "]";
      //// if ( secondary ==  secondariesAndPrimary[leadingIDX] ) G4cout << " ***";
      // ******************
      if ( survivingMap[ secondary ] )  fParticleChange.AddSecondary( secondary );
      else                              delete secondary;
    }
  ///  G4cout << G4endl;
  
  // -- clean the wrapped process particle change:
  wrappedProcessParticleChange->Clear();

  if ( finalStatePrimary ) delete finalStatePrimary;  

  // -- finally, returns the trimmed particle change:
  return &fParticleChange;
  
}
