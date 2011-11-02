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
//---------------------------------------------------------------------------
//
// ClassName:   G4FTFCaptureAtRest
//
// Author:    Alberto Ribon
//
// Date:      18 October 2011
//
// Modified:  
//            02 November 2011, A. Ribon : migration to the new exceptions.
//
//----------------------------------------------------------------------------
//

#include "G4FTFCaptureAtRest.hh"


#include "G4ParticleDefinition.hh"
#include "G4HadProjectile.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Nucleus.hh"
#include "G4HadFinalState.hh"
#include "G4NucleiProperties.hh"

#include "G4ProcessManager.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4TheoFSGenerator.hh"


G4FTFCaptureAtRest::G4FTFCaptureAtRest( const G4String& processName )
  : G4VRestProcess( processName, fHadronic ) {

  // Create the FTFP final-state model. 
  // (We follow what it is done in the class G4FTFPAntiBarionBuilder
  // except quasi-elastic which is not needed at rest.)

  theMin =   0.0*GeV;
  theMax = 100.0*TeV;
  theModel = new G4TheoFSGenerator( "FTFP" );

  theStringModel = new G4FTFModel;
  theStringDecay = new G4ExcitedStringDecay( theLund = new G4LundStringFragmentation );
  theStringModel->SetFragmentationModel( theStringDecay );

  theCascade = new G4GeneratorPrecompoundInterface; // Not a cascade - goes straight to Preco 
  thePreEquilib = new G4PreCompoundModel( theHandler = new G4ExcitationHandler );
  theCascade->SetDeExcitation( thePreEquilib );  

  theModel->SetHighEnergyGenerator( theStringModel );
  theModel->SetTransport( theCascade );
  theModel->SetMinEnergy( theMin );
  theModel->SetMaxEnergy( 100*TeV );
}


G4FTFCaptureAtRest::~G4FTFCaptureAtRest() {
  delete theCascade;
  delete theStringDecay;
  delete theStringModel;
  delete theModel;
  delete thePreEquilib;
  delete theHandler;
  delete theLund;
}


G4bool G4FTFCaptureAtRest::IsApplicable( const G4ParticleDefinition& particle )  {
  // For the time being, we use Fritiof annihilation at rest only for
  // anti-protons, but it could apply as well for anti-Sigma+ .
  // For the other anti-baryons that Fritiof is able to annihilate on a
  // nucleus, i.e. anti-neutron, anti-Lambda0, anti-Sigma-, anti-Sigma0,
  // anti-Csi-, anti-Csi0, and anti-Omega-, they cannot have "at rest"
  // capture in a nucleus because either they are neutrals and therefore
  // never at rest, or they are positively charged and therefore cannot
  // be captured in a nucleus.
  if ( particle == *( G4AntiProton::AntiProton() ) ) return true;
  return false;
}


G4VParticleChange* G4FTFCaptureAtRest::AtRestDoIt( const G4Track& track, const G4Step& step ) {
  
  // Check applicability
  if ( ! IsApplicable( *(track.GetDynamicParticle()->GetDefinition()) ) ) {
    G4ExceptionDescription ed;
    ed << "Error: particle is: " << track.GetDynamicParticle()->GetDefinition()->GetParticleName()
       << "\t ; it must be an anti-proton ! " << G4endl;
    G4Exception( "G4FTFCaptureAtRest::AtRestDoIt()", "HAD_FTF_0000", 
                 FatalException, ed );
    return 0;
  }
  
  // Select the target nucleus
  G4Material * material = track.GetMaterial();
  G4Nucleus * targetNucleus = 0;
  do {
      targetNucleus = new G4Nucleus( material );
      if ( targetNucleus->GetA_asInt() < 1 ) {
        delete targetNucleus;
        targetNucleus = 0;
      }
  } while ( targetNucleus == 0 );    
  G4int targetNucleusZ = targetNucleus->GetZ_asInt();
  G4int targetNucleusA = targetNucleus->GetA_asInt();
  //G4cout << " targetNucleus Z=" << targetNucleusZ << " A=" << targetNucleusA << G4endl;

  // Prepare to call the interaction
  G4HadProjectile projectile( track );  
  G4HadFinalState* result = 0;
  G4int reentryCount = 0;
  do {
    try {
      // Call the interaction
      result = theModel->ApplyYourself( projectile, *targetNucleus );
      ++reentryCount;
    }
    catch( G4HadronicException aR ) {
      DumpState( track, "G4FTFCaptureAtRest::AtRestDoIt()" );
      G4ExceptionDescription ed;
      ed << "Call for " << theModel->GetModelName() << G4endl
         << "Target nucleus Z=" << targetNucleusZ 
         << "  A=" << targetNucleusA;
      G4Exception( "G4FTFCaptureAtRest::AtRestDoIt()", "HAD_FTF_0001", 
                   FatalException, ed );
    }
    if ( reentryCount > 100 ) {
      DumpState( track, "G4FTFCaptureAtRest::AtRestDoIt()" );
      G4ExceptionDescription ed;
      ed << "Reentering AtRestDoIt too often." << G4endl
         << "Call for " << theModel->GetModelName() << G4endl
         << "Target nucleus Z=" << targetNucleusZ 
         << "  A=" << targetNucleusA;
      G4Exception( "G4FTFCaptureAtRest::AtRestDoIt()", "HAD_FTF_0002", 
                   FatalException, ed );
    } 
  } while ( !result );

  // Transform from G4HadFinalState to G4ParticleChange .
  // (We follow the (private) method  G4HadronicProcess::FillResult
  // used by the method  G4HadronicProcess::PostStepDoIt .)

  aParticleChange.Clear();
  aParticleChange.Initialize( track );
  aParticleChange.ProposeLocalEnergyDeposit( result->GetLocalEnergyDeposit() );  
  // Check status of primary: it should have been killed
  if( result->GetStatusChange() == stopAndKill ) {
    aParticleChange.ProposeTrackStatus( fStopAndKill );
    aParticleChange.ProposeEnergy( 0.0 );
  } else {
    DumpState( track, "G4FTFCaptureAtRest::AtRestDoIt()" );
    G4ExceptionDescription ed;
    ed << "AtRestDoIt did not kill the absorbed particle." << G4endl
       << "Call for " << theModel->GetModelName() << G4endl
       << "Target nucleus Z=" << targetNucleusZ 
       << "  A=" << targetNucleusA;
    G4Exception( "G4FTFCaptureAtRest::AtRestDoIt()", "HAD_FTF_0003", 
                 FatalException, ed );
  }

  G4int nSec = result->GetNumberOfSecondaries();
  //G4cout << "#### ParticleDebug : number of secondaries = " << nSec
  //       << " ;  local energy deposit = " << result->GetLocalEnergyDeposit() 
  //       << " MeV" << G4endl;
  aParticleChange.SetNumberOfSecondaries( nSec );
  G4LorentzVector final4mom;
  if ( nSec > 0 ) {
    G4double time0 = track.GetGlobalTime();
    // Loop over the secondaries (which include the remnant nucleus)
    for ( G4int i=0; i < nSec; i++ ) {
      final4mom += result->GetSecondary(i)->GetParticle()->Get4Momentum();
      G4double time = result->GetSecondary(i)->GetTime();
      if ( time < time0) time = time0;
      G4Track* secTrack = new G4Track( result->GetSecondary(i)->GetParticle(),
                                       time, track.GetPosition() );
      G4double newWeight = track.GetWeight() * result->GetSecondary(i)->GetWeight();
      //G4cout << "#### ParticleDebug "
      //       << result->GetSecondary(i)->GetParticle()->GetDefinition()->GetParticleName()
      //       << " ;  weight=" << result->GetSecondary(i)->GetWeight()
      //       << " ;  4-momentum=" << result->GetSecondary(i)->GetParticle()->Get4Momentum()
      //       << G4endl;
      secTrack->SetWeight( newWeight );
      secTrack->SetTouchableHandle( track.GetTouchableHandle() );
      aParticleChange.AddSecondary( secTrack );
    }
  }

  // Check energy-momentum conservation
  G4LorentzVector projectile4mom = track.GetDynamicParticle()->Get4Momentum();
  G4double targetMass = G4NucleiProperties::GetNuclearMass( targetNucleusA, targetNucleusZ );
  G4LorentzVector target4mom( 0, 0, 0, targetMass ); // Neglect thermal motion
  G4LorentzVector initial4mom = projectile4mom + target4mom;
  G4LorentzVector diff = initial4mom - final4mom;
  const G4double threshold = 1.0*MeV;
  //G4cout << "===ANTI-PROTON CAPTURE AT REST=== : Ekin = " 
  //       << ( projectile4mom.e() - projectile4mom.mag() ) / MeV << G4endl;  // Debug
  if ( std::abs( diff.e() )  > threshold  ||
       std::abs( diff.px() ) > threshold  ||
       std::abs( diff.py() ) > threshold  ||
       std::abs( diff.pz() ) > threshold ) {
    //G4cout << "*** G4FTFCaptureAtRest::AtRestDoIt : 4-momentum non conservation " 
    //       << diff << G4endl
    //       << "    initial4mom = " << initial4mom << " ;  final4mom = " << final4mom  
    //       << G4endl;
  }

  result->Clear();

  //return &aParticleChange;                           // This is not enough
  return G4VRestProcess::AtRestDoIt( track, step );
}


void G4FTFCaptureAtRest::DumpState( const G4Track& aTrack , const G4String& method ) {
  G4cout << "Unrecoverable error in method " << method << G4endl
         << "TrackID= " << aTrack.GetTrackID() << "  ParentID= " << aTrack.GetParentID()
         << "  " << aTrack.GetParticleDefinition()->GetParticleName() << G4endl
         << "Ekin(GeV)= " << aTrack.GetKineticEnergy()/CLHEP::GeV 
         << ";  direction= " << aTrack.GetMomentumDirection() << G4endl
         << "Position(mm)= " << aTrack.GetPosition()/CLHEP::mm << ";";
  if ( aTrack.GetMaterial() ) G4cout << "  material " << aTrack.GetMaterial()->GetName();
  G4cout << G4endl;
  if ( aTrack.GetVolume() )
    G4cout << "PhysicalVolume  <" << aTrack.GetVolume()->GetName() << ">" << G4endl;
} 

