// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PAIclusterModel.cc,v 1.3 2001-02-01 11:21:22 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4Timer.hh"

#include "G4PAIclusterModel.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "globals.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VTouchable.hh"


////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4PAIclusterModel::G4PAIclusterModel(G4Envelope *anEnvelope) :
  G4VClusterModel("G4PAIclusterModel",anEnvelope)
{
  fPAIonisation = new G4PAIonisation(anEnvelope->GetMaterial()->GetName()) ;
}

///////////////////////////////////////////////////////////////////////////

G4PAIclusterModel::~G4PAIclusterModel()
{
  delete fPAIonisation ;
}

///////////////////////////////////////////////////////////////////////////////
//
// Returns condition for application of the model depending on particle type


G4bool G4PAIclusterModel::IsApplicable(const G4ParticleDefinition& particle)
{
  return  ( particle.GetPDGCharge() != 0.0 && particle.GetPDGMass() != 0 ) ; 
}

/////////////////////////////////////////////////////////////////////
//
// UserTrigger() method: method which has to decide if
// the parameterisation has to be applied.
// Here ModelTrigger() asks the user (ie you) a 0/1 answer.
//
// Note that quantities like the local/global position/direction etc..
// are available at this level via the fastTrack parameter (allowing 
// to check distance from boundaries, see below  to allow the decision)
//

G4bool G4PAIclusterModel::ModelTrigger(const G4FastTrack& fastTrack) 
{
  G4double kinEnergy, mass, gamma ;
  kinEnergy  = fastTrack.GetPrimaryTrack()->GetKineticEnergy() ;
  mass       = fastTrack.GetPrimaryTrack()->GetDefinition()->GetPDGMass() ;
  gamma      = 1.0 + kinEnergy/mass ;
  if (gamma >= 1.2) return true  ;
  else              return false ;
}

//////////////////////////////////////////////////////////////////////////////
//
// 

void G4PAIclusterModel::DoIt( const G4FastTrack& fastTrack , 
		                    G4FastStep&  fastStep         )
{
  G4double charge, charge2, kinEnergy, mass, massRatio, scaledTkin ;
  G4double distance, energyTransfer, energyLoss, lambda, step, stepSum = 0.0 ;
  G4ThreeVector clusterPosition ;

  fClusterPositionVector.clear() ;
  fClusterEnergyVector.clear() ;

  charge     = fastTrack.GetPrimaryTrack()->GetDefinition()->GetPDGCharge() ;
  charge2    = charge*charge ;
  kinEnergy  = fastTrack.GetPrimaryTrack()->GetKineticEnergy() ;
  mass       = fastTrack.GetPrimaryTrack()->GetDefinition()->GetPDGMass() ;
  massRatio  = proton_mass_c2/mass ;
  scaledTkin = kinEnergy*massRatio ;

  G4ParticleMomentum direction(fastTrack.GetPrimaryTrackLocalDirection());

  distance = fastTrack.GetEnvelopeSolid()->
             DistanceToOut(fastTrack.GetPrimaryTrackLocalPosition(),direction) ;

  G4ThreeVector position = fastTrack.GetPrimaryTrackLocalPosition() + 
                           distance*direction ;

  // Set final position:

  fastStep.SetPrimaryTrackFinalPosition(position);

  // Cluster counting loop

  lambda = fPAIonisation->GetFreePath( scaledTkin, charge2 ) ;
  step   = RandExponential::shoot(lambda) ;
  //  if (step < 0.0) step = 0.0 ;
  stepSum += step ;
  //  distance -= stepSum ;
  //  if(distance < 0.0) // no change, return 
  if(stepSum > distance ) // no change, return 
  {
    return ;  
  }
  else
  {
    G4ThreeVector globalStartPosition  = fastTrack.GetPrimaryTrack()->
                                         GetPosition() ;
    G4ParticleMomentum globalDirection = fastTrack.GetPrimaryTrack()->
                                         GetMomentumDirection() ; 

    //  while(distance >= 0.0)  
    while(stepSum <= distance )  // global (or local ?) cluster coordinates
    {
      //  clusterPosition = fastTrack.GetPrimaryTrackLocalPosition() + 
      //                  stepSum*direction ;  
  
      clusterPosition = globalStartPosition  + stepSum*globalDirection ;    
      energyTransfer  = fPAIonisation->GetRandomEnergyTransfer(scaledTkin) ;

      fClusterPositionVector.insert(clusterPosition) ;      
      fClusterEnergyVector.insert(energyTransfer) ;

      step = RandExponential::shoot(lambda) ;
      // if (step < 0.0) step = 0.0 ;
      stepSum    += step ;
      //  distance   -= step ;
      energyLoss += energyTransfer ;     
    } 
    kinEnergy -= energyLoss ;
    fastStep.SetPrimaryTrackFinalKineticEnergy(kinEnergy) ;  
    // fastStep.SetTotalEnergyDeposited(energyLoss);

    BuildDetectorResponse() ; 
  }
  return ;
}


//
//
///////////////////////////////////////////////////////////////////////












