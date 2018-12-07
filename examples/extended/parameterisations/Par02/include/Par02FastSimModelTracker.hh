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
/// \file Par02FastSimModelTracker.hh
/// \brief Definition of the Par02FastSimModelTracker class

#ifndef PAR02_TRACKER_FAST_SIM_MODEL_H
#define PAR02_TRACKER_FAST_SIM_MODEL_H

#include "G4VFastSimulationModel.hh"
#include "Par02DetectorParametrisation.hh"
#include "G4Step.hh"
#include "G4Navigator.hh"

/// Shortcut to the ordinary tracking for tracking detectors.
///
/// The fast simulation model describes what should be done instead of a
/// normal tracking. Instead of the ordinary tracking, a particle momentum
/// at the entrance of the tracking detector is smeared 
/// (by Par02Smearer::SmearMomentum()) and the particle is placed at the
/// tracking detector exit, at the place it would reach without the change
/// of its momentum. Based on G4 
/// examples/extended/parametrisations/Par01/include/Par01EMShowerModel.hh .
/// @author Anna Zaborowska

class Par02FastSimModelTracker : public G4VFastSimulationModel {
  public:
    
    /// A constructor.
    /// @param aModelName A name of the fast simulation model.
    /// @param aEnvelope A region where the model can take over the ordinary tracking.
    /// @param aParamType A parametrisation type.
    Par02FastSimModelTracker( G4String aModelName, G4Region* aEnvelope, 
                         Par02DetectorParametrisation::Parametrisation aParamType );
    
    /// A constructor.
    /// @param aModelName A name of the fast simulation model.
    /// @param aEnvelope A region where the model can take over the ordinary tracking.
    Par02FastSimModelTracker( G4String aModelName, G4Region* aEnvelope );
    
    /// A constructor.
    /// @param aModelName A name of the fast simulation model.
    Par02FastSimModelTracker( G4String aModelName );

    ~Par02FastSimModelTracker();
    
    /// Checks if this model should be applied to this particle type.
    /// @param aParticle A particle definition (type).
    virtual G4bool IsApplicable( const G4ParticleDefinition& aParticle );
    
    /// Checks if the model should be applied taking into account the kinematics
    /// of a track.
    /// @param aFastTrack A track.
    virtual G4bool ModelTrigger( const G4FastTrack& aFastTrack );
    
    /// Calculates the final position (at the outer boundary of the tracking detector)
    /// of a particle with the momentum at the entrance of the tracking detector.
    /// Smears the particle momentum and saves it, together with the tracking detector
    /// resolution and efficiency to the Par02PrimaryParticleInformation.
    /// @param aFastTrack A track.
    /// @param aFastStep A step.
    virtual void DoIt( const G4FastTrack& aFastTrack, G4FastStep& aFastStep );

  private:
    
    /// A pointer to Par02DetectorParametrisation used to get the efficiency and
    /// resolution of the tracking detector for a given particle and
    /// parametrisation type.
    Par02DetectorParametrisation* fCalculateParametrisation;
    
    /// A parametrisation type.
    Par02DetectorParametrisation::Parametrisation fParametrisation;
};

#endif

