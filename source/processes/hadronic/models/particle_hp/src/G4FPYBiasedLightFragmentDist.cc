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
/*
 * File:   G4FPYBiasedLightFragmentDist.cc
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on June 2, 2011, 11:02 AM
 */
 
#include "G4Ions.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "G4FFGDebuggingMacros.hh"
#include "G4FFGEnumerations.hh"
#include "G4FPYBiasedLightFragmentDist.hh"
#include "G4FissionProductYieldDist.hh"

G4FPYBiasedLightFragmentDist::
G4FPYBiasedLightFragmentDist( G4int WhichIsotope,
                              G4FFGEnumerations::MetaState WhichMetaState,
                              G4FFGEnumerations::FissionCause WhichCause,
                              G4FFGEnumerations::YieldType WhichYieldType,
                              std::istringstream& dataStream )
:   G4FissionProductYieldDist( WhichIsotope,
                               WhichMetaState,
                               WhichCause,
                               WhichYieldType,
                               dataStream)
{
    // Initialize the class
    Initialize();
}

G4FPYBiasedLightFragmentDist::
G4FPYBiasedLightFragmentDist( G4int WhichIsotope,
                              G4FFGEnumerations::MetaState WhichMetaState,
                              G4FFGEnumerations::FissionCause WhichCause,
                              G4FFGEnumerations::YieldType WhichYieldType,
                              G4int Verbosity,
                              std::istringstream& dataStream)
:   G4FissionProductYieldDist( WhichIsotope,
                               WhichMetaState,
                               WhichCause,
                               WhichYieldType,
                               Verbosity,
                               dataStream )
{
    // Initialize the class
    Initialize();
}

void G4FPYBiasedLightFragmentDist::
Initialize( void )
{
G4FFG_FUNCTIONENTER__

    // Initialize the half-weight for fission product sampling
    HalfWeight_ = (G4int)floor((Isotope_ % 1000) / 2.0);

G4FFG_FUNCTIONLEAVE__
}

G4Ions* G4FPYBiasedLightFragmentDist::
GetFissionProduct( void )
{
G4FFG_FUNCTIONENTER__

    G4Ions* Particle;
    G4bool IsHeavy;
    G4bool IsNotFeasable;
    G4int Counter;

    Counter = 0;
    do
    {
        // If we have sampled too many times then we may need to do something else
        if(Counter == 1000)
        {
            Particle = NULL;
            break;
        }
        
        // Generate a (0, 1] random number, then find the corresponding particle
        Particle =  FindParticle(RandomEngine_->G4SampleUniform());

        // If the daughter is heavier than half of the original nucleus weight
        // or not physically feasable then set the flags to sample again
        IsHeavy = (Particle->GetAtomicMass() > HalfWeight_);
        IsNotFeasable = (Particle->GetAtomicMass() > RemainingA_ + 1
            || Particle->GetAtomicNumber() > RemainingZ_ + 1);
        
        Counter++;
    } while (IsHeavy || IsNotFeasable); // Loop checking, 11.05.2015, T. Koi
    
    // Something went wrong, so figure out how to fix it
    if(Particle == NULL)
    {
        /// \todo   Figure out what to do if a suitable particle is not found in 1000 iterations
    }

    // Return the G4ParticleDefintion pointer to the random light fragment
G4FFG_FUNCTIONLEAVE__
    return Particle;
}

G4FPYBiasedLightFragmentDist::~G4FPYBiasedLightFragmentDist( void )
{
G4FFG_FUNCTIONENTER__

    // Empty - all the data elements to be deconstructed are removed by
    // ~G4FissionProductYieldDist()
G4FFG_FUNCTIONLEAVE__
}
