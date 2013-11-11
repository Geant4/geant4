/*
 * File:   G4FPYNormalFragmentDist.cc
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on July 26, 2011, 12:26 PM
 */
 
#include "G4Ions.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "G4FFGDebuggingMacros.hh"
#include "G4FFGEnumerations.hh"
#include "G4FPYNormalFragmentDist.hh"
#include "G4FissionProductYieldDist.hh"

G4FPYNormalFragmentDist::
G4FPYNormalFragmentDist( G4int WhichIsotope,
                         G4FFGEnumerations::MetaState WhichMetaState,
                         G4FFGEnumerations::FissionCause WhichCause,
                         G4FFGEnumerations::YieldType WhichYieldType,
                         std::istringstream& dataFile)
:   G4FissionProductYieldDist( WhichIsotope,
                               WhichMetaState,
                               WhichCause,
                               WhichYieldType,
                               dataFile)
{
    // Initialize the class
    Initialize();
}

G4FPYNormalFragmentDist::
G4FPYNormalFragmentDist( G4int WhichIsotope,
                         G4FFGEnumerations::MetaState WhichMetaState,
                         G4FFGEnumerations::FissionCause WhichCause,
                         G4FFGEnumerations::YieldType WhichYieldType,
                         G4int Verbosity,
                         std::istringstream& dataFile)
:   G4FissionProductYieldDist( WhichIsotope,
                               WhichMetaState,
                               WhichCause,
                               WhichYieldType,
                               Verbosity,
                               dataFile)
{
    // Initialize the class
    Initialize();
}

void G4FPYNormalFragmentDist::
Initialize( void )
{
G4FFG_FUNCTIONENTER__

    // Nothing here yet

G4FFG_FUNCTIONLEAVE__
}

G4Ions* G4FPYNormalFragmentDist::
GetFissionProduct( void )
{
G4FFG_FUNCTIONENTER__

    G4Ions* Particle;

    // Generate a (0, 1] random number and return the respective particle.
    // The ENDF data tables lists 72172 as the largest fission fragment produced
    // for any fission event. The maximum alpha production is 10 and the
    // smallest fissile isotope is 90227. This means that if isotope 72172 were
    // selected as the first daughter product, then at 10 alpha particles only
    // 15 nucleons and -2 protons would remain for the second daughter product.
    // Although the actual probability of this occurring is very small, or 0 in
    // this case, a check should still be made to ensure that the second
    // daughter product can be physically realized. This would prevent a
    // situation such as this extreme example which results in a nucleus of 13
    // neutrons and 2 anti-protons.
    // This quick sanity check may become even more valid if the ENDF data
    // tables are expanded in the future and include larger fission products.
    do
    {
        Particle = FindParticle(RandomEngine_->G4SampleUniform());
    } while(Particle->GetAtomicMass() > RemainingA_ + 1
            || Particle->GetAtomicNumber() > RemainingZ_ + 1);

G4FFG_FUNCTIONLEAVE__
    return Particle;
}

G4FPYNormalFragmentDist::~G4FPYNormalFragmentDist( void )
{
G4FFG_FUNCTIONENTER__

    // Empty - all the data elements to be deconstructed are removed by
    // ~G4FissionProductYieldDist()
    
G4FFG_FUNCTIONLEAVE__
}
