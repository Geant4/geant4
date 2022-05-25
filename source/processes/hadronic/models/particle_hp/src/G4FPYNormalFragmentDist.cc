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

    G4Ions* Particle=nullptr;

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

    G4int icounter=0;
    G4int icounter_max=1024;
    do
    {
      icounter++;
      if ( icounter > icounter_max ) {
	 G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
         break;
      }
        Particle = FindParticle(RandomEngine_->G4SampleUniform());
    } while(Particle->GetAtomicMass() > RemainingA_ + 1
            || Particle->GetAtomicNumber() > RemainingZ_ + 1);
    // Loop checking, 11.05.2015, T. Koi

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
