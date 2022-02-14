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
 * File:   G4FPYBiasedLightFragmentDist.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on June 2, 2011, 11:02 AM
 */

#ifndef G4FPYBIASEDLIGHTFRAGMENTDIST_HH
#define	G4FPYBIASEDLIGHTFRAGMENTDIST_HH

#include "G4Ions.hh"
#include "globals.hh"

#include "G4FFGEnumerations.hh"
#include "G4FissionProductYieldDist.hh"

/** G4FPYBiasedLightFragmentDist is an inherited class of G4FissionProductYield
 *  that only samples the 'light' fission fragments, defined by A \< 115
 *      - This inherited class of G4FissionProductYield samples only the lighter
 *      fission fragments, defined by A \< 115
 *      - This biasing was implemented because of an artifact that is introduced
 *      due to random sampling of fission fragments. Typically small fission
 *      fragments (neutrons, alphas, gammas) are sampled after the first
 *      fragment is sampled. If a heavy fragment (A \>= 115) is sampled first
 *      then the resulting lighter fission fragment, after all the other
 *      particles have been removed from the available mass, will most likely
 *      land far off the neutron drip line.
 *      - This implementation reduces the probability that such an improbable
 *      nucleus will be created by first sampling the lighter fission fragment
 *      and allowing the heavy fission fragment, which has a lot more
 *      flexibility for varying neutron populations, to make up the slack.
 */
class G4FPYBiasedLightFragmentDist
    : public G4FissionProductYieldDist
{
public:
// Constructor definition
    /** Default constructor
     *  - Usage:
     *      - \p WhichIsotope: Isotope number of the element in ZZZAAA form
     *      - \p WhichMetaState: \p GROUND_STATE, \p META_1, or \p META_2
     *      - \p WhichCause: \p SPONTANEOUS or \p N_INDUCED
     *      - \p WhichYieldType: \p INDEPENDENT or \p CUMULATIVE
     *
     *  - Notes:
     */
    G4FPYBiasedLightFragmentDist( G4int WhichIsotope,
                                  G4FFGEnumerations::MetaState WhichMetaState,
                                  G4FFGEnumerations::FissionCause WhichCause,
                                  G4FFGEnumerations::YieldType WhichYieldType,
                                  std::istringstream& dataStream);
    /** Overloaded constructor
     *  - Usage:
     *      - \p WhichIsotope: Isotope number of the element in ZZZAAA form
     *      - \p WhichMetaState: \p GROUND_STATE, \p META_1, or \p META_2
     *      - \p WhichCause: \p SPONTANEOUS or \p N_INDUCED
     *      - \p WhichYieldType: \p INDEPENDENT or \p CUMULATIVE
     *      - \p Verbosity: Verbosity level
     *
     *  - Notes:
     */
    G4FPYBiasedLightFragmentDist( G4int WhichIsotope,
                                  G4FFGEnumerations::MetaState WhichMetaState,
                                  G4FFGEnumerations::FissionCause WhichCause,
                                  G4FFGEnumerations::YieldType WhichYieldType,
                                  G4int Verbosity,
                                  std::istringstream& dataStream);
protected:
    /** Initialize is a common function called by all constructors. */
    void Initialize( void );

protected:
// Variables
    /** Defines the half-weight of the fission isotope */
    G4int HalfWeight_;
// Functions
    /** Selects a fission product from the probability tree, limited by the
     *  number of nucleons available to the system
     */
    virtual G4Ions* GetFissionProduct( void );

// Destructor function(s)
public:
    /** Default deconstructor. It is a virtual function since
     *  G4FPYBiasedLightFragmentDist inherits from G4FissionProductYieldDist
     */
    virtual ~G4FPYBiasedLightFragmentDist( void );
};

#endif	/* G4FPYBIASEDLIGHTFRAGMENTDIST_HH */

