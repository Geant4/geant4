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
 * File:   G4FissionFragmentGenerator.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on May 11, 2011, 12:04 PM
 */

#ifndef G4FISSIONFRAGMENTGENERATOR_HH
#define	G4FISSIONFRAGMENTGENERATOR_HH

#include "G4Ions.hh"
#include "globals.hh"
#include "G4HadFinalState.hh"
#include "G4HadProjectile.hh"

#include "G4FFGEnumerations.hh"
#include "G4FissionProductYieldDist.hh"
#include "G4TableTemplate.hh"

/** G4FissionFragmentGenerator is the front end class to be used by the user for
 *  handling all fission event generation.
 *
 *  This class is intended to be instantiated for one type of fission event for
 *  as specific isotope/isomer, fission type, and incident neutron energy. For
 *  this reason no functions exist to change or modify these values once the
 *  class in constructed. A new class must be created by the user for each type
 *  of fission event, if such functionality is desired.
 */
class G4FissionFragmentGenerator{
public:
// Constructor definition
    /** Default constructor
     *  - Usage: No arguments required
     *
     *  - Notes:
     *      - There are methods that should be called to set operating
     *        parameters before generating any fission events with
     *        G4FissionFragmentGenerator. These are:
     *          - G4SetIsotope()
     *          - G4SetMetaState()
     *          - G4SetCause()
     *          - G4SetIncidentEnergy()
     *          - G4SetYieldType()
     *          - G4SetAlphaProduction()
     *          - G4SetAlphaProductionProbability()
     *          - G4SetSamplingScheme()
     *      - If any or all of these parameters are not set by the user, then
     *        default values will be used.
     *          - Isotope: \p 92238
     *          - Metastable state: \p GROUND_STATE
     *          - Cause: \p SPONTANEOUS
     *          - Incident energy: \p 0.025 eV
     *          - Yield type: \p INDEPENDENT
     *          - Alpha production: \p 0
     *          - Alpha production probability: \p 0
     *          - Sampling scheme: \p NORMAL
     */
    G4FissionFragmentGenerator( void );
    /** Overloaded constructor
     *  - Usage:
     *      - \p Verbosity: Verbosity level
     *
     *  - Notes:
            - Refer to the documentation for the default constructor for
     *        setting up the operating parameters.
     */
    G4FissionFragmentGenerator( G4int Verbosity );
protected:
    /** Initialize is a common function called by all constructors. */
    void Initialize( void );

public:
// Functions
    /** Generates a single fission event
     *  - Usage: No arguments required
     *
     *  - Notes:
     *      - Generates a single fission event by calling the overloaded function
     *        and passing an argument of '1'
     */
    G4DynamicParticleVector* G4GenerateFission( void );
    /** Generates a single fission event
     *  - Usage:
	 *		-\p Projectile: G4HadProjectile of the fission-inducing particle
     *
     *  - Notes:
     *      - Generates a single fission event by calling the overloaded function
     *        and passing an argument of '1'
     */
    G4DynamicParticleVector* G4GenerateFission( const G4HadProjectile& Projectile );
    /** Generates NumberOfFissions fission events
     *  - Usage:
     *      -\p NumberOfFissions: The number of fission events to generate
     *
     *  - Notes:
     *      - Generates \p NumberOfFissions fission events
     */
    const std::vector< G4DynamicParticleVector* > G4GenerateFission( G4long NumberOfFissions,
                                                                     const G4HadProjectile& Projectile );
    /** Returns a randomly sampled fission product */
    G4Ions* G4GenerateFissionProduct( void );
    /** Returns the production rate of alpha particles for fission events */
    G4double G4GetAlphaProduction( void );
    /** Returns the probability of ternary fission */
    G4double G4GetTernaryProbability( void );
    /** Returns the FissionCause of the fission event. */
    G4FFGEnumerations::FissionCause G4GetCause( void );
    /** Returns the energy of the fission inducing particle. */
    G4double G4GetIncidentEnergy( void );
    /** Returns the code of the fission isotope in ZZZAAA format. */
    G4int G4GetIsotope( void );
    /** Returns the MetaState of the fission isotope. */
    G4FFGEnumerations::MetaState G4GetMetaState( void );
    /** Returns the FissionSamplingScheme that is currently in use. */
    G4FFGEnumerations::FissionSamplingScheme G4GetSamplingScheme( void );
    /** Returns the yield type that is currently in use */
    G4FFGEnumerations::YieldType G4GetYieldType( void );
    /** Initializes a new \p G4FPY...Dist class based on the class descriptor
     *  variables of G4FissionFragmentGenerator.
     */
    bool InitializeFissionProductYieldClass( std::istringstream& dataFile );
    /** Converts the Z, A and M of an isotope into an integer representation **/
    static G4int G4MakeIsotopeCode(G4int Z, G4int A, G4int M);
    /** Sets the number of alpha particles produced in fission.
     *  - Usage:
     *      - if \p AlphaProduction is negative then alpha particles are sampled
     *        on a Gaussian with a mean of \p abs(AlphaProduction).
     *
     *  - Notes:
     *      - The maximum number of alpha particles that may be created is
     *        physically limited by the nucleons present in the parent nucleus.
     *        Setting the AlphaProduction too high will have unpredictable
     *        results on the sampling of the fission products.
     */
    void G4SetAlphaProduction( G4double WhatAlphaProduction );
    /** Sets the probability of ternary fission
     *  - Usage:
     *      - \p WhatAlphaProductionProbability: Probability of generating alpha
     *      particles for a fission event. 1 = 100% chance of alpha production
     *
     *  - Notes:
     */
    void G4SetTernaryProbability( G4double WhatTernaryProbability );
    /** Sets the cause of fission event.
     *  - Usage:
     *      - \p WhichCause: \p SPONTANEOUS, \p N_INDUCED, \p P_INDUCED, or
     *      \p G_INDUCED
     *
     *  - Notes:
     */
    void G4SetCause( G4FFGEnumerations::FissionCause WhichCause );
    /** Sets the incident energy, if any, of the particle that cause fission.
     *  - Usage:
     *      - \p WhatIncidentEnergy: Kinetic energy of the particle with units applied;
     *
     *  - Notes:
     */
    void G4SetIncidentEnergy( G4double WhatIncidentEnergy );
    /** Sets the fission isotope
     *  - Usage:
     *      - \p WhichIsotope: Code of the isotope in ZZZAAA format
     *
     *  - Notes:
     */
    void G4SetIsotope( G4int WhichIsotope );
    /** Sets the metastable state of the fission isotope.
     *  - Usage:
     *      - \p WhichMetaState: \p GROUND_STATE, \p META_1, or \p META_2
     *
     *  - Notes:
     */
    void G4SetMetaState( G4FFGEnumerations::MetaState WhichMetaState );
    /** Set the sampling scheme.
     *  - Usage:
     *      - NewScheme: The G4FissionSamplingScheme value for the sampling
     *        scheme to use.
     *
     *  - Notes:
     *      - \p NORMAL: Sets the parameters of this class to sample fission
     *           events without any biasing.
     *      - \p LIGHT_FRAGMENT: Sets the parameters of this class to bias the
     *           fragment generation by always selecting a light fragment
     *           (A \< 115) first.
     *      - \p WENDT: Sets the parameters of this class to sample fission
     *           events according to the Wendt sampling scheme. Please refer to
     *           the code documentation for G4FPYWendtSamplingDist for a more
     *           detailed explanation.
     */
    void G4SetSamplingScheme( G4FFGEnumerations::FissionSamplingScheme NewScheme );
    /** Sets the ENDF yield type to be used for the data
     *  - Usage:
     *      - \p WhichYieldType: \p INDEPENDENT or \p COMULATIVE
     *
     *  - Notes:
     */
    void G4SetYieldType( G4FFGEnumerations::YieldType WhichYieldType );
    /** Sets the verbosity levels
     *  - Usage:
     *      - \p WhichVerbosity: Combination of  levels
     *
     *  - Notes:
     *      - \p SILENT: All verbose output is repressed
     *      - \p UPDATES: Only high-level internal changes are reported
     *      - \p DAUGHTER_INFO: Displays information about daughter product sampling
     *      - \p NEUTRON_INFO: Displays information about neutron sampling
     *      - \p GAMMA_INFO: Displays information about gamma sampling
     *      - \p ALPHA_INFO: Displays information about alpha sampling
     *      - \p MOMENTUM_INFO: Displays information about momentum balancing
     *      - \p EXTRAPOLATION_INTERPOLATION_INFO: Displays information about any data extrapolation or interpolation that occurs
     *      - \p DEBUG: Reports program flow as it steps through functions
     *      - \p PRINT_ALL: Displays any and all output
     */
    void G4SetVerbosity( G4int WhatVerbosity );

protected:
// Variables
    // Class descriptor variables
        /** Number in ZZZAAA format of the isotope that
         *  G4FissionFragmentGenerator references
         */
        G4int Isotope_;
        /** MetaState information of the isotope that G4FissionFragmentGenerator
         *  references
         *  \n A value of 0 refers to the ground state
         */
        G4FFGEnumerations::MetaState MetaState_;
        /** The cause of fission: \p SPONTANEOUS or \p N_INDUCED. */
        G4FFGEnumerations::FissionCause Cause_;
        /** Kinetic energy, if any, of the incident particle in GeV. */
        G4double IncidentEnergy_;
        /** The type of yield to be used: \p INDEPENDET or \p CUMULATIVE */
        G4FFGEnumerations::YieldType YieldType_;
        /** Sets the ternary fission probability. Valid ranges are [0, 1] */
        G4double TernaryProbability_;
        /** Controls whether alpha particles are emitted, and how many */
        G4double AlphaProduction_;
        /** If Isotope_, MetaState_, Cause_, or IncidentEnergy_ are changed in
         *  the middle of a run then the class pointed at by YieldData_ will
         *  need to be reconstructed
         */
        G4bool IsReconstructionNeeded_;
        /** Verbosity level */
        G4int Verbosity_;

    // Defines the current sampling scheme and the respective class
        /** The sampling scheme that is used: \p NORMAL, \p LIGHT_FRAGMENT, or
         *  \p WENDT.
         */
        G4FFGEnumerations::FissionSamplingScheme SamplingScheme_;
        /** Pointer to G4FissionProductYieldDist class that holds all the
         *  probabilistic yield data
         */
        G4FissionProductYieldDist* YieldData_;

// Destructor function(s)
public:
    /** Default deconstructor */
    ~G4FissionFragmentGenerator();
};

#endif	/* G4FISSIONFRAGMENTGENERATOR_HH */

