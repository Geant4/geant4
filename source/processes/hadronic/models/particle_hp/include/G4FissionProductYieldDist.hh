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
 * File:   G4FissionProductYieldDist.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on May 11, 2011, 12:04 PM
 */

#ifndef G4FISSIONPRODUCTYIELDDIST_HH
#define	G4FISSIONPRODUCTYIELDDIST_HH

#include "G4Ions.hh"
#include "G4Gamma.hh"
#include "G4IonTable.hh"
#include "G4ParticleHPNames.hh"
#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4DynamicParticleVector.hh"
#include "G4ReactionProduct.hh"

#include "G4ENDFTapeRead.hh"
#include "G4ENDFYieldDataContainer.hh"
#include "G4FFGEnumerations.hh"
#include "G4FPYNubarValues.hh"
#include "G4FPYSamplingOps.hh"
#include "G4FPYTreeStructures.hh"

/** G4FissionProductYieldDist is the base class for storing all the fission
 *  data and generating fission events. */
class G4FissionProductYieldDist
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
    G4FissionProductYieldDist( G4int WhichIsotope,
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
    G4FissionProductYieldDist( G4int WhichIsotope,
                               G4FFGEnumerations::MetaState WhichMetaState,
                               G4FFGEnumerations::FissionCause WhichCause,
                               G4FFGEnumerations::YieldType WhichYieldType,
                               G4int Verbosity,
                               std::istringstream& dataStream);
private:
    /** Initialize is a common function called by all constructors. */
    void Initialize( std::istringstream& dataStream );

public:
// Functions
    /** Generates a fission event using default sampling and returns the pointer
     *  to that fission event.
     *  - Usage: No arguments required
     *
     *  - Notes:
     *      - The fission products are loaded into FissionContainer in the
     *        following order:
     *          - First daughter product
     *          - Second daughter product
     *          - Alpha particles
     *          - Neutrons
     *          - Gamma rays
     *      - The last particle will have a NULL NextFragment pointer
     */
    G4DynamicParticleVector* G4GetFission( void );
    /** Selects a fission fragment at random from the probability tree and
     *  returns the \p G4Ions pointer.
     *  - Usage: No arguments required
     *
     *  - Notes:
     */
    G4Ions* G4GetFissionProduct( void );
    /** Set the alpha production behavior for fission event generation.
     *  - Usage:
     *      - if \p AlphaProduction is negative then alpha particles are sampled
     *        randomly.
     *
     *  - Notes:
     *      - The maximum number of alpha particles that may be created is
     *        physically limited by the nucleons present in the parent nucleus.
     *        Setting the AlphaProduction too high will have unpredictable
     *        results on the sampling of the fission products.
     */
    void G4SetAlphaProduction( G4double WhatAlphaProduction );
    /** Sets the energy of the incident particle
     *  - Usage:
     *      - \p WhatIncidentEnergy: Kinetic energy, if any, of the incident
     *      neutron in GeV
     *
     *  - Notes:
     */
    void G4SetEnergy( G4double WhatIncidentEnergy );
    /** Sets the probability of ternary fission
     *  - Usage:
     *      - \p WhatTernaryProbability: Probability of generating a ternary
     *        fission event.
     *
     *  - Notes:
     */
    void G4SetTernaryProbability( G4double TernaryProbability );
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
         *  G4FissionProductYieldDist references
         */
        const G4int Isotope_;
        /** MetaState information of the isotope that G4FissionProductYieldDist
         *  references
         *  \n Possible values are \p GROUND_STATE, \p META_1, or \p META_2
         */
        const G4FFGEnumerations::MetaState MetaState_;
        /** The cause of fission: \p SPONTANEOUS or \p N_INDUCED. */
        const G4FFGEnumerations::FissionCause Cause_;
        /** The type of yield to be used: \p INDEPENDET or \p CUMULATIVE */
        const G4FFGEnumerations::YieldType YieldType_;

    // Datafile variables
        /** Name of the fission yield product data file that
         *  G4FissionProductYieldDist references
         */
        G4ENDFTapeRead* ENDFData_;

    // Fission generation variables
        /** Contains the \p G4Ions pointer to an alpha particle */
        G4Ions* AlphaDefinition_;
        /** Controls whether alpha particles are emitted, and how many */
        G4double AlphaProduction_;
        /** Sets the ternary fission probability. Valid ranges are [0, 1] */
        G4double TernaryProbability_;
        /** Contains the \p g4ParticleDefinition pointer to a gamma particle */
        G4Gamma* GammaDefinition_;
        /** Kinetic energy, if any, of the incident particle in GeV. */
        G4double IncidentEnergy_;
        /** Sets the mean gamma energy, in MeV, produced by the fission of the
         *  isotope described by Isotope_
         */
        G4double MeanGammaEnergy_;
        /** Contains the G4ParticleDefinition pointer to a neutron, cast as a
         *  G4Ion for compatibility*/
        G4Ions* NeutronDefinition_;
        /** Nubar for the isotope and incident neutron energy that
         *  G4FissionProductYieldDist references.
         */
        G4double Nubar_;
        /** Width of the gaussian distribution that samples nubar for the
         *  isotope and incident neutron energy that G4FissionProductYieldDist
         *  references.
         */
        G4double NubarWidth_;
        /** Counter for the number of protons available to the fission event */
        G4int RemainingZ_;
        /** Counter for the number of nucleons available to the fission event */
        G4int RemainingA_;
        /** Container for the energy remaining to be assigned in the fission generation */
        G4double RemainingEnergy_;
        /** Verbosity level */
        G4int Verbosity_;

    // Pointers to the field of trees and relevant normalization data
        /** An array, or 'field', of the probability trees */
        ProbabilityTree* Trees_;
        /** Defines the smallest Z particle in the field of trees */
        G4Ions* SmallestZ_;
        /** Defines the smallest A particle in the field of trees */
        G4Ions* SmallestA_;
        /** Defines the largest Z particle in the field of trees. */
        G4Ions* LargestZ_;
        /** Defines the largest Z particle in the field of trees */
        G4Ions* LargestA_;
        /** Number of specific energy groups */
        G4int YieldEnergyGroups_;
        /** Energy values of each energy */
        G4double* YieldEnergies_;
        /** Variable for ensuring that the input data is normalized */
        G4double* MaintainNormalizedData_;
        /** A running total of all the probabilities */
        G4double* DataTotal_;
        /** The number of trees in the field */
        G4int TreeCount_;
        /** A run-time counter for the total number of branches stored */
        G4int BranchCount_;

    // Pointers to runtime classes that G4FissionProductYieldDist needs to
    // function properly
        /** Pointer to \p G4IonTable
         *  \n All \p G4Ions are created using
         *  \p G4IonTable
         */
        G4IonTable* IonTable_;
        /** Pointer to \p G4NeutronHPNames
         * \n Provides access to the list of element names included in Geant4
         */
        G4ParticleHPNames* ElementNames_;
        /** Pointer to the \p CLHEP library random engine */
        G4FPYSamplingOps* RandomEngine_;

//Functions
    /** Checks to make sure that alpha overpopulation will not occur, which
     *  could result in an unsolvable zero momentum in the LAB system.
     */
    void CheckAlphaSanity( void );
    /** Returns the \p G4Ions definitions pointer for the particle whose
     *  probability segment contains the (0, 1] random number \p RandomParticle
     */
    G4Ions* FindParticle( G4double RandomParticle );
    /** Returns the \p G4Ions definitions pointer for the particle whose
     *  probability segment contains the (0, 1] random number \p RandomParticle
     *  by extrapolating values using the current data set.
     *  This function exists so that that different models of extrapolation
     *  may be more easily implemented in the future.
     */
    G4Ions* FindParticleExtrapolation( G4double RandomParticle,
                                       G4bool LowerEnergyGroupExists );
    /** Returns the \p G4Ions definitions pointer for the particle whose
     *  probability segment contains the (0, 1] random number \p RandomParticle
     *  by interpolating values in the current data set.
     *  This function exists so that that different models of interpolation
     *  may be more easily implemented in the future.
     */
    G4Ions* FindParticleInterpolation( G4double RandomParticle,
                                       G4int LowerEnergyGroup );
    /** Returns the \p G4Ions definitions pointer for the particle whose
     *  probability segment contains the (0, 1] random number \p RandomParticle
     *  by searching through a branch. Both the extrapolation and interpolation
     *  schemes currently use this function to identify the particle.
     */
    G4Ions* FindParticleBranchSearch( ProbabilityBranch* Branch,
                                      G4double RandomParticle,
                                      G4int EnergyGroup1,
                                      G4int EnergyGroup2 );
    /** Generates a \p G4DynamicParticleVector with the fission alphas
     */
    virtual void GenerateAlphas( std::vector< G4ReactionProduct* >* Alphas );
    /** Generate a linked chain of neutrons and return the pointer to the last
     * neutron in the chain.
     */
    virtual void GenerateNeutrons( std::vector< G4ReactionProduct* >* Neutrons );
    /** Selects a fission product from the probability tree, limited by the
     *  number of nucleons available to the system
     */
    virtual G4Ions* GetFissionProduct( void ) = 0;
    /** Returns the \p G4Ions definition pointer to the isotope defined by
     *  \p Product and \p MetaState.
     *  Searches the \p ParticleTable for the particle defined by \p Product
     *  (ZZZAAA) and \p MetaState and returns the \p G4Ions
     *  pointer to that particle. If the particle does not exist then it is
     *  created in \p G4ParticleTable and the pointer to the new particle is
     *  returned.
     */
    G4Ions* GetParticleDefinition( G4int Product,
                                   G4FFGEnumerations::MetaState MetaState );
    /** Generates the directory location for the data file referenced by
     *  G4FissionProductYieldDist
     */
    G4String MakeDirectoryName( void );
    /** Generates the appropriate file name for the isotope requested */
    G4String MakeFileName( G4int Isotope,
                           G4FFGEnumerations::MetaState MetaState );
    /** Creates a \p G4DynamicParticle from an existing \p G4ReactionProduct */
    G4DynamicParticle* MakeG4DynamicParticle( G4ReactionProduct* );
    /** Generates the unique name for an isotope/isomer defined by \p Isotope\
     *  and \p MetaState in the following format: ZZZ_AAAmX_NAME
     */
    G4String MakeIsotopeName( G4int Isotope,
                              G4FFGEnumerations::MetaState MetaState );
    /** Dynamically allocates and initializes the 'field' of 'trees' with the
     *  'trunks'
     */
    virtual void MakeTrees( void );
    /** Reads in the probability data from the data file */
    virtual void ReadProbabilities( void );
    /** Renormalizes the data in a ProbabilityTree.
     *  Traverses the tree structure and renormalizes all the probability data
     *  into probability segments, ensuring that no segment overlaps the
     *  other.
     */
    void Renormalize( ProbabilityBranch* Branch );
    /** Sample the energy of the alpha particles. The energy used by the alpha
     *  particles is subtracted from the available energy
     */
    void SampleAlphaEnergies( std::vector< G4ReactionProduct* >* Alphas );
    /** Samples the energy of the gamma rays */
    void SampleGammaEnergies( std::vector< G4ReactionProduct* >* Gammas );
    /** Sample the energy of the neutrons using the Watt fission spectrum. The
     *  kinetic energy consumed is returned.
     */
    void SampleNeutronEnergies( std::vector< G4ReactionProduct* >* Neutrons );
    /** Sets the nubar values for the isotope referenced by
     *  G4FissionProductYieldDistdefined from the data sets defined in
     *  SpecialOps.hh
     */
    void SetNubar( void );
    /** Sorts information for a potential new particle into the correct tree */
    virtual void SortProbability( G4ENDFYieldDataContainer* YieldData );

// Destructor function(s)
public:
    /** Default deconstructor. It is a virtual function since
     *  G4FissionProductYieldDist is a parent class
     */
    virtual ~G4FissionProductYieldDist( void );
protected:
    /** Recursively burns each branch in a probability tree. */
    void BurnTree( ProbabilityBranch* Branch );
};

#endif	/* G4FISSIONPRODUCTYIELDDIST_HH */

