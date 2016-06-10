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
 * File:   G4FPYSamplingOps.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on June 30, 2011, 11:10 AM
 */

/* * * * * * * * * * * * * * * *   References   * * * * * * * * * * * * * * * *
 *                                                                            *
 *  1.  "Sampling ENDL Watt Fission Spectra, D. E. Cullen, LLNL, April, 2004  *
 *                                                                            *
 * * * * * * * * * * * * * * * *   References   * * * * * * * * * * * * * * * */

#ifndef G4FPYSAMPLINGOPS_HH
#define	G4FPYSAMPLINGOPS_HH

#include "Randomize.hh"
#include "globals.hh"

#include "G4FFGEnumerations.hh"
#include "G4ShiftedGaussian.hh"
#include "G4WattFissionSpectrumValues.hh"

/** G4FPYSamplingOps performs all the uniform and Gaussian distribution sampling
 *  operations
 */
class G4FPYSamplingOps
{
public:
// Constructor definition
    /** Default constructor
     *  - Usage: No arguments required
     *  - Notes:
     */
    G4FPYSamplingOps( void );
    /** Overloaded constructor
     *  - Usage:
     *      - \p Verbosity: Verbosity level
     *  - Notes:
     */
    G4FPYSamplingOps( G4int Verbosity );
protected:
    /** Initialize is a common function called by all constructors. */
    void Initialize( void );

public:
// Functions
    /** Returns an integer value taken from a Gaussian distribution.
     *  This overloaded version assumes that the range is not restricted to
     *  positive values only.
     *  - Usage:
     *      - \p Mean: Mean about which the Gaussian distribution will be
     *        sampled
     *      - \p StdDev: Standard deviation of the Gaussian distribution. 68.3%
     *        of the values will lie within the first standard deviation, 95.4%
     *        within the second standard deviation, etc...
     *  - Notes:
     */
    G4int G4SampleIntegerGaussian( G4double Mean,
                                   G4double StdDev );
    /** Returns an integer value taken from a Gaussian distribution about
     *  \p Mean and with a standard deviation of \p StdDev.
     *  - Usage:
     *      - \p Mean: Mean about which the Gaussian distribution will be
     *        sampled
     *      - \p StdDev: Standard deviation of the Gaussian distribution. 68.3%
     *        of the values will lie within the first standard deviation, 95.4%
     *        within the second standard deviation, etc...
     *      - \p Range: \p POSITIVE or \p ALL
     *  - Notes:
     */
    G4int G4SampleIntegerGaussian( G4double Mean,
                                   G4double StdDev,
                                   G4FFGEnumerations::GaussianRange Range );
    /** Returns a double value taken from a Gaussian distribution about \p Mean
     *  and with a standard deviation of \p StdDev.
     *  - Usage:
     *      - \p Mean: Mean about which the Gaussian distribution will be
     *        sampled
     *      - \p StdDev: Standard deviation of the Gaussian distribution. 68.3%
     *        of the values will lie within the first standard deviation, 95.4%
     *        within the second standard deviation, etc...
     *  - Notes:
     */
    G4double G4SampleGaussian( G4double Mean,
                               G4double StdDev );
    /** Returns a double value taken from a Gaussian distribution about \p Mean
     *  and with a standard deviation of \p StdDev.
     *  - Usage:
     *      - \p Mean: Mean about which the Gaussian distribution will be
     *        sampled
     *      - \p StdDev: Standard deviation of the Gaussian distribution. 68.3%
     *        of the values will lie within the first standard deviation, 95.4%
     *        within the second standard deviation, etc...
     *      - \p Range: \p POSITIVE or \p ALL
     *  - Notes:
     */
    G4double G4SampleGaussian( G4double Mean,
                               G4double StdDev,
                               G4FFGEnumerations::GaussianRange Range );
    /** Returns a double value evenly distributed in the range (0, 1].
     *  - Usage: No arguments required
     *  - Notes:
     */
    G4double G4SampleUniform( void );
    /** Returns a double value evenly distributed in the range
     *  (\p Lower, \p Upper].
     *  - Usage:
     *      - \p Lower: Lower bounds of the distribution
     *      - \p Upper: Upper bounds of the distribution
     *
     *  - Notes:
     */
    G4double G4SampleUniform( G4double Lower,
                              G4double Upper );
    /** Samples the Watt fission spectrum for the selected isotope, using an
     *  algorithm adopted from Ref. 1
     *  - Usage:
     *      - \p WhatIsotope: The isotope that is to be sampled
     *      - \p WhatCause: The cause of the isotope to be sampled
     *      - \p WhatEnergy: The energy, in MeV of the incident particle
     *  - Notes:
     *      - All variables needed for this function are grouped together in
     *      WattConstants_.
     *      -
     */
    G4double G4SampleWatt( G4int WhatIsotope,
                           G4FFGEnumerations::FissionCause WhatCause,
                           G4double WhatEnergy );
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
        /** Mean for sampling a Gaussian distribution */
        G4double Mean_;
        /** Standard deviation for sampling a GaussianDistribution */
        G4double StdDev_;
        /** Structure chain that contains the all the previous values used
         *  for sampling a Gaussian distribution
         */
        G4ShiftedGaussian* ShiftedGaussianValues_;
        /** Verbosity level */
        G4int Verbosity_;
        /** Structure that contains the values for sampling the Watt fission
         *  spectrum
         */
        WattSpectrumConstants* WattConstants_;

    // Pointers to external classes
        /** Pointer to the CLHEP random number generator. */
        CLHEP::HepRandomEngine* RandomEngine_;

    // Internal variables for use with sampling a Gaussian distribution.
        /** Declares whether the second paired random number has been already
         *  returned.
         */
        G4bool NextGaussianIsStoredInMemory_;
        /** Contains the first of the two paired random numbers from the
         *  Gaussian distribution sampling.
         */
        G4double GaussianOne_;
        /** Contains the second of the two paired random numbers from the
         *  Gaussian distribution sampling.
         */
        G4double GaussianTwo_;
        /** Defines the tolerance that ShiftParameters() must match. */
        G4double Tolerance_;
// Functions
    /** Check to see if the user requested parameters have already been
     *  calculated. If they have, it recalls the stored parameters and sets
     *  them as the current values.
     */
    G4bool CheckAndSetParameters( void );
    /** Evaluates the constants that are required for the Watt fission spectrum
     *  sampling.
     */
    void EvaluateWattConstants( void );
    /** Samples a Gaussian distribution defined by the internal class variables
     *  NewMean_ and NewStdDev_.
     */
    G4double SampleGaussian( void );
    /** Sets the mean and standard deviation of the Gaussian distribution
     *  sampled by this class when \p POSITIVE values are requested.
     *  ShiftMean() performs two different operations based on the requested
     *  data type.
     *  - \p INTEGER: Iteratively searches for an adjusted mean that produces
     *    the same result as the mean requested by the implementor. In this
     *    instance the standard deviation is not adjusted.
     *  - \p DOUBLE: Adjusts the standard deviation of the Gaussian distribution
     *    so that the first seven standard deviations occur are all positive.
     *    The chance that a negative value will result using this method is
     *    2.56E<sup>-12</sup>
     */
    void ShiftParameters( G4FFGEnumerations::GaussianReturnType Type );

// Destructor function(s)
public:
    /** Default deconstructor. */
    ~G4FPYSamplingOps( void );
};

#endif	/* G4FPYSAMPLINGOPS_HH */

