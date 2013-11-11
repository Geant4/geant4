/*
 * File:   G4ShiftedGaussian.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on July 20, 2011, 11:55 AM
 */

#ifndef G4SHIFTEDGAUSSIAN_HH
#define	G4SHIFTEDGAUSSIAN_HH

#include <utility>
#include <vector>

#include "globals.hh"

/** G4ShiftedGaussian is a class for storing the shifted values used for
 *  sampling a Gaussian distribution and returning only positive values; it is
 *  integrated into G4FPYSamplingOps
 */
class G4ShiftedGaussian
{
public:
// Constructor definition
    /** Default constructor
     *  - Usage: No arguments required
     *  - Notes:
     */
    G4ShiftedGaussian( void );
    /** Overloaded constructor
     *  - Usage:
     *      - \p Verbosity: Verbosity level
     *  - Notes:
     */
    G4ShiftedGaussian( G4int Verbosity );
protected:
    /** Initialize is a common function called by all constructors. */
    void Initialize( void );

public:
// Functions
    /** Returns the shifted mean that correlates to a \p RequestedMean and
     *  \p RequestedStdDev pair. 0 is returned if there is no associated value.
     */
    G4double G4FindShiftedMean( G4double RequestedMean,
                                G4double RequestedStdDev );
    /** Inserts a \p ShiftedMean indexed by the \p RequestedMean and
     * \p RequestedStdDev
     */
    void G4InsertShiftedMean( G4double ShiftedMean,
                              G4double RequestedMean,
                              G4double RequestedStdDev );
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
    /** Contains the adjusted mean of the POSITIVE only Gaussian distribution
     *  associated with a \p RequestedMean and \p RequestedStdDev pair.
     */
    std::vector<
        std::pair<
            std::pair<
                G4double,
                G4double
            >,
            G4double
        >
    > ShiftedMean_;
    /** Verbosity level */
    G4int Verbosity_;

// Destructor function(s)
public:
    /** Default deconstructor. */
    ~G4ShiftedGaussian( void );
};

#endif	/* G4SHIFTEDGAUSSIAN_HH */

