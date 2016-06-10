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
 * File:   G4ENDFTapeRead.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on September 6, 2011, 10:01 AM
 */

#ifndef G4ENDFTAPEREAD_HH
#define	G4ENDFTAPEREAD_HH

#include "globals.hh"

#include "G4ENDFYieldDataContainer.hh"
#include "G4FFGEnumerations.hh"
#include "G4TableTemplate.hh"

/** G4ENDFTapeRead is a class designed to read in data from unformatted ENDF data
 *  tapes for MT = 454 or MT = 459, which correspond to independent fission
 *  yields and cumulative fission yields, respectively. The data is stored
 *  internally and can be recalled one product at a time by calling
 *  G4GetNextYield().
 */
class G4ENDFTapeRead
{
public:
// Constructor definition
    /** Default constructor
     *  - Usage:
     *      - \p FileLocation: the absolute path to the file
     *      - \p FileName: the name of the data file
     *      - \p WhichYield: \p INDEPENDENT or \p CUMULATIVE
     *      - \p WhichCause: \p SPONTANEOUS or \p N_INDUCED
     *
     *  - Notes: The data will be read in immediately upon construction.
     */
    G4ENDFTapeRead( G4String FileLocation,
                    G4String FileName,
                    G4FFGEnumerations::YieldType WhichYield,
                    G4FFGEnumerations::FissionCause WhichCause );
    /** Overloaded constructor
     *  - Usage:
     *      - \p FileLocation: the absolute path to the file
     *      - \p FileName: the name of the data file
     *      - \p WhichYield: \p INDEPENDENT or \p CUMULATIVE
     *      - \p WhichCause: \p SPONTANEOUS or \p N_INDUCED
     *      - \p Verbosity: Verbosity level
     *
     *  - Notes: The data will be read in immediately upon construction.
     */
    G4ENDFTapeRead( G4String FileLocation,
                    G4String FileName,
                    G4FFGEnumerations::YieldType WhichYield,
                    G4FFGEnumerations::FissionCause WhichCause,
                    G4int Verbosity );
    /** Overloaded constructor
     *  - Usage:
     *      - \p DataFile: The absolute path to the data file
     *      - \p WhichYield: \p INDEPENDENT or \p CUMULATIVE
     *      - \p WhichCause: \p SPONTANEOUS or \p N_INDUCED
     *      - \p Verbosity: Verbosity level
     *
     *  - Notes: The data will be read in immediately upon construction.
     */
    G4ENDFTapeRead( std::istringstream& dataStream,
                    G4FFGEnumerations::YieldType WhichYield,
                    G4FFGEnumerations::FissionCause WhichCause,
                    G4int Verbosity );
protected:
    /** Initialize is a common function called by all constructors. */
    void Initialize( G4String dataFile );
    /** Initialize is a common function calles by all constructors */
    void Initialize( std::istringstream& dataStream );

public:
// Functions
    /** Returns and array containing the values of each of the energy groups
     *  - Usage: No arguments required
     *
     *  - Notes:
     */
    G4double* G4GetEnergyGroupValues( void );
    /** Returns the number of energy yield groups that were extracted from the
     *  ENDF tape file
     *  - Usage: No arguments required
     *
     *  - Notes:
     */
    G4int G4GetNumberOfEnergyGroups( void );
    /** Returns the number of fission products that were extracted from the
     *  ENDF tape file
     *  - Usage: No arguments required
     *
     *  - Notes:
     */
    G4int G4GetNumberOfFissionProducts( void );
    /** Returns the data for the requested fission product
     *  - Usage:
     *      - \p WhichYield: 0-based index of the fission product for which to
     *           get the yield data
     *
     *  - Notes:
     *      - This will return a pointer to the next G4FissionYieldContainer.
     *        NULL will be returned if no more fission containers exist.
     */
    G4ENDFYieldDataContainer* G4GetYield( G4int WhichYield );
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
    
private:
// Functions
    /** Read in the data from an ENDF data tape. */
    void ReadInData( std::istringstream& dataStream );

// Data members
    /** Stores the number corresponding to the fission cause that will be extracted */
    //const G4FFGEnumerations::FissionCause Cause_;
    /** Counter for the number of energy groups that were extracted */
    G4int EnergyGroups_;
    /** Array containing the values of the extracted energy groups */
    G4double* EnergyGroupValues_;
    /** Verbosity level */
    G4int Verbosity_;
    /** Storage for the extracted data */
    G4TableTemplate< G4ENDFYieldDataContainer >* YieldContainerTable_;
    /** Stores the number corresponding to the yield type that will be extracted */
    const G4FFGEnumerations::YieldType YieldType_;

// Destructor function(s)
public:
    /** Default Deconstructor */
    ~G4ENDFTapeRead( void );
};

#endif	/* G4ENDFTAPEREAD_HH */

