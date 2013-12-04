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
 * File:   G4FFGDebuggingMacros.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on August 17, 2012, 12:54
 */

#ifndef G4FFGDEBUGGINGMACROS_HH
#define G4FFGDEBUGGINGMACROS_HH

#include "globals.hh"

#include "G4FFGEnumerations.hh"
#include "G4FFGVerboseMacros.hh"

// Define the function as available by the compiler
#if defined(__GNUC__)
    #define G4FFG_FUNCTION_SIGNATURE__ G4String(__func__) + "()"
#elif defined(_MSC_VER)
    // I'm not sure if this is the correct syntax for MS VC
    #define G4FFG_FUNCTION_SIGNATURE__ G4String(__FUNCTION__) + "()"
#else
    #define G4FFG_FUNCTION_SIGNATURE__ "a function"
#endif

// Only define the variables and macros if G4DEBUG_VERBOSE is set
#if defined(G4DEBUG_VERBOSE)
    /** G4FFGDEBUG_RECURSIVE_REPRESSION is used to aid in the repression of
     *  debugging messages from recursive function calls.
     */
    extern G4long G4FFGDEBUG_RECURSIVE_REPRESSION;
    /** G4FFGDEBUG_RECURSIVE_REPRESSION_COUNTER tracks the number of recursive
     *  function debugging messages were repressed.
     */
    extern G4long G4FFGDEBUG_RECURSIVE_REPRESSION_COUNTER;
    /** G4FFGDEBUG_DATA_STRUCTURE_REPRESSION is used to aid in the repression of
     *  debugging messages from functions that sort/access data elements.
     */
    extern G4long G4FFGDEBUG_DATA_STRUCTURE_REPRESSION;
    /** G4FFGDEBUG_DATA_STRUCTURE_REPRESSION_COUNTER tracks the number of recursive
     *  function debugging messages were repressed.
     */
    extern G4long G4FFGDEBUG_DATA_STRUCTURE_REPRESSION_COUNTER;

// Entering functions
    /** G4FFG_FUNCTIONENTER__ is blank if G4DEBUG_VERBOSE is not defined at compile time.
     *  Otherwise, it is used by the fission fragment code to notify when a function is
     *  first entered.
     */
    #define G4FFG_FUNCTIONENTER__ \
    if((Verbosity_ & G4FFGEnumerations::DEBUG) && !(Verbosity_ & G4FFGEnumerations::REPRESS_FUNCTION_ENTER_LEAVE_MESSAGES)) \
    { \
        G4FFG_SPACING__ \
        G4cout << "Entering ";\
        G4FFG_LOCATION__ \
        G4cout << G4endl; \
    } \
    G4FFG_DEPTH++;
    
    /** G4FFG_SAMPLING_FUNCTIONENTER__ wraps around G4FFG_FUNCTIONENTER__ and
     *  can be used in conjunctions with
     *  G4FFGEnumeration::REPRESS_RANDOM_SAMPLING_MESSAGES to repress debugging output
     *  for psuedorandom number generation functions
     */
    #define G4FFG_SAMPLING_FUNCTIONENTER__ \
    if(!(Verbosity_ & G4FFGEnumerations::REPRESS_RANDOM_SAMPLING_MESSAGES)) \
    { \
        G4FFG_FUNCTIONENTER__ \
    }
    
    /** G4FFG_RECURSIVE_FUNCTIONENTER__ wraps around G4FFG_FUNCTIONENTER__ and
     *  can be used in conjunctions with
     *  G4FFGEnumeration::REPRESS_RECURSIVE_DEBUG_MESSAGES to repress debugging output
     *  recursive function calls.
     */
    #define G4FFG_RECURSIVE_FUNCTIONENTER__ \
    if(Verbosity_ & G4FFGEnumerations::REPRESS_RECURSIVE_DEBUG_MESSAGES)\
    { \
        if(G4FFGDEBUG_RECURSIVE_REPRESSION == 0) \
        { \
            G4FFG_FUNCTIONENTER__ \
        } else \
        { \
            G4FFGDEBUG_RECURSIVE_REPRESSION_COUNTER++; \
        } \
        G4FFGDEBUG_RECURSIVE_REPRESSION++; \
    } else \
    { \
        G4FFG_FUNCTIONENTER__ \
    }
    
    /** G4FFG_DATA_FUNCTIONENTER__ wraps around G4FFG_FUNCTIONENTER__ and
     *  can be used in conjunctions with
     *  G4FFGEnumeration::REPRESS_DATA_STRUCTURE_DEBUG_MESSAGES to repress debugging output
     *  recursive function calls.
     */
    #define G4FFG_DATA_FUNCTIONENTER__ \
    if(Verbosity_ & G4FFGEnumerations::REPRESS_DATA_STRUCTURE_DEBUG_MESSAGES)\
    { \
        if(G4FFGDEBUG_RECURSIVE_REPRESSION == 0) \
        { \
            G4FFG_FUNCTIONENTER__ \
        } else \
        { \
            G4FFGDEBUG_DATA_STRUCTURE_REPRESSION_COUNTER++; \
        } \
        G4FFGDEBUG_RECURSIVE_REPRESSION++; \
    } else \
    { \
        G4FFG_FUNCTIONENTER__ \
    }

// Leaving functions
    /** G4FFG_FUNCTIONLEAVE__ is blank if G4DEBUG_VERBOSE is not defined at compile time.
     *  Otherwise, it is used by the fission fragment code to notify when a function is
     *  exited. It will also be found before \p return statements, since those exit a
     *  funtion as well.
     */
    #define G4FFG_FUNCTIONLEAVE__ \
    G4FFG_DEPTH--; \
    if((Verbosity_ & G4FFGEnumerations::DEBUG) && !(Verbosity_ & G4FFGEnumerations::REPRESS_FUNCTION_ENTER_LEAVE_MESSAGES)) \
    { \
        G4FFG_SPACING__ \
        G4cout << "Leaving ";\
        G4FFG_LOCATION__ \
        G4cout << G4endl; \
    }
    
    /** G4FFG_SAMPLING_FUNCTIONLEAVE__ wraps around G4FFG_FUNCTIONLEAVE__ and
     *  can be used in conjunctions with
     *  G4FFGEnumeration::REPRESS_RANDOM_SAMPLING_MESSAGES to repress debugging output
     *  for psuedorandom number generation functions
     */
    #define G4FFG_SAMPLING_FUNCTIONLEAVE__ \
    if(!(Verbosity_ & G4FFGEnumerations::REPRESS_RANDOM_SAMPLING_MESSAGES)) \
    { \
        G4FFG_FUNCTIONLEAVE__ \
    }
    
    /** G4FFG_RECURSIVE_FUNCTIONLEAVE__ wraps around G4FFG_FUNCTIONLEAVE__ and
     *  can be used in conjunctions with
     *  G4FFGEnumeration::REPRESS_RECURSIVE_DEBUG_MESSAGES to repress debugging output
     *  recursive function calls.
     */
    #define G4FFG_RECURSIVE_FUNCTIONLEAVE__ \
    if(Verbosity_ & G4FFGEnumerations::REPRESS_RECURSIVE_DEBUG_MESSAGES)\
    { \
        G4FFGDEBUG_RECURSIVE_REPRESSION--; \
        if(G4FFGDEBUG_RECURSIVE_REPRESSION == 0) \
        { \
            if(G4FFGDEBUG_RECURSIVE_REPRESSION_COUNTER > 0) \
            { \
                G4FFG_SPACING__ \
                G4cout << "==== " << G4FFGDEBUG_RECURSIVE_REPRESSION_COUNTER  * 2 << " recursive function messages suppressed ====" << G4endl;  \
                G4FFGDEBUG_RECURSIVE_REPRESSION_COUNTER = 0; \
            } \
            G4FFG_FUNCTIONLEAVE__ \
        } \
    } else \
    { \
        G4FFG_FUNCTIONLEAVE__ \
    }
    
    /** G4FFG_DATA_FUNCTIONLEAVE__ wraps around G4FFG_FUNCTIONLEAVE__ and
     *  can be used in conjunctions with
     *  G4FFGEnumeration::REPRESS_DATA_STRUCTURE_DEBUG_MESSAGES to repress debugging output
     *  recursive function calls.
     */
    #define G4FFG_DATA_FUNCTIONLEAVE__ \
    if(Verbosity_ & G4FFGEnumerations::REPRESS_DATA_STRUCTURE_DEBUG_MESSAGES)\
    { \
        G4FFGDEBUG_RECURSIVE_REPRESSION--; \
        if(G4FFGDEBUG_RECURSIVE_REPRESSION == 0) \
        { \
            if(G4FFGDEBUG_DATA_STRUCTURE_REPRESSION_COUNTER > 0) \
            { \
                G4FFG_SPACING__ \
                G4cout << "==== " << G4FFGDEBUG_DATA_STRUCTURE_REPRESSION_COUNTER * 2 << " data structure function messages suppressed ====" << G4endl;  \
                G4FFGDEBUG_DATA_STRUCTURE_REPRESSION_COUNTER = 0; \
            } \
            G4FFG_FUNCTIONLEAVE__ \
        } \
    } else \
    { \
        G4FFG_FUNCTIONLEAVE__ \
    }
#else /* G4DEBUG_VERBOSE */
// If G4DEBUG_VERBOSE is not defined then we will need to define these macros but leave them empty
// Except for G4FFG_FUNCTIONENTER__ and G4FFG_FUNCTIONLEAVE__, which will be used to track G4FFG_DEPTH
	#define G4FFG_FUNCTIONENTER__ G4FFG_DEPTH++;
	#define G4FFG_SAMPLING_FUNCTIONENTER__
	#define G4FFG_RECURSIVE_FUNCTIONENTER__
	#define G4FFG_DATA_FUNCTIONENTER__
	#define G4FFG_FUNCTIONLEAVE__ G4FFG_DEPTH--;
	#define G4FFG_SAMPLING_FUNCTIONLEAVE__
	#define G4FFG_RECURSIVE_FUNCTIONLEAVE__
	#define G4FFG_DATA_FUNCTIONLEAVE__
#endif /* G4DEBUG_VERBOSE */

#endif /* G4FFGDEBUGGINGMACROS_HH */

