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
 * ============================================================================
 *
 *       Filename:  CexmcException.cc
 *
 *    Description:  cexmc exceptions
 *
 *        Version:  1.0
 *        Created:  23.12.2009 13:42:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include "CexmcException.hh"
#include "CexmcRunManager.hh"
#include "CexmcCommon.hh"


CexmcException::CexmcException( CexmcExceptionType  type_ ) : type( type_ )
{
}


CexmcException::~CexmcException() throw()
{
}


const char *  CexmcException::what( void ) const throw()
{
    switch ( type )
    {
    case CexmcSystemException :
        return CEXMC_LINE_START "Unspecified system exception.";
    case CexmcEventActionIsNotInitialized :
        return CEXMC_LINE_START "Event action was not initialized.";
    case CexmcCmdLineParseException :
        return CEXMC_LINE_START "Exception in command line arguments. "
               "Check your command line arguments.";
    case CexmcPreinitException :
        return CEXMC_LINE_START "Exception in PreInit phase. "
               "Probably you did not specify preinit macro or production model "
               "in it.";
    case CexmcFileCompressException :
        return CEXMC_LINE_START "File compress exception. "
               "Check if the file exists and output directory is writable.";
    case CexmcReadProjectIncomplete :
        return CEXMC_LINE_START "Read project is incomplete. "
               "Check if the read project has complete set of data files.";
    case CexmcProjectExists :
        return CEXMC_LINE_START "Project you try to write exists. "
               "Use command line option '-y' to force project override.";
    case CexmcCmdIsNotAllowed :
        return CEXMC_LINE_START "Command is not allowed here. "
                "Probably you tried to issue a command which would override "
                "read-only settings in the read project.";
    case CexmcBadAngularRange :
        return CEXMC_LINE_START "Bad angular range. "
                "Check specified angular ranges (probably they mismatch ranges "
                "in the read project).";
    case CexmcBadThreshold :
        return CEXMC_LINE_START "Bad threshold. "
                "Check specified thresholds (probably they mismatch thresholds "
                "in the read project).";
    case CexmcBadCalorimeterTriggerAlgorithm :
        return CEXMC_LINE_START "Bad calorimeter trigger algorithm. "
                "Check specified calorimeter trigger algorithm (probably it "
                "mismatches algorithm in the read project).";
    case CexmcBadOCVetoAlgorithm :
        return CEXMC_LINE_START "Bad outer crystals veto algorithm. "
                "Check specified outer crystals veto algorithm (probably it "
                "mismatches algorithm in the read project).";
    case CexmcBadOCVetoFraction :
        return CEXMC_LINE_START "Bad outer crystals veto fraction. "
                "Check specified outer crystals veto fraction (probably it "
                "has higher value than in the read project).";
    case CexmcCalorimeterRegionNotInitialized :
        return CEXMC_LINE_START "Calorimeter sensitive region was not "
                "initialized. Check gdml source file.";
    case CexmcCalorimeterGeometryDataNotInitialized :
        return CEXMC_LINE_START "Calorimeter geometry data was not "
                "initialized. Check gdml source file.";
    case CexmcMultipleDetectorRoles :
        return CEXMC_LINE_START "Multiple detector roles for one logical "
                "volume. Check gdml source file.";
    case CexmcKinematicsException :
        return CEXMC_LINE_START "Kinematics exception. "
                "Check your production model.";
    case CexmcPoorEventData :
        return CEXMC_LINE_START "Cannot create rich event data from poor "
                "original data. Set less value of event data verbose level.";
    case CexmcIncompatibleGeometry :
        return CEXMC_LINE_START "Incompatible geometry definition. "
                "Check your geometry data.";
    case CexmcIncompleteProductionModel :
        return CEXMC_LINE_START "Production model provided is incomplete. "
                "Check that your production model defines all particles.";
    case CexmcIncompatibleProductionModel :
        return CEXMC_LINE_START "Production model provided is incompatible "
                "with physical process that tries to bind it.";
    case CexmcBeamAndIncidentParticlesMismatch :
        return CEXMC_LINE_START "Beam and incident particles mismatch. "
                "Instantiated reconstructor will not give reliable results in "
                "this case. Make sure that beam particle and incident particle "
                "of provided production model are same.";
    case CexmcInvalidAngularRange :
        return CEXMC_LINE_START "An angular range is not valid. "
                "Check specified angular ranges.";
#ifdef CEXMC_USE_CUSTOM_FILTER
    case CexmcCFBadSource :
        return CEXMC_LINE_START "Custom filter source file does not exist or "
                "damaged. Check your command line arguments or source file.";
    case CexmcCFParseError :
        return CEXMC_LINE_START "Custom filter source file contains errors. "
                "Check your custom filter source file.";
    case CexmcCFUninitialized :
        return CEXMC_LINE_START "Custom filter was not properly initialized.";
    case CexmcCFUninitializedVector :
        return CEXMC_LINE_START "Custom filter: A vector variable was not "
                "initialized.";
    case CexmcCFUnexpectedContext :
        return CEXMC_LINE_START "Custom filter: a variable or function call is "
                "in wrong type context (scalar when a vector is expected or "
                "vice versa). Check your custom filter script.";
    case CexmcCFUnexpectedFunction :
        return CEXMC_LINE_START "Custom filter: a function is unknown or in "
                "wrong type context. Check your custom filter script.";
    case CexmcCFUnexpectedVariable :
        return CEXMC_LINE_START "Custom filter: a variable is unknown or in "
                "wrong type context. Check your custom filter script.";
    case CexmcCFUnexpectedVariableUsage :
        return CEXMC_LINE_START "Custom filter: a variable is used in wrong "
                "type context (scalar when a vector is expected or vice "
                "versa). Check your custom filter script.";
    case CexmcCFUnexpectedVectorIndex :
        return CEXMC_LINE_START "Custom filter: a vector variable with wrong "
                "index. Indices of vectors should start from 1. Check your "
                "custom filter script.";
#endif
    case CexmcWeirdException :
        return CEXMC_LINE_START "A weird exception occurred. "
                "The world must collapse now.";
    default :
        return CEXMC_LINE_START "Unknown exception.";
    }
}


void  ThrowExceptionIfProjectIsRead( CexmcExceptionType  type,
                                     G4bool  extraCond )
{
    CexmcRunManager *  runManager( static_cast< CexmcRunManager * >(
                                            G4RunManager::GetRunManager() ) );
    if ( runManager->ProjectIsRead() && extraCond )
        throw CexmcException( type );
}

