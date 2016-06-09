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
 * =============================================================================
 *
 *       Filename:  CexmcException.hh
 *
 *    Description:  cexmc exceptions
 *
 *        Version:  1.0
 *        Created:  04.11.2009 00:00:58
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_EXCEPTION_HH
#define CEXMC_EXCEPTION_HH

#include <stdexcept>
#include <G4Types.hh>


enum  CexmcExceptionType
{
    CexmcUnknownException,
    CexmcSystemException,
    CexmcEventActionIsNotInitialized,
    CexmcCmdLineParseException,
    CexmcPreinitException,
    CexmcFileCompressException,
    CexmcReadProjectIncomplete,
    CexmcProjectExists,
    CexmcCmdIsNotAllowed,
    CexmcBadAngularRange,
    CexmcBadThreshold,
    CexmcBadCalorimeterTriggerAlgorithm,
    CexmcBadOCVetoAlgorithm,
    CexmcBadOCVetoFraction,
    CexmcCalorimeterRegionNotInitialized,
    CexmcCalorimeterGeometryDataNotInitialized,
    CexmcMultipleDetectorRoles,
    CexmcKinematicsException,
    CexmcPoorEventData,
    CexmcIncompatibleGeometry,
    CexmcIncompleteProductionModel,
    CexmcIncompatibleProductionModel,
    CexmcBeamAndIncidentParticlesMismatch,
    CexmcInvalidAngularRange,
#ifdef CEXMC_USE_CUSTOM_FILTER
    CexmcCFBadSource,
    CexmcCFParseError,
    CexmcCFUninitialized,
    CexmcCFUninitializedVector,
    CexmcCFUnexpectedContext,
    CexmcCFUnexpectedFunction,
    CexmcCFUnexpectedVariable,
    CexmcCFUnexpectedVariableUsage,
    CexmcCFUnexpectedVectorIndex,
#endif
    CexmcWeirdException
};


class  CexmcException : public std::exception
{
    public:
        explicit CexmcException( CexmcExceptionType  type );

        ~CexmcException() throw();

    public:
        const char *  what( void ) const throw();

    private:
        CexmcExceptionType  type;
};


void  ThrowExceptionIfProjectIsRead( CexmcExceptionType  type,
                                     G4bool  extraCond = true );


#endif

