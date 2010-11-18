/*
 * =============================================================================
 *
 *       Filename:  CexmcSensitiveDetectorsAttributes.hh
 *
 *    Description:  sensitive detectors attributes (types, roles and regions)
 *
 *        Version:  1.0
 *        Created:  06.10.2010 13:07:17
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_SENSITIVE_DETECTORS_ATTRIBUTES
#define CEXMC_SENSITIVE_DETECTORS_ATTRIBUTES

#include <G4String.hh>


enum  CexmcDetectorType
{
    CexmcEDDetector,
    CexmcTPDetector,
    CexmcNumberOfDetectorTypes
};


enum  CexmcDetectorRole
{
    CexmcMonitorDetectorRole,
    CexmcVetoCounterDetectorRole,
    CexmcCalorimeterDetectorRole,
    CexmcTargetDetectorRole,
    CexmcNumberOfDetectorRoles
};


const G4String  CexmcDetectorTypeName[ CexmcNumberOfDetectorTypes ] =
{
    "ED", "TP"
};


const G4String  CexmcDetectorRoleName[ CexmcNumberOfDetectorRoles ] =
{
    "Monitor", "VetoCounter", "Calorimeter", "Target"
};


const G4String  CexmcCalorimeterRegionName( "Calorimeter" );


#endif

