/*
 * ============================================================================
 *
 *       Filename:  CexmcPhysicsManager.cc
 *
 *    Description:  interface for external access to physics aspects
 *
 *        Version:  1.0
 *        Created:  03.11.2009 14:22:21
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include "CexmcPhysicsManager.hh"
#include "CexmcPhysicsManagerMessenger.hh"
#include "CexmcCommon.hh"

CexmcPhysicsManager::CexmcPhysicsManager() : basicMaxIL( CexmcDblMax ),
    maxILCorrection( 0 ), proposedMaxIL( CexmcDblMax ),
    numberOfTriggeredStudiedInteractions( 0 ),
    onlyBeamParticleCanTriggerStudiedProcess( false ), messenger( NULL )
{
    messenger = new CexmcPhysicsManagerMessenger( this );
}


CexmcPhysicsManager::~CexmcPhysicsManager()
{
    delete messenger;
}

