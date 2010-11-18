/*
 * ============================================================================
 *
 *       Filename:  CexmcCommon.cc
 *
 *    Description:  global objects etc.
 *
 *        Version:  1.0
 *        Created:  03.12.2009 22:19:20
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <G4Allocator.hh>
#include "CexmcTrackPointInfo.hh"
#include "CexmcEnergyDepositStore.hh"
#include "CexmcTrackPointsStore.hh"

G4Allocator< CexmcTrackPointInfo >      trackPointInfoAllocator;
G4Allocator< CexmcEnergyDepositStore >  energyDepositStoreAllocator;
G4Allocator< CexmcTrackPointsStore >    trackPointsStoreAllocator;

