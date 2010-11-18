/*
 * =============================================================================
 *
 *       Filename:  CexmcEnergyDepositStore.hh
 *
 *    Description:  store energy deposit data and const references to
 *                  energy deposit collections in calorimeters
 *
 *        Version:  1.0
 *        Created:  25.11.2009 15:32:51
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_ENERGY_DEPOSIT_STORE_HH
#define CEXMC_ENERGY_DEPOSIT_STORE_HH

#include <G4Allocator.hh>
#include "CexmcCommon.hh"


struct  CexmcEnergyDepositStore
{
    CexmcEnergyDepositStore( G4double  monitorED,
                             G4double  vetoCounterEDLeft,
                             G4double  vetoCounterEDRight,
                             G4double  calorimeterEDLeft,
                             G4double  calorimeterEDRight,
                             G4int     calorimeterEDLeftMaxX,
                             G4int     calorimeterEDLeftMaxY,
                             G4int     calorimeterEDRightMaxX,
                             G4int     calorimeterEDRightMaxY,
                             const CexmcEnergyDepositCalorimeterCollection &
                                       calorimeterEDLeftCollection,
                             const CexmcEnergyDepositCalorimeterCollection &
                                       calorimeterEDRightCollection ) :
        monitorED( monitorED ), vetoCounterEDLeft( vetoCounterEDLeft ),
        vetoCounterEDRight( vetoCounterEDRight ),
        calorimeterEDLeft( calorimeterEDLeft ),
        calorimeterEDRight( calorimeterEDRight ),
        calorimeterEDLeftMaxX( calorimeterEDLeftMaxX ),
        calorimeterEDLeftMaxY( calorimeterEDLeftMaxY ),
        calorimeterEDRightMaxX( calorimeterEDRightMaxX ),
        calorimeterEDRightMaxY( calorimeterEDRightMaxY ),
        calorimeterEDLeftCollection( calorimeterEDLeftCollection ),
        calorimeterEDRightCollection( calorimeterEDRightCollection )
    {}

    void *  operator new( size_t  size );

    void    operator delete( void *  obj );

    G4double  monitorED;

    G4double  vetoCounterEDLeft;

    G4double  vetoCounterEDRight;

    G4double  calorimeterEDLeft;

    G4double  calorimeterEDRight;

    G4int     calorimeterEDLeftMaxX;

    G4int     calorimeterEDLeftMaxY;

    G4int     calorimeterEDRightMaxX;

    G4int     calorimeterEDRightMaxY;

    const CexmcEnergyDepositCalorimeterCollection &
              calorimeterEDLeftCollection;

    const CexmcEnergyDepositCalorimeterCollection &
              calorimeterEDRightCollection;
};


extern G4Allocator< CexmcEnergyDepositStore >  energyDepositStoreAllocator;


inline void *  CexmcEnergyDepositStore::operator new( size_t )
{
  return energyDepositStoreAllocator.MallocSingle();
}


inline void  CexmcEnergyDepositStore::operator delete( void *  obj )
{
    energyDepositStoreAllocator.FreeSingle(
                        reinterpret_cast< CexmcEnergyDepositStore * >( obj ) );
}


#endif

