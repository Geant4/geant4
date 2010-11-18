/*
 * =============================================================================
 *
 *       Filename:  CexmcEnergyDepositInLeftRightSet.hh
 *
 *    Description:  energy deposit scorer in left/right detector sets
 *                  (e.g. veto counters and calorimeters)
 *
 *        Version:  1.0
 *        Created:  14.11.2009 12:45:53
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_ENERGY_DEPOSIT_IN_LEFT_RIGHT_SET_HH
#define CEXMC_ENERGY_DEPOSIT_IN_LEFT_RIGHT_SET_HH

#include "CexmcSimpleEnergyDeposit.hh"
#include "CexmcCommon.hh"

class  CexmcSetup;


class  CexmcEnergyDepositInLeftRightSet : public CexmcSimpleEnergyDeposit
{
    public:
        CexmcEnergyDepositInLeftRightSet( const G4String &  name,
                                          const CexmcSetup *  setup );

    public:
        void  PrintAll( void );

    protected:
        G4int  GetIndex( G4Step *  step );

    protected:
        const CexmcSetup *  setup;

    public:
        static CexmcSide  GetSide( G4int  index );

        static G4int  GetLeftRightBitsOffset( void );

    protected:
        static G4int  leftRightBitsOffset;
};


inline CexmcSide  CexmcEnergyDepositInLeftRightSet::GetSide( G4int  index )
{
    if ( index >> leftRightBitsOffset == 1 )
        return CexmcRight;

    return CexmcLeft;
}


inline G4int  CexmcEnergyDepositInLeftRightSet::GetLeftRightBitsOffset( void )
{
    return leftRightBitsOffset;
}


#endif

