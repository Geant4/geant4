/*
 * =============================================================================
 *
 *       Filename:  CexmcEnergyDepositInCalorimeter.hh
 *
 *    Description:  energy deposit scorer in calorimeters
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

#ifndef CEXMC_ENERGY_DEPOSIT_IN_CALORIMETER_HH
#define CEXMC_ENERGY_DEPOSIT_IN_CALORIMETER_HH

#include "CexmcEnergyDepositInLeftRightSet.hh"

class  CexmcSetup;


class  CexmcEnergyDepositInCalorimeter : public CexmcEnergyDepositInLeftRightSet
{
    public:
        CexmcEnergyDepositInCalorimeter( const G4String &  name,
                                         const CexmcSetup *  setup );

    public:
        void  PrintAll( void );

    protected:
        G4int  GetIndex( G4Step *  step );

    public:
        static G4int  GetRow( G4int  index );

        static G4int  GetColumn( G4int  index );

        static G4int  GetCopyDepth1BitsOffset( void );

    protected:
        static G4int  copyDepth1BitsOffset;
};


inline G4int  CexmcEnergyDepositInCalorimeter::GetRow( G4int  index )
{
    index &= ( ( 1 << ( leftRightBitsOffset - 1 ) ) |
                          ( ( 1 << ( leftRightBitsOffset - 1 ) ) - 1 ) );
    return index >> copyDepth1BitsOffset;
}


inline G4int  CexmcEnergyDepositInCalorimeter::GetColumn( G4int  index )
{
    return index & ( ( 1 << ( copyDepth1BitsOffset - 1 ) ) |
                              ( ( 1 << ( copyDepth1BitsOffset - 1 ) ) - 1 ) );
}


inline G4int  CexmcEnergyDepositInCalorimeter::GetCopyDepth1BitsOffset( void )
{
    return copyDepth1BitsOffset;
}


#endif

