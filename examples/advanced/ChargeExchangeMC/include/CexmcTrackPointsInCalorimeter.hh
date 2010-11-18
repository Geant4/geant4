/*
 * =============================================================================
 *
 *       Filename:  CexmcTrackPointsInCalorimeter.hh
 *
 *    Description:  track points in calorimeter
 *
 *        Version:  1.0
 *        Created:  22.11.2009 21:57:15
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_TRACK_POINTS_IN_CALORIMETER_HH
#define CEXMC_TRACK_POINTS_IN_CALORIMETER_HH

#include "CexmcTrackPointsInLeftRightSet.hh"

class  CexmcSetup;


class  CexmcTrackPointsInCalorimeter : public CexmcTrackPointsInLeftRightSet
{
    public:
        CexmcTrackPointsInCalorimeter( const G4String &  name,
                                       const CexmcSetup *  setup );

    public:
        void  PrintAll( void );

    protected:
        G4int  GetIndex( G4Step *  step );

    public:
        static G4int  GetRow( G4int  index );

        static G4int  GetColumn( G4int  index );

        static G4int  GetCopyDepth0BitsOffset( void );

        static G4int  GetCopyDepth1BitsOffset( void );

    protected:
        static G4int  copyDepth0BitsOffset;

        static G4int  copyDepth1BitsOffset;
};


inline G4int  CexmcTrackPointsInCalorimeter::GetRow( G4int  index )
{
    index &= ( ( 1 << ( leftRightBitsOffset - 1 ) ) |
                          ( ( 1 << ( leftRightBitsOffset - 1 ) ) - 1 ) );
    return index >> copyDepth1BitsOffset;
}


inline G4int  CexmcTrackPointsInCalorimeter::GetColumn( G4int  index )
{
    index &= ( ( 1 << ( copyDepth1BitsOffset - 1 ) ) |
                          ( ( 1 << ( copyDepth1BitsOffset - 1 ) ) - 1 ) );
    return index >> copyDepth0BitsOffset;
}


inline G4int  CexmcTrackPointsInCalorimeter::GetCopyDepth0BitsOffset( void )
{
    return copyDepth0BitsOffset;
}


inline G4int  CexmcTrackPointsInCalorimeter::GetCopyDepth1BitsOffset( void )
{
    return copyDepth1BitsOffset;
}


#endif

