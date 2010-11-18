/*
 * =============================================================================
 *
 *       Filename:  CexmcTrackPointsInLeftRightSet.hh
 *
 *    Description:  track points in left/right detector sets
 *                  (e.g. veto counters and calorimeters)
 *
 *        Version:  1.0
 *        Created:  22.11.2009 21:12:32
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_TRACK_POINTS_IN_LEFT_RIGHT_SET_HH
#define CEXMC_TRACK_POINTS_IN_LEFT_RIGHT_SET_HH

#include "CexmcTrackPoints.hh"
#include "CexmcCommon.hh"

class  CexmcSetup;


class  CexmcTrackPointsInLeftRightSet : public CexmcTrackPoints
{
    public:
        CexmcTrackPointsInLeftRightSet( const G4String &  name,
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


inline CexmcSide  CexmcTrackPointsInLeftRightSet::GetSide( G4int  index )
{
    if ( index >> leftRightBitsOffset == 1 )
        return CexmcRight;

    return CexmcLeft;
}


inline G4int  CexmcTrackPointsInLeftRightSet::GetLeftRightBitsOffset( void )
{
    return leftRightBitsOffset;
}


#endif

