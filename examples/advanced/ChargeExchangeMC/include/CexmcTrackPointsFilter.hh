/*
 * =============================================================================
 *
 *       Filename:  CexmcTrackPointsFilter.hh
 *
 *    Description:  track points of interest
 *
 *        Version:  1.0
 *        Created:  16.11.2009 22:23:00
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_TRACK_POINTS_FILTER_HH
#define CEXMC_TRACK_POINTS_FILTER_HH

#include <G4VSDFilter.hh>

class  G4String;
class  G4Step;


class  CexmcTrackPointsFilter : public G4VSDFilter
{
    public:
        explicit CexmcTrackPointsFilter( const G4String &  name );

    public:
        G4bool  Accept( const G4Step *  step ) const;
};


#endif

