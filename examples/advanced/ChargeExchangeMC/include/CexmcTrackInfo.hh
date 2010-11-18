/*
 * =============================================================================
 *
 *       Filename:  CexmcTrackInfo.hh
 *
 *    Description:  track info
 *
 *        Version:  1.0
 *        Created:  22.11.2009 18:42:40
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_TRACK_INFO_HH
#define CEXMC_TRACK_INFO_HH

#include <G4VUserTrackInformation.hh>
#include "CexmcCommon.hh"


class  CexmcTrackInfo : public G4VUserTrackInformation
{
    public:
        explicit CexmcTrackInfo( CexmcTrackType  trackType = CexmcInsipidTrack,
                                 G4int copyNumber = 0 );

    public:
        void  Print( void ) const;

    public:
        virtual G4int   GetTypeInfo( void ) const;

    public:
        CexmcTrackType  GetTrackType( void ) const;

        void            SetTrackType( CexmcTrackType  value );

        G4int           GetCopyNumber( void ) const;

    private:
        CexmcTrackType  trackType;

        G4int           copyNumber;
};


inline CexmcTrackType  CexmcTrackInfo::GetTrackType( void ) const
{
    return trackType;
}


inline void  CexmcTrackInfo::SetTrackType( CexmcTrackType  value )
{
    trackType = value;
}


inline G4int  CexmcTrackInfo::GetCopyNumber( void ) const
{
    return copyNumber;
}


#endif

