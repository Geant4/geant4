/*
 * =============================================================================
 *
 *       Filename:  CexmcSimpleTrackPointInfoStore.hh
 *
 *    Description:  serialization helper for track point info objects
 *
 *        Version:  1.0
 *        Created:  31.12.2009 13:55:51
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_SIMPLE_TRACK_POINT_INFO_STORE_HH
#define CEXMC_SIMPLE_TRACK_POINT_INFO_STORE_HH

#ifdef CEXMC_USE_PERSISTENCY

#include <boost/serialization/access.hpp>
#include "CexmcSimpleThreeVectorStore.hh"
#include "CexmcCommon.hh"

class  CexmcTrackPointInfo;


class  CexmcSimpleTrackPointInfoStore
{
    friend class  boost::serialization::access;
#ifdef CEXMC_USE_CUSTOM_FILTER
    friend class  CexmcASTEval;
#endif

    public:
        CexmcSimpleTrackPointInfoStore();

        CexmcSimpleTrackPointInfoStore( const CexmcTrackPointInfo &  tpInfo );

    public:
        operator CexmcTrackPointInfo() const;

    private:
        template  < typename  Archive >
        void  serialize( Archive &  archive, const unsigned int  version );

    private:
        CexmcSimpleThreeVectorStore  positionLocal;

        CexmcSimpleThreeVectorStore  positionWorld;

        CexmcSimpleThreeVectorStore  directionLocal;

        CexmcSimpleThreeVectorStore  directionWorld;

        G4double                     momentumAmp;

        G4int                        particlePDGEncoding;

        G4int                        trackId;

        CexmcTrackType               trackType;
};


template  < typename  Archive >
void  CexmcSimpleTrackPointInfoStore::serialize( Archive &  archive,
                                                 const unsigned int )
{
    archive & positionLocal;
    archive & positionWorld;
    archive & directionLocal;
    archive & directionWorld;
    archive & momentumAmp;
    archive & particlePDGEncoding;
    archive & trackId;
    archive & trackType;
}

#endif

#endif

