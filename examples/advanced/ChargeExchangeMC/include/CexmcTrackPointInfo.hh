/*
 * =============================================================================
 *
 *       Filename:  CexmcTrackPointInfo.hh
 *
 *    Description:  single track point information
 *
 *        Version:  1.0
 *        Created:  16.11.2009 12:51:50
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_TRACK_POINT_INFO_HH
#define CEXMC_TRACK_POINT_INFO_HH

#include <iosfwd>
#include <G4ThreeVector.hh>
#include <G4ParticleDefinition.hh>
#include <G4Allocator.hh>
#include <G4UnitsTable.hh>
#include <G4Types.hh>
#include "CexmcCommon.hh"


struct  CexmcTrackPointInfo
{
    CexmcTrackPointInfo() : trackId( CexmcInvalidTrackId )
    {}

    CexmcTrackPointInfo( const G4ThreeVector &  positionLocal,
                         const G4ThreeVector &  positionWorld,
                         const G4ThreeVector &  directionLocal,
                         const G4ThreeVector &  directionWorld,
                         G4double  momentumAmp,
                         const G4ParticleDefinition *  particle,
                         G4int  trackId, CexmcTrackType  trackType ) :
        positionLocal( positionLocal ), positionWorld( positionWorld ),
        directionLocal( directionLocal ), directionWorld( directionWorld ),
        momentumAmp( momentumAmp ), particle( particle ), trackId( trackId ),
        trackType( trackType )
    {}

    G4bool  IsValid( void ) const
    {
        return trackId != CexmcInvalidTrackId;
    }

    void *  operator new( size_t  size );

    void    operator delete( void *  obj );

    G4ThreeVector                 positionLocal;

    G4ThreeVector                 positionWorld;

    G4ThreeVector                 directionLocal;

    G4ThreeVector                 directionWorld;

    G4double                      momentumAmp;

    const G4ParticleDefinition *  particle;

    G4int                         trackId;

    CexmcTrackType                trackType;

    // following type cast operator is only needed by G4THitsMap template
    // (in PrintAll()), it has no actual use here
    operator G4double()
    {
        return 0;
    }
};


extern G4Allocator< CexmcTrackPointInfo >  trackPointInfoAllocator;


inline void *  CexmcTrackPointInfo::operator new( size_t )
{
  return trackPointInfoAllocator.MallocSingle();
}


inline void  CexmcTrackPointInfo::operator delete( void *  obj )
{
    trackPointInfoAllocator.FreeSingle(
                            reinterpret_cast< CexmcTrackPointInfo * >( obj ) );
}


std::ostream &  operator<<( std::ostream &  out,
                            const CexmcTrackPointInfo &  trackPointInfo );

#endif

