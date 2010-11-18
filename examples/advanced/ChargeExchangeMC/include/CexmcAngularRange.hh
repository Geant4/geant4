/*
 * =============================================================================
 *
 *       Filename:  CexmcAngularRange.hh
 *
 *    Description:  angular range object
 *
 *        Version:  1.0
 *        Created:  01.12.2009 16:29:25
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_ANGULAR_RANGE_HH
#define CEXMC_ANGULAR_RANGE_HH

#include <vector>
#include <iosfwd>
#include <G4Types.hh>


struct  CexmcAngularRange
{
    CexmcAngularRange()
    {}

    CexmcAngularRange( G4double  top, G4double  bottom, G4int  index ) :
        top( top ), bottom( bottom ), index( index )
    {}

    G4double  top;

    G4double  bottom;

    G4int     index;

    template  < typename  Archive >
    void  serialize( Archive &  archive, const unsigned int  version );
};


typedef std::vector< CexmcAngularRange >  CexmcAngularRangeList;


template  < typename  Archive >
inline void  CexmcAngularRange::serialize( Archive &  archive,
                                           const unsigned int )
{
    archive & top;
    archive & bottom;
    archive & index;
}


inline bool  operator<( const CexmcAngularRange &  left,
                        const CexmcAngularRange &  right )
{
    if ( left.top != right.top )
        return left.top > right.top;
    if ( left.bottom != right.bottom )
        return left.bottom < right.bottom;

    return false;
}


void  GetNormalizedAngularRange( const CexmcAngularRangeList &  src,
                                 CexmcAngularRangeList &  dst );


void  GetAngularGaps( const CexmcAngularRangeList &  src,
                      CexmcAngularRangeList &  dst );


std::ostream &  operator<<( std::ostream &  out,
                            const CexmcAngularRange &  angularRange );


std::ostream &  operator<<( std::ostream &  out,
                            const CexmcAngularRangeList &  angularRanges );


#endif

