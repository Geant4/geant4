/*
 * =============================================================================
 *
 *       Filename:  CexmcSimpleRangeWithValue.hh
 *
 *    Description:  simple range with value (can be serialized)
 *
 *        Version:  1.0
 *        Created:  17.02.2010 14:47:55
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_SIMPLE_RANGE_WITH_VALUE_HH
#define CEXMC_SIMPLE_RANGE_WITH_VALUE_HH

#include <vector>
#include <iostream>
#include <G4Types.hh>


enum  CexmcValueCategory
{
    CexmcPlainValueCategory,
    CexmcEnergyValueCategory
};


template  < CexmcValueCategory  RangeCategory = CexmcPlainValueCategory,
            CexmcValueCategory  ValueCategory = CexmcPlainValueCategory >
struct  CexmcSimpleRangeWithValue
{
    CexmcSimpleRangeWithValue()
    {}

    CexmcSimpleRangeWithValue( G4double  bottom, G4double  top,
                               G4double  value ) :
        bottom( bottom ), top( top ), value( value )
    {}

    G4double  bottom;

    G4double  top;

    G4double  value;

    template  < typename  Archive >
    void  serialize( Archive &  archive, const unsigned int  version );
};


template  < CexmcValueCategory  RangeCategory,
            CexmcValueCategory  ValueCategory >
template  < typename  Archive >
inline void
        CexmcSimpleRangeWithValue< RangeCategory, ValueCategory >::serialize(
                                        Archive &  archive, const unsigned int )
{
    archive & bottom;
    archive & top;
    archive & value;
}


template  < CexmcValueCategory  RangeCategory,
            CexmcValueCategory  ValueCategory >
inline bool  operator<(
    const CexmcSimpleRangeWithValue< RangeCategory, ValueCategory > &  left,
    const CexmcSimpleRangeWithValue< RangeCategory, ValueCategory > &  right )
{
    if ( left.bottom != right.bottom )
        return left.bottom < right.bottom;
    if ( left.top != right.top )
        return left.top > right.top;

    return false;
}


typedef CexmcSimpleRangeWithValue< CexmcEnergyValueCategory >
                                                CexmcEnergyRangeWithDoubleValue;
typedef std::vector< CexmcEnergyRangeWithDoubleValue >
                                            CexmcEnergyRangeWithDoubleValueList;


std::ostream &  operator<<( std::ostream &  out,
                        const CexmcEnergyRangeWithDoubleValue &  range );

std::ostream &  operator<<( std::ostream &  out,
                        const CexmcEnergyRangeWithDoubleValueList &  range );


#endif

