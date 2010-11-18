/*
 * =============================================================================
 *
 *       Filename:  CexmcSimpleRangeWithValue.cc
 *
 *    Description:  auxiliary functions for simple range instances
 *
 *        Version:  1.0
 *        Created:  17.02.2010 22:46:01
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#include <iomanip>
#include <G4UnitsTable.hh>
#include "CexmcSimpleRangeWithValue.hh"


std::ostream &  operator<<( std::ostream &  out,
                        const CexmcEnergyRangeWithDoubleValue &  range )
{
    out << "[ " << std::setw( 3 ) << G4BestUnit( range.bottom, "Energy" ) <<
           ", " << std::setw( 3 ) << G4BestUnit( range.top, "Energy" ) <<
           " )   " << range.value;

    return out;
}


std::ostream &  operator<<( std::ostream &  out,
                        const CexmcEnergyRangeWithDoubleValueList &  ranges )
{
    out << std::endl;
    for ( CexmcEnergyRangeWithDoubleValueList::const_iterator
                                  k( ranges.begin() ); k != ranges.end(); ++k )
    {
        out << "          " << *k << std::endl;
    }

    return out;
}

