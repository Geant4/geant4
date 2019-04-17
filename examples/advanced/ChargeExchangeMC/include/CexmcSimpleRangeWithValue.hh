//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
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
#include <iosfwd>
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

    CexmcSimpleRangeWithValue( G4double  bottom_, G4double  top_,
                               G4double  value_ ) :
        bottom( bottom_ ), top( top_ ), value( value_ )
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
inline G4bool  operator<(
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

