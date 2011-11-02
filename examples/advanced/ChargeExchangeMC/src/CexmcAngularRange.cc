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
 * ============================================================================
 *
 *       Filename:  CexmcAngularRange.cc
 *
 *    Description:  angular range object
 *
 *        Version:  1.0
 *        Created:  28.12.2009 22:35:07
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "CexmcAngularRange.hh"


void  GetNormalizedAngularRange( const CexmcAngularRangeList &  src,
                                 CexmcAngularRangeList &  dst )
{
    dst = src;
    if ( dst.size() < 2 )
        return;

    std::sort( dst.begin(), dst.end() );

    const G4double  epsilon( 1E-7 );

    for ( CexmcAngularRangeList::iterator  k( dst.begin() + 1 );
                                                            k != dst.end(); )
    {
        if ( std::fabs( k->top - ( k - 1 )->top ) < epsilon ||
             k->bottom + epsilon >= ( k - 1 )->bottom )
        {
            k = dst.erase( k );
            continue;
        }
        if ( k->top + epsilon >= ( k - 1 )->bottom )
        {
            ( k - 1 )->bottom = k->bottom;
            k = dst.erase( k );
            continue;
        }
        ++k;
    }
}


void  GetAngularGaps( const CexmcAngularRangeList &  src,
                      CexmcAngularRangeList &  dst )
{
    if ( src.empty() )
    {
        dst.push_back( CexmcAngularRange( 1.0, -1.0 , 0 ) );
        return;
    }

    CexmcAngularRangeList  normalizedAngularRanges;
    GetNormalizedAngularRange( src, normalizedAngularRanges );

    G4int  index( 0 );
    if ( normalizedAngularRanges[ 0 ].top < 1.0 )
        dst.push_back( CexmcAngularRange(
                        1.0, normalizedAngularRanges[ 0 ].top, index++ ) );

    for ( CexmcAngularRangeList::iterator
            k( normalizedAngularRanges.begin() );
            k != normalizedAngularRanges.end(); ++k )
    {
        if ( k + 1 == normalizedAngularRanges.end() )
            break;
        dst.push_back( CexmcAngularRange(
                        k->bottom, ( k + 1 )->top, index++ ) );
    }

    if ( normalizedAngularRanges.back().bottom > -1.0 )
        dst.push_back( CexmcAngularRange(
                        normalizedAngularRanges.back().bottom, -1.0, index ) );
}


std::ostream &  operator<<( std::ostream &  out,
                            const CexmcAngularRange &  angularRange )
{
    std::ostream::fmtflags  savedFlags( out.flags() );
    std::streamsize         prec( out.precision() );

    out.precision( 4 );
    out.flags( std::ios::fixed );

    out << std::setw( 2 ) << angularRange.index + 1 << " [" << std::setw( 7 ) <<
           angularRange.top << ", " << std::setw( 7 ) << angularRange.bottom <<
           ")";

    out.precision( prec );
    out.flags( savedFlags );

    return out;
}


std::ostream &  operator<<( std::ostream &  out,
                            const CexmcAngularRangeList &  angularRanges )
{
    out << std::endl;
    for ( CexmcAngularRangeList::const_iterator  k( angularRanges.begin() );
                                                k != angularRanges.end(); ++k )
    {
        out << "                 " << *k << std::endl;
    }

    return out;
}

