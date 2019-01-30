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

    CexmcAngularRange( G4double  top_, G4double  bottom_, G4int  index_ ) :
        top( top_ ), bottom( bottom_ ), index( index_ )
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


inline G4bool  operator<( const CexmcAngularRange &  left,
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

