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
  // explicit argument is only needed by G4THitsMap template,
  // it has no actual use here
  explicit CexmcTrackPointInfo( G4double  unused = 0. ) :
    momentumAmp ( unused ), trackId( CexmcInvalidTrackId )		       
  {}

    CexmcTrackPointInfo( const G4ThreeVector &  positionLocal_,
                         const G4ThreeVector &  positionWorld_,
                         const G4ThreeVector &  directionLocal_,
                         const G4ThreeVector &  directionWorld_,
                         G4double  momentumAmp_,
                         const G4ParticleDefinition *  particle_,
                         G4int  trackId_, CexmcTrackType  trackType_ ) :
        positionLocal( positionLocal_ ), positionWorld( positionWorld_ ),
        directionLocal( directionLocal_ ), directionWorld( directionWorld_ ),
        momentumAmp( momentumAmp_ ), particle( particle_ ), trackId( trackId_ ),
        trackType( trackType_ )
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

