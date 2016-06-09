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
 *       Filename:  CexmcTrackPointsInCalorimeter.hh
 *
 *    Description:  track points in calorimeter
 *
 *        Version:  1.0
 *        Created:  22.11.2009 21:57:15
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_TRACK_POINTS_IN_CALORIMETER_HH
#define CEXMC_TRACK_POINTS_IN_CALORIMETER_HH

#include "CexmcTrackPointsInLeftRightSet.hh"

class  CexmcSetup;


class  CexmcTrackPointsInCalorimeter : public CexmcTrackPointsInLeftRightSet
{
    public:
        CexmcTrackPointsInCalorimeter( const G4String &  name,
                                       const CexmcSetup *  setup );

    public:
        void   PrintAll( void );

    protected:
        G4int  GetIndex( G4Step *  step );

    public:
        static G4int  GetRow( G4int  index );

        static G4int  GetColumn( G4int  index );

        static G4int  GetCopyDepth0BitsOffset( void );

        static G4int  GetCopyDepth1BitsOffset( void );

    protected:
        static const G4int  copyDepth0BitsOffset = 8;

        static const G4int  copyDepth1BitsOffset = 16;
};


inline G4int  CexmcTrackPointsInCalorimeter::GetRow( G4int  index )
{
    index &= ( ( 1 << ( leftRightBitsOffset - 1 ) ) |
                          ( ( 1 << ( leftRightBitsOffset - 1 ) ) - 1 ) );
    return index >> copyDepth1BitsOffset;
}


inline G4int  CexmcTrackPointsInCalorimeter::GetColumn( G4int  index )
{
    index &= ( ( 1 << ( copyDepth1BitsOffset - 1 ) ) |
                          ( ( 1 << ( copyDepth1BitsOffset - 1 ) ) - 1 ) );
    return index >> copyDepth0BitsOffset;
}


inline G4int  CexmcTrackPointsInCalorimeter::GetCopyDepth0BitsOffset( void )
{
    return copyDepth0BitsOffset;
}


inline G4int  CexmcTrackPointsInCalorimeter::GetCopyDepth1BitsOffset( void )
{
    return copyDepth1BitsOffset;
}


#endif

