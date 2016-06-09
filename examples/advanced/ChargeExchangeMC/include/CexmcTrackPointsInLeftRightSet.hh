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
 *       Filename:  CexmcTrackPointsInLeftRightSet.hh
 *
 *    Description:  track points in left/right detector sets
 *                  (e.g. veto counters and calorimeters)
 *
 *        Version:  1.0
 *        Created:  22.11.2009 21:12:32
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_TRACK_POINTS_IN_LEFT_RIGHT_SET_HH
#define CEXMC_TRACK_POINTS_IN_LEFT_RIGHT_SET_HH

#include "CexmcTrackPoints.hh"
#include "CexmcCommon.hh"

class  CexmcSetup;


class  CexmcTrackPointsInLeftRightSet : public CexmcTrackPoints
{
    public:
        CexmcTrackPointsInLeftRightSet( const G4String &  name,
                                        const CexmcSetup *  setup );

    public:
        void   PrintAll( void );

    protected:
        G4int  GetIndex( G4Step *  step );

    protected:
        const CexmcSetup *  setup;

    public:
        static CexmcSide  GetSide( G4int  index );

        static G4int      GetLeftRightBitsOffset( void );

    protected:
        static const G4int  leftRightBitsOffset = 24;
};


inline CexmcSide  CexmcTrackPointsInLeftRightSet::GetSide( G4int  index )
{
    if ( index >> leftRightBitsOffset == 1 )
        return CexmcRight;

    return CexmcLeft;
}


inline G4int  CexmcTrackPointsInLeftRightSet::GetLeftRightBitsOffset( void )
{
    return leftRightBitsOffset;
}


#endif

