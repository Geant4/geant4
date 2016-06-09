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
 *       Filename:  CexmcTrackInfo.hh
 *
 *    Description:  track info
 *
 *        Version:  1.0
 *        Created:  22.11.2009 18:42:40
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_TRACK_INFO_HH
#define CEXMC_TRACK_INFO_HH

#include <G4VUserTrackInformation.hh>
#include "CexmcCommon.hh"


class  CexmcTrackInfo : public G4VUserTrackInformation
{
    public:
        explicit CexmcTrackInfo( CexmcTrackType  trackType = CexmcInsipidTrack,
                                 G4int  copyNumber = 0 );

    public:
        void            Print( void ) const;

    public:
        virtual G4int   GetTypeInfo( void ) const;

    public:
        CexmcTrackType  GetTrackType( void ) const;

        void            SetTrackType( CexmcTrackType  value );

        G4int           GetCopyNumber( void ) const;

    private:
        CexmcTrackType  trackType;

        G4int           copyNumber;
};


inline CexmcTrackType  CexmcTrackInfo::GetTrackType( void ) const
{
    return trackType;
}


inline void  CexmcTrackInfo::SetTrackType( CexmcTrackType  value )
{
    trackType = value;
}


inline G4int  CexmcTrackInfo::GetCopyNumber( void ) const
{
    return copyNumber;
}


#endif

