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
 *       Filename:  CexmcTrackPoints.hh
 *
 *    Description:  track points collection
 *
 *        Version:  1.0
 *        Created:  16.11.2009 12:41:54
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_TRACK_POINTS_HH
#define CEXMC_TRACK_POINTS_HH

#include <G4THitsMap.hh>
#include "CexmcPrimitiveScorer.hh"
#include "CexmcTrackPointInfo.hh"

class  G4HCofThisEvent;
class  G4Step;


typedef G4THitsMap< CexmcTrackPointInfo >  CexmcTrackPointsCollection;

typedef std::map< G4int, CexmcTrackPointInfo * >
                                           CexmcTrackPointsCollectionData;


class  CexmcTrackPoints : public CexmcPrimitiveScorer
{
    public:
        explicit CexmcTrackPoints( const G4String &  name );

    public:
        void   Initialize( G4HCofThisEvent *  hcOfThisEvent );

        void   EndOfEvent( G4HCofThisEvent *  hcOfThisEvent );

        void   DrawAll( void );

        void   PrintAll( void );

        void   clear( void );

    protected:
        G4int   GetTrackId( G4Step *  step );

        G4int   GetIndex( G4Step *  step );

        G4bool  ProcessHits( G4Step *  step, G4TouchableHistory *  tHistory );

    protected:
        CexmcTrackPointsCollection *  eventMap;

    private:
        G4int                         hcId;
};


#endif

