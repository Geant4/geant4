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
 *       Filename:  CexmcRun.hh
 *
 *    Description:  run data (acceptances etc.)
 *
 *        Version:  1.0
 *        Created:  19.12.2009 23:52:51
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_RUN_HH
#define CEXMC_RUN_HH

#include <map>
#include <G4Run.hh>


typedef std::map< G4int, G4int >            CexmcNmbOfHitsInRanges;

typedef CexmcNmbOfHitsInRanges::value_type  CexmcNmbOfHitsInRangesData;


class  CexmcRun : public G4Run
{
    public:
        CexmcRun();

    public:
        void  IncrementNmbOfHitsSampled( G4int  index );

        void  IncrementNmbOfHitsSampledFull( G4int  index );

        void  IncrementNmbOfHitsTriggeredRealRange( G4int  index );

        void  IncrementNmbOfHitsTriggeredRecRange( G4int  index );

        void  IncrementNmbOfOrphanHits( G4int  index );

        void  IncrementNmbOfFalseHitsTriggeredEDT( void );

        void  IncrementNmbOfFalseHitsTriggeredRec( void );

        void  IncrementNmbOfSavedEvents( void );

        void  IncrementNmbOfSavedFastEvents( void );

    public:
        const CexmcNmbOfHitsInRanges &  GetNmbOfHitsSampled( void ) const;

        const CexmcNmbOfHitsInRanges &  GetNmbOfHitsSampledFull( void ) const;

        const CexmcNmbOfHitsInRanges &  GetNmbOfHitsTriggeredRealRange( void )
                                                                        const;

        const CexmcNmbOfHitsInRanges &  GetNmbOfHitsTriggeredRecRange( void )
                                                                        const;

        const CexmcNmbOfHitsInRanges &  GetNmbOfOrphanHits( void ) const;

        G4int  GetNmbOfFalseHitsTriggeredEDT( void ) const;

        G4int  GetNmbOfFalseHitsTriggeredRec( void ) const;

        G4int  GetNmbOfSavedEvents( void ) const;

        G4int  GetNmbOfSavedFastEvents( void ) const;

    private:
        CexmcNmbOfHitsInRanges  nmbOfHitsSampled;

        CexmcNmbOfHitsInRanges  nmbOfHitsSampledFull;

        CexmcNmbOfHitsInRanges  nmbOfHitsTriggeredRealRange;

        CexmcNmbOfHitsInRanges  nmbOfHitsTriggeredRecRange;

        CexmcNmbOfHitsInRanges  nmbOfOrphanHits;

        G4int                   nmbOfFalseHitsTriggeredEDT;

        G4int                   nmbOfFalseHitsTriggeredRec;

        G4int                   nmbOfSavedEvents;

        G4int                   nmbOfSavedFastEvents;
};


inline const CexmcNmbOfHitsInRanges &
                            CexmcRun::GetNmbOfHitsSampled( void ) const
{
    return nmbOfHitsSampled;
}


inline const CexmcNmbOfHitsInRanges &
                            CexmcRun::GetNmbOfHitsSampledFull( void ) const
{
    return nmbOfHitsSampledFull;
}


inline const CexmcNmbOfHitsInRanges &
                        CexmcRun::GetNmbOfHitsTriggeredRealRange( void ) const
{
    return nmbOfHitsTriggeredRealRange;
}


inline const CexmcNmbOfHitsInRanges &
                        CexmcRun::GetNmbOfHitsTriggeredRecRange( void ) const
{
    return nmbOfHitsTriggeredRecRange;
}


inline const CexmcNmbOfHitsInRanges &
                            CexmcRun::GetNmbOfOrphanHits( void ) const
{
    return nmbOfOrphanHits;
}


inline G4int  CexmcRun::GetNmbOfFalseHitsTriggeredEDT( void ) const
{
    return nmbOfFalseHitsTriggeredEDT;
}


inline G4int  CexmcRun::GetNmbOfFalseHitsTriggeredRec( void ) const
{
    return nmbOfFalseHitsTriggeredRec;
}


inline G4int  CexmcRun::GetNmbOfSavedEvents( void ) const
{
    return nmbOfSavedEvents;
}


inline G4int  CexmcRun::GetNmbOfSavedFastEvents( void ) const
{
    return nmbOfSavedFastEvents;
}


#endif

