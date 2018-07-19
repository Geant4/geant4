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
 *       Filename:  CexmcEventAction.hh
 *
 *    Description:  event action
 *
 *        Version:  1.0
 *        Created:  27.10.2009 22:41:29
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_EVENT_ACTION_HH
#define CEXMC_EVENT_ACTION_HH

#include <G4UserEventAction.hh>
#include "CexmcAngularRange.hh"

class  G4Event;
class  CexmcPhysicsManager;
class  CexmcEnergyDepositDigitizer;
struct  CexmcEnergyDepositStore;
class  CexmcTrackPointsDigitizer;
struct  CexmcTrackPointsStore;
class  CexmcEventActionMessenger;
struct  CexmcProductionModelData;
class  CexmcChargeExchangeReconstructor;


class  CexmcEventAction : public G4UserEventAction
{
    public:
        explicit CexmcEventAction( CexmcPhysicsManager *  physicsManager,
                                   G4int  verbose = 0 );

        virtual ~CexmcEventAction();

    public:
        void      BeginOfEventAction( const G4Event *  event );

        void      EndOfEventAction( const G4Event *  event );

    public:
        void      BeamParticleChangeHook( void );

        void      SetVerboseOnCexmcLevel( G4int  value );

        void      SetVerboseDrawLevel( G4int  value );

        CexmcChargeExchangeReconstructor *  GetReconstructor( void );

    private:
        void  PrintReconstructedData(
                        const CexmcAngularRangeList &  angularRanges,
                        const CexmcAngularRange &  angularGap ) const;

#ifdef CEXMC_USE_ROOT
        void  FillEDTHistos( const CexmcEnergyDepositStore *  edStore,
                const CexmcAngularRangeList &  triggeredAngularRanges ) const;

        void  FillTPTHistos( const CexmcTrackPointsStore *  tpStore,
                const CexmcProductionModelData &  pmData,
                const CexmcAngularRangeList &  triggeredAngularRanges ) const;

        void  FillRTHistos( G4bool  reconstructorHasFullTrigger,
                const CexmcEnergyDepositStore *  edStore,
                const CexmcTrackPointsStore *  tpStore,
                const CexmcProductionModelData &  pmData,
                const CexmcAngularRangeList &  triggeredAngularRanges ) const;
#endif

        void  DrawTrajectories( const G4Event *  event );

        void  DrawTrackPoints( const CexmcTrackPointsStore *  tpStore ) const;

        void  DrawReconstructionData( void );

        void  UpdateRunHits( const CexmcAngularRangeList &  aRangesReal,
                             const CexmcAngularRangeList &  aRangesRec,
                             G4bool  tpDigitizerHasTriggered,
                             G4bool  edDigitizerHasTriggered,
                             G4bool  edDigitizerMonitorHasTriggered,
                             G4bool  reconstructorHasTriggered,
                             const CexmcAngularRange &  aGap );

#ifdef CEXMC_USE_PERSISTENCY
        void  SaveEvent( const G4Event *  event,
                         G4bool  edDigitizerMonitorHasTriggered,
                         const CexmcEnergyDepositStore *  edStore,
                         const CexmcTrackPointsStore *  tpStore,
                         const CexmcProductionModelData &  pmData );

        void  SaveEventFast( const G4Event *  event,
                             G4bool  tpDigitizerHasTriggered,
                             G4bool  edDigitizerHasTriggered,
                             G4bool  edDigitizerMonitorHasTriggered,
                             G4double  opCosThetaSCM );
#endif

    public:
        static CexmcEnergyDepositStore *  MakeEnergyDepositStore(
                            const CexmcEnergyDepositDigitizer *  digitizer );

        static CexmcTrackPointsStore *  MakeTrackPointsStore(
                            const CexmcTrackPointsDigitizer *  digitizer );

        static void  PrintEnergyDeposit(
                            const CexmcEnergyDepositStore *  edStore );

        static void  PrintTrackPoints(
                            const CexmcTrackPointsStore *  tpStore );

        static void  PrintProductionModelData(
                            const CexmcAngularRangeList &  angularRanges,
                            const CexmcProductionModelData &  pmData );

    private:
        CexmcPhysicsManager *               physicsManager;

        CexmcChargeExchangeReconstructor *  reconstructor;

#ifdef CEXMC_USE_ROOT
  G4double                            opKinEnergy;
#endif

    private:
        G4int                               verbose;

        G4int                               verboseDraw;

        CexmcEventActionMessenger *         messenger;
};


inline void  CexmcEventAction::SetVerboseOnCexmcLevel( G4int  value )
{
    verbose = value;
}


inline void  CexmcEventAction::SetVerboseDrawLevel( G4int  value )
{
    verboseDraw = value;
}


inline CexmcChargeExchangeReconstructor *
                                    CexmcEventAction::GetReconstructor( void )
{
    return reconstructor;
}


#endif

