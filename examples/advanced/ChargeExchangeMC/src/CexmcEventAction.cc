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
 *       Filename:  CexmcEventAction.cc
 *
 *    Description:  event action
 *
 *        Version:  1.0
 *        Created:  27.10.2009 22:48:08
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <cmath>
#ifdef CEXMC_USE_PERSISTENCY
#include <boost/archive/binary_oarchive.hpp>
#endif
#include <G4DigiManager.hh>
#include <G4Event.hh>
#include <G4Circle.hh>
#include <G4VisAttributes.hh>
#include <G4VisManager.hh>
#include <G4VTrajectory.hh>
#include <G4TrajectoryContainer.hh>
#include <G4Colour.hh>
#include <G4SystemOfUnits.hh>
#include "CexmcEventAction.hh"
#include "CexmcEventActionMessenger.hh"
#include "CexmcEventInfo.hh"
#include "CexmcEventSObject.hh"
#include "CexmcEventFastSObject.hh"
#include "CexmcTrackingAction.hh"
#include "CexmcChargeExchangeReconstructor.hh"
#include "CexmcRunManager.hh"
#include "CexmcHistoManager.hh"
#include "CexmcRun.hh"
#include "CexmcPhysicsManager.hh"
#include "CexmcProductionModel.hh"
#include "CexmcProductionModelData.hh"
#include "CexmcEnergyDepositDigitizer.hh"
#include "CexmcEnergyDepositStore.hh"
#include "CexmcTrackPointsDigitizer.hh"
#include "CexmcTrackPointsStore.hh"
#include "CexmcTrackPointInfo.hh"
#include "CexmcException.hh"
#include "CexmcCommon.hh"


namespace
{
    G4double  CexmcSmallCircleScreenSize( 5.0 );
    G4double  CexmcBigCircleScreenSize( 10.0 );
    G4Colour  CexmcTrackPointsMarkerColour( 0.0, 1.0, 0.4 );
    G4Colour  CexmcRecTrackPointsMarkerColour( 1.0, 0.4, 0.0 );

#ifdef CEXMC_USE_ROOT
    inline G4double  CexmcGetKinEnergy( G4double  momentumAmp, G4double  mass )
    {
        return std::sqrt( momentumAmp * momentumAmp + mass * mass ) - mass;
    }
#endif
}


CexmcEventAction::CexmcEventAction( CexmcPhysicsManager *  physicsManager_,
                                    G4int  verbose_ ) :
    physicsManager( physicsManager_ ), reconstructor( NULL ),
#ifdef CEXMC_USE_ROOT
    opKinEnergy( 0. ),
#endif
    verbose( verbose_ ), verboseDraw( 4 ), messenger( NULL )
{
    G4DigiManager *  digiManager( G4DigiManager::GetDMpointer() );
    digiManager->AddNewModule( new CexmcEnergyDepositDigitizer(
                                                    CexmcEDDigitizerName ) );
    digiManager->AddNewModule( new CexmcTrackPointsDigitizer(
                                                    CexmcTPDigitizerName ) );
    reconstructor = new CexmcChargeExchangeReconstructor(
                                        physicsManager->GetProductionModel() );
    messenger = new CexmcEventActionMessenger( this );
}


CexmcEventAction::~CexmcEventAction()
{
    delete reconstructor;
    delete messenger;
}


void  CexmcEventAction::BeamParticleChangeHook( void )
{
    reconstructor->SetupBeamParticle();
}


void  CexmcEventAction::BeginOfEventAction( const G4Event * )
{
    G4RunManager *         runManager( G4RunManager::GetRunManager() );
    CexmcTrackingAction *  trackingAction
            ( static_cast< CexmcTrackingAction * >(
                        const_cast< G4UserTrackingAction * >(
                                    runManager->GetUserTrackingAction() ) ) );
    trackingAction->BeginOfEventAction();

    physicsManager->ResetNumberOfTriggeredStudiedInteractions();
}


CexmcEnergyDepositStore *  CexmcEventAction::MakeEnergyDepositStore(
                                const CexmcEnergyDepositDigitizer *  digitizer )
{
    G4double  monitorED( digitizer->GetMonitorED() );
    G4double  vetoCounterEDLeft( digitizer->GetVetoCounterEDLeft() );
    G4double  vetoCounterEDRight( digitizer->GetVetoCounterEDRight() );
    G4double  calorimeterEDLeft( digitizer->GetCalorimeterEDLeft() );
    G4double  calorimeterEDRight( digitizer->GetCalorimeterEDRight() );
    G4int     calorimeterEDLeftMaxX( digitizer->GetCalorimeterEDLeftMaxX() );
    G4int     calorimeterEDLeftMaxY( digitizer->GetCalorimeterEDLeftMaxY() );
    G4int     calorimeterEDRightMaxX( digitizer->GetCalorimeterEDRightMaxX() );
    G4int     calorimeterEDRightMaxY( digitizer->GetCalorimeterEDRightMaxY() );

    const CexmcEnergyDepositCalorimeterCollection &
              calorimeterEDLeftCollection(
                            digitizer->GetCalorimeterEDLeftCollection() );
    const CexmcEnergyDepositCalorimeterCollection &
              calorimeterEDRightCollection(
                            digitizer->GetCalorimeterEDRightCollection() );

    /* ATTENTION: return object in heap - must be freed by caller! */
    return new CexmcEnergyDepositStore( monitorED, vetoCounterEDLeft,
                    vetoCounterEDRight, calorimeterEDLeft, calorimeterEDRight,
                    calorimeterEDLeftMaxX, calorimeterEDLeftMaxY,
                    calorimeterEDRightMaxX, calorimeterEDRightMaxY,
                    calorimeterEDLeftCollection, calorimeterEDRightCollection );
}


CexmcTrackPointsStore *  CexmcEventAction::MakeTrackPointsStore(
                                const CexmcTrackPointsDigitizer *  digitizer )
{
    const CexmcTrackPointInfo &
                monitorTP( digitizer->GetMonitorTP() );
    const CexmcTrackPointInfo &
                targetTPBeamParticle(
                    digitizer->GetTargetTPBeamParticle() );
    const CexmcTrackPointInfo &
                targetTPOutputParticle(
                    digitizer->GetTargetTPOutputParticle() );
    const CexmcTrackPointInfo &
                targetTPNucleusParticle(
                    digitizer->GetTargetTPNucleusParticle() );
    const CexmcTrackPointInfo &
                targetTPOutputParticleDecayProductParticle1(
                    digitizer->
                    GetTargetTPOutputParticleDecayProductParticle( 0 ) );
    const CexmcTrackPointInfo &
                targetTPOutputParticleDecayProductParticle2(
                    digitizer->
                    GetTargetTPOutputParticleDecayProductParticle( 1 ) );
    const CexmcTrackPointInfo &
                vetoCounterTPLeft(
                    digitizer->GetVetoCounterTPLeft() );
    const CexmcTrackPointInfo &
                vetoCounterTPRight(
                    digitizer->GetVetoCounterTPRight() );
    const CexmcTrackPointInfo &
                calorimeterTPLeft(
                    digitizer->GetCalorimeterTPLeft() );
    const CexmcTrackPointInfo &
                calorimeterTPRight(
                    digitizer->GetCalorimeterTPRight() );

    /* ATTENTION: return object in heap - must be freed by caller! */
    return new CexmcTrackPointsStore( monitorTP, targetTPBeamParticle,
                  targetTPOutputParticle, targetTPNucleusParticle,
                  targetTPOutputParticleDecayProductParticle1,
                  targetTPOutputParticleDecayProductParticle2,
                  vetoCounterTPLeft, vetoCounterTPRight,
                  calorimeterTPLeft, calorimeterTPRight );
}


void  CexmcEventAction::PrintEnergyDeposit(
                                    const CexmcEnergyDepositStore *  edStore )
{
    G4cout << " --- Energy Deposit" << G4endl;
    G4cout << "       monitor : " <<
            G4BestUnit( edStore->monitorED, "Energy" ) << G4endl;
    G4cout << "        vc (l) : " <<
            G4BestUnit( edStore->vetoCounterEDLeft, "Energy" ) << G4endl;
    G4cout << "        vc (r) : " <<
            G4BestUnit( edStore->vetoCounterEDRight, "Energy" ) << G4endl;
    G4cout << "       cal (l) : " <<
            G4BestUnit( edStore->calorimeterEDLeft, "Energy" );
    G4cout << edStore->calorimeterEDLeftCollection;
    G4cout << "       cal (r) : " <<
            G4BestUnit( edStore->calorimeterEDRight, "Energy" );
    G4cout << edStore->calorimeterEDRightCollection;
}


void  CexmcEventAction::PrintTrackPoints(
                                    const CexmcTrackPointsStore *  tpStore )
{
    if ( ! tpStore )
        return;

    G4cout << " --- Track Points" << G4endl;
    G4cout << "       monitor : " << tpStore->monitorTP << G4endl;
    G4cout << "        target : " << tpStore->targetTPBeamParticle << G4endl;
    G4cout << "               : " << tpStore->targetTPOutputParticle << G4endl;
    G4cout << "               : " << tpStore->targetTPNucleusParticle << G4endl;
    G4cout << "               : " <<
            tpStore->targetTPOutputParticleDecayProductParticle1 << G4endl;
    G4cout << "               : " <<
            tpStore->targetTPOutputParticleDecayProductParticle2 << G4endl;
    G4cout << "        vc (l) : " << tpStore->vetoCounterTPLeft << G4endl;
    G4cout << "        vc (r) : " << tpStore->vetoCounterTPRight << G4endl;
    G4cout << "       cal (l) : " << tpStore->calorimeterTPLeft << G4endl;
    G4cout << "       cal (r) : " << tpStore->calorimeterTPRight << G4endl;
    G4cout << "      ---" << G4endl;
    G4cout << "      angle between the " <<
        tpStore->targetTPOutputParticle.particle->GetParticleName() <<
        " decay products : " <<
        tpStore->targetTPOutputParticleDecayProductParticle1.directionWorld.
            angle(
        tpStore->targetTPOutputParticleDecayProductParticle2.directionWorld ) /
            deg << " deg" << G4endl;
}


void  CexmcEventAction::PrintProductionModelData(
                                const CexmcAngularRangeList &  angularRanges,
                                const CexmcProductionModelData &  pmData )
{
    G4cout << " --- Triggered angular ranges: " << angularRanges;
    G4cout << " --- Production model data: " << pmData;
}


void  CexmcEventAction::PrintReconstructedData(
                    const CexmcAngularRangeList &  triggeredRecAngularRanges,
                    const CexmcAngularRange &  angularGap ) const
{
    G4cout << " --- Reconstructed data: " << G4endl;
    G4cout << "       -- entry points:" << G4endl;
    G4cout << "              left: " << G4BestUnit(
        reconstructor->GetCalorimeterEPLeftPosition(), "Length" ) << G4endl;
    G4cout << "             right: " << G4BestUnit(
        reconstructor->GetCalorimeterEPRightPosition(), "Length" ) << G4endl;
    G4cout << "            target: " << G4BestUnit(
        reconstructor->GetTargetEPPosition(), "Length" ) << G4endl;
    G4cout << "       -- the angle: " << reconstructor->GetTheAngle() / deg <<
        " deg" << G4endl;
    G4cout << "       -- mass of the output particle: " << G4BestUnit(
        reconstructor->GetOutputParticleMass(), "Energy" ) << G4endl;
    G4cout << "       -- mass of the nucleus output particle: " << G4BestUnit(
        reconstructor->GetNucleusOutputParticleMass(), "Energy" ) << G4endl;
    if ( reconstructor->IsMassCutUsed() )
    {
        if ( reconstructor->HasMassCutTriggered() )
            G4cout << "            < mass cut passed >" << G4endl;
        else
            G4cout << "            < mass cut failed >" << G4endl;
    }
    if ( reconstructor->IsAbsorbedEnergyCutUsed() )
    {
        if ( reconstructor->HasAbsorbedEnergyCutTriggered() )
            G4cout << "            < absorbed energy cut passed >" << G4endl;
        else
            G4cout << "            < absorbed energy cut failed >" << G4endl;
    }
    const CexmcProductionModelData &  pmData(
                                    reconstructor->GetProductionModelData() );
    G4cout << "       -- production model data: " << pmData;
    G4cout << "       -- triggered angular ranges: ";
    if ( triggeredRecAngularRanges.empty() )
        G4cout << "< orphan detected, gap " << angularGap << " >" << G4endl;
    else
        G4cout << triggeredRecAngularRanges;
}


#ifdef CEXMC_USE_ROOT

void  CexmcEventAction::FillEDTHistos( const CexmcEnergyDepositStore *  edStore,
                const CexmcAngularRangeList &  triggeredAngularRanges ) const
{
    CexmcHistoManager *  histoManager( CexmcHistoManager::Instance() );

    histoManager->Add( CexmcAbsorbedEnergy_EDT_Histo, 0,
                       edStore->calorimeterEDLeft,
                       edStore->calorimeterEDRight );

    for ( CexmcAngularRangeList::const_iterator
                        k( triggeredAngularRanges.begin() );
                                    k != triggeredAngularRanges.end(); ++k )
    {
        histoManager->Add( CexmcAbsEnInLeftCalorimeter_ARReal_EDT_Histo,
                           k->index, edStore->calorimeterEDLeft );
        histoManager->Add( CexmcAbsEnInRightCalorimeter_ARReal_EDT_Histo,
                           k->index, edStore->calorimeterEDRight );
    }
}


void  CexmcEventAction::FillTPTHistos( const CexmcTrackPointsStore *  tpStore,
                const CexmcProductionModelData &  pmData,
                const CexmcAngularRangeList &  triggeredAngularRanges ) const
{
    CexmcHistoManager *  histoManager( CexmcHistoManager::Instance() );

    if ( tpStore->monitorTP.IsValid() )
    {
        histoManager->Add( CexmcMomentumBP_TPT_Histo, 0,
                           tpStore->monitorTP.momentumAmp );
        histoManager->Add( CexmcTPInMonitor_TPT_Histo, 0,
                           tpStore->monitorTP.positionLocal.x(),
                           tpStore->monitorTP.positionLocal.y() );
    }

    if ( tpStore->targetTPOutputParticle.IsValid() )
    {
        histoManager->Add( CexmcTPInTarget_TPT_Histo, 0,
                           tpStore->targetTPOutputParticle.positionLocal.x(),
                           tpStore->targetTPOutputParticle.positionLocal.y(),
                           tpStore->targetTPOutputParticle.positionLocal.z() );
        if ( histoManager->GetVerboseLevel() > 0 )
        {
            histoManager->Add( CexmcMomentumIP_TPT_Histo, 0,
                               pmData.incidentParticleLAB.rho() );
        }
    }

    for ( CexmcAngularRangeList::const_iterator
                        k( triggeredAngularRanges.begin() );
                                    k != triggeredAngularRanges.end(); ++k )
    {
        if ( tpStore->calorimeterTPLeft.IsValid() )
        {
            /* kinetic energy and momentum of gamma are equal */
            G4double  kinEnergy( tpStore->calorimeterTPLeft.momentumAmp );
            histoManager->Add( CexmcKinEnAtLeftCalorimeter_ARReal_TPT_Histo,
                               k->index, kinEnergy );
        }
        if ( tpStore->calorimeterTPRight.IsValid() )
        {
            G4double  kinEnergy( tpStore->calorimeterTPRight.momentumAmp );
            histoManager->Add( CexmcKinEnAtRightCalorimeter_ARReal_TPT_Histo,
                               k->index, kinEnergy );
        }
        if ( tpStore->targetTPOutputParticle.IsValid() )
        {
            histoManager->Add( CexmcTPInTarget_ARReal_TPT_Histo, k->index,
                           tpStore->targetTPOutputParticle.positionLocal.x(),
                           tpStore->targetTPOutputParticle.positionLocal.y(),
                           tpStore->targetTPOutputParticle.positionLocal.z() );
            histoManager->Add( CexmcKinEnOP_LAB_ARReal_TPT_Histo, k->index,
                               opKinEnergy );
            histoManager->Add( CexmcAngleOP_SCM_ARReal_TPT_Histo, k->index,
                               pmData.outputParticleSCM.cosTheta() );
        }
        if ( tpStore->targetTPOutputParticleDecayProductParticle1.IsValid() &&
             tpStore->targetTPOutputParticleDecayProductParticle2.IsValid() )
        {
            G4double  openAngle(
                        tpStore->targetTPOutputParticleDecayProductParticle1.
                        directionWorld.angle( tpStore->
                                targetTPOutputParticleDecayProductParticle2.
                                        directionWorld ) / deg );
            histoManager->Add( CexmcOpenAngle_ARReal_TPT_Histo, k->index,
                               openAngle );
        }
    }
}


void  CexmcEventAction::FillRTHistos( G4bool  reconstructorHasFullTrigger,
                const CexmcEnergyDepositStore *  edStore,
                const CexmcTrackPointsStore *  tpStore,
                const CexmcProductionModelData &  pmData,
                const CexmcAngularRangeList &  triggeredAngularRanges ) const
{
    CexmcHistoManager *  histoManager( CexmcHistoManager::Instance() );

    G4double    opMass( reconstructor->GetOutputParticleMass() );
    G4double    nopMass( reconstructor->GetNucleusOutputParticleMass() );

    histoManager->Add( CexmcRecMasses_EDT_Histo, 0, opMass, nopMass );

    for ( CexmcAngularRangeList::const_iterator
                        k( triggeredAngularRanges.begin() );
                                    k != triggeredAngularRanges.end(); ++k )
    {
        if ( tpStore->calorimeterTPLeft.IsValid() )
        {
            histoManager->Add( CexmcOPDPAtLeftCalorimeter_ARReal_EDT_Histo,
                               k->index,
                               tpStore->calorimeterTPLeft.positionLocal.x(),
                               tpStore->calorimeterTPLeft.positionLocal.y() );
        }
        if ( tpStore->calorimeterTPRight.IsValid() )
        {
            histoManager->Add( CexmcOPDPAtRightCalorimeter_ARReal_EDT_Histo,
                               k->index,
                               tpStore->calorimeterTPRight.positionLocal.x(),
                               tpStore->calorimeterTPRight.positionLocal.y() );
        }
        histoManager->Add( CexmcRecOPDPAtLeftCalorimeter_ARReal_EDT_Histo,
                           k->index,
                           reconstructor->GetCalorimeterEPLeftPosition().x(),
                           reconstructor->GetCalorimeterEPLeftPosition().y() );
        histoManager->Add( CexmcRecOPDPAtRightCalorimeter_ARReal_EDT_Histo,
                           k->index,
                           reconstructor->GetCalorimeterEPRightPosition().x(),
                           reconstructor->GetCalorimeterEPRightPosition().y() );
    }

    if ( ! reconstructorHasFullTrigger )
        return;

    if ( tpStore->monitorTP.IsValid() )
    {
        histoManager->Add( CexmcMomentumBP_RT_Histo, 0,
                           tpStore->monitorTP.momentumAmp );
    }

    if ( tpStore->targetTPOutputParticle.IsValid() )
    {
        histoManager->Add( CexmcTPInTarget_RT_Histo, 0,
                           tpStore->targetTPOutputParticle.positionLocal.x(),
                           tpStore->targetTPOutputParticle.positionLocal.y(),
                           tpStore->targetTPOutputParticle.positionLocal.z() );
    }

    histoManager->Add( CexmcRecMasses_RT_Histo, 0,
                       reconstructor->GetOutputParticleMass(),
                       reconstructor->GetNucleusOutputParticleMass() );

    histoManager->Add( CexmcAbsorbedEnergy_RT_Histo, 0,
                       edStore->calorimeterEDLeft,
                       edStore->calorimeterEDRight );

    G4double  recCosTheta( reconstructor->GetProductionModelData().
                                               outputParticleSCM.cosTheta());

    for ( CexmcAngularRangeList::const_iterator
                        k( triggeredAngularRanges.begin() );
                                    k != triggeredAngularRanges.end(); ++k )
    {
        histoManager->Add( CexmcRecMassOP_ARReal_RT_Histo, k->index, opMass );
        histoManager->Add( CexmcRecMassNOP_ARReal_RT_Histo, k->index, nopMass );
        if ( tpStore->calorimeterTPLeft.IsValid() )
        {
            G4double  kinEnergy( tpStore->calorimeterTPLeft.momentumAmp );
            histoManager->Add( CexmcKinEnAtLeftCalorimeter_ARReal_RT_Histo,
                               k->index, kinEnergy );
            histoManager->Add( CexmcMissEnFromLeftCalorimeter_ARReal_RT_Histo,
                               k->index,
                               kinEnergy - edStore->calorimeterEDLeft );
        }
        if ( tpStore->calorimeterTPRight.IsValid() )
        {
            G4double  kinEnergy( tpStore->calorimeterTPRight.momentumAmp );
            histoManager->Add( CexmcKinEnAtRightCalorimeter_ARReal_RT_Histo,
                               k->index, kinEnergy );
            histoManager->Add( CexmcMissEnFromRightCalorimeter_ARReal_RT_Histo,
                               k->index,
                               kinEnergy - edStore->calorimeterEDRight );
        }
        if ( tpStore->targetTPOutputParticle.IsValid() )
        {
            histoManager->Add( CexmcTPInTarget_ARReal_RT_Histo, k->index,
                           tpStore->targetTPOutputParticle.positionLocal.x(),
                           tpStore->targetTPOutputParticle.positionLocal.y(),
                           tpStore->targetTPOutputParticle.positionLocal.z() );
            histoManager->Add( CexmcKinEnOP_LAB_ARReal_RT_Histo, k->index,
                               opKinEnergy );
            histoManager->Add( CexmcAngleOP_SCM_ARReal_RT_Histo, k->index,
                               pmData.outputParticleSCM.cosTheta() );
            G4double  diffCosTheta( pmData.outputParticleSCM.cosTheta() -
                                    recCosTheta );
            histoManager->Add( CexmcDiffAngleOP_SCM_ARReal_RT_Histo, k->index,
                               diffCosTheta );
        }
        if ( tpStore->targetTPOutputParticleDecayProductParticle1.IsValid() &&
             tpStore->targetTPOutputParticleDecayProductParticle2.IsValid() )
        {
            G4double  openAngle(
                        tpStore->targetTPOutputParticleDecayProductParticle1.
                        directionWorld.angle( tpStore->
                                targetTPOutputParticleDecayProductParticle2.
                                        directionWorld ) / deg );
            histoManager->Add( CexmcOpenAngle_ARReal_RT_Histo, k->index,
                               openAngle );
            G4double  diffOpenAngle( openAngle - reconstructor->GetTheAngle() /
                                     deg );
            histoManager->Add( CexmcDiffOpenAngle_ARReal_RT_Histo, k->index,
                               diffOpenAngle );
        }
        if ( tpStore->calorimeterTPLeft.IsValid() )
        {
            histoManager->Add( CexmcOPDPAtLeftCalorimeter_ARReal_RT_Histo,
                               k->index,
                               tpStore->calorimeterTPLeft.positionLocal.x(),
                               tpStore->calorimeterTPLeft.positionLocal.y() );
        }
        if ( tpStore->calorimeterTPRight.IsValid() )
        {
            histoManager->Add( CexmcOPDPAtRightCalorimeter_ARReal_RT_Histo,
                               k->index,
                               tpStore->calorimeterTPRight.positionLocal.x(),
                               tpStore->calorimeterTPRight.positionLocal.y() );
        }
        histoManager->Add( CexmcAbsEnInLeftCalorimeter_ARReal_RT_Histo,
                           k->index, edStore->calorimeterEDLeft );
        histoManager->Add( CexmcAbsEnInRightCalorimeter_ARReal_RT_Histo,
                           k->index, edStore->calorimeterEDRight );
        histoManager->Add( CexmcRecAngleOP_SCM_ARReal_RT_Histo,
                           k->index, recCosTheta );
        histoManager->Add( CexmcRecOpenAngle_ARReal_RT_Histo,
                           k->index, reconstructor->GetTheAngle() / deg );
        histoManager->Add( CexmcRecOPDPAtLeftCalorimeter_ARReal_RT_Histo,
                           k->index,
                           reconstructor->GetCalorimeterEPLeftPosition().x(),
                           reconstructor->GetCalorimeterEPLeftPosition().y() );
        histoManager->Add( CexmcRecOPDPAtRightCalorimeter_ARReal_RT_Histo,
                           k->index,
                           reconstructor->GetCalorimeterEPRightPosition().x(),
                           reconstructor->GetCalorimeterEPRightPosition().y() );
    }
}

#endif


void  CexmcEventAction::DrawTrajectories( const G4Event *  event )
{
    G4VisManager *  visManager( static_cast< G4VisManager * >(
                                    G4VVisManager::GetConcreteInstance() ) );
    if ( ! visManager || ! visManager->GetCurrentGraphicsSystem() )
        return;

    G4int                    nTraj( 0 );
    G4TrajectoryContainer *  trajContainer( event->GetTrajectoryContainer() );

    if ( ! trajContainer )
        return;

    nTraj = trajContainer->entries();

    for ( int  i( 0 ); i < nTraj; ++i )
    {
        G4VTrajectory *  traj( ( *trajContainer )[ i ] );
        traj->DrawTrajectory();
    }
}


void  CexmcEventAction::DrawTrackPoints(
                                const CexmcTrackPointsStore *  tpStore ) const
{
    G4VisManager *  visManager( static_cast< G4VisManager * >(
                                    G4VVisManager::GetConcreteInstance() ) );
    if ( ! visManager || ! visManager->GetCurrentGraphicsSystem() )
        return;

    G4Circle         circle;
    G4VisAttributes  visAttributes( CexmcTrackPointsMarkerColour );
    circle.SetScreenSize( CexmcSmallCircleScreenSize );
    circle.SetFillStyle( G4Circle::filled );
    circle.SetVisAttributes( visAttributes );

    if ( tpStore->monitorTP.IsValid() )
    {
        circle.SetPosition( tpStore->monitorTP.positionWorld );
        visManager->Draw( circle );
    }

    if ( tpStore->targetTPBeamParticle.IsValid() )
    {
        circle.SetPosition( tpStore->targetTPBeamParticle.positionWorld );
        visManager->Draw( circle );
    }

    if ( tpStore->targetTPOutputParticle.IsValid() )
    {
        circle.SetPosition( tpStore->targetTPOutputParticle.positionWorld );
        visManager->Draw( circle );
    }

    if ( tpStore->vetoCounterTPLeft.IsValid() )
    {
        circle.SetPosition( tpStore->vetoCounterTPLeft.positionWorld );
        visManager->Draw( circle );
    }

    if ( tpStore->vetoCounterTPRight.IsValid() )
    {
        circle.SetPosition( tpStore->vetoCounterTPRight.positionWorld );
        visManager->Draw( circle );
    }

    if ( tpStore->calorimeterTPLeft.IsValid() )
    {
        circle.SetPosition( tpStore->calorimeterTPLeft.positionWorld );
        visManager->Draw( circle );
    }

    if ( tpStore->calorimeterTPRight.IsValid() )
    {
        circle.SetPosition( tpStore->calorimeterTPRight.positionWorld );
        visManager->Draw( circle );
    }
}


void  CexmcEventAction::DrawReconstructionData( void )
{
    G4VisManager *  visManager( static_cast< G4VisManager * >(
                                    G4VVisManager::GetConcreteInstance() ) );
    if ( ! visManager || ! visManager->GetCurrentGraphicsSystem() )
        return;

    G4Circle  circle( reconstructor->GetTargetEPWorldPosition() );
    circle.SetScreenSize( CexmcSmallCircleScreenSize );
    circle.SetFillStyle( G4Circle::filled );
    G4VisAttributes  visAttributes( CexmcRecTrackPointsMarkerColour );
    circle.SetVisAttributes( visAttributes );
    visManager->Draw( circle );

    circle.SetScreenSize( CexmcBigCircleScreenSize );
    circle.SetPosition( reconstructor->GetCalorimeterEPLeftWorldPosition() );
    visManager->Draw( circle );

    circle.SetPosition( reconstructor->GetCalorimeterEPRightWorldPosition() );
    visManager->Draw( circle );
}


void  CexmcEventAction::UpdateRunHits(
                                    const CexmcAngularRangeList &  aRangesReal,
                                    const CexmcAngularRangeList &  aRangesRec,
                                    G4bool  tpDigitizerHasTriggered,
                                    G4bool  edDigitizerHasTriggered,
                                    G4bool  edDigitizerMonitorHasTriggered,
                                    G4bool  reconstructorHasFullTrigger,
                                    const CexmcAngularRange &  aGap )
{
    G4RunManager *    runManager( G4RunManager::GetRunManager() );
    const CexmcRun *  run( static_cast< const CexmcRun * >(
                                                runManager->GetCurrentRun() ) );
    CexmcRun *        theRun( const_cast< CexmcRun * >( run ) );

    if ( tpDigitizerHasTriggered )
    {
        for ( CexmcAngularRangeList::const_iterator  k( aRangesReal.begin() );
                                                k != aRangesReal.end(); ++k )
        {
            theRun->IncrementNmbOfHitsSampledFull( k->index );
            if ( edDigitizerMonitorHasTriggered )
                theRun->IncrementNmbOfHitsSampled( k->index );
            if ( reconstructorHasFullTrigger )
                theRun->IncrementNmbOfHitsTriggeredRealRange( k->index );
        }
        if ( reconstructorHasFullTrigger )
        {
            if ( aRangesRec.empty() )
            {
                theRun->IncrementNmbOfOrphanHits( aGap.index );
            }
            else
            {
                for ( CexmcAngularRangeList::const_iterator
                        k( aRangesRec.begin() ); k != aRangesRec.end(); ++k )
                {
                    theRun->IncrementNmbOfHitsTriggeredRecRange( k->index );
                }
            }
        }
    }
    else
    {
        if ( edDigitizerHasTriggered )
            theRun->IncrementNmbOfFalseHitsTriggeredEDT();
        if ( reconstructorHasFullTrigger )
            theRun->IncrementNmbOfFalseHitsTriggeredRec();
    }
}


#ifdef CEXMC_USE_PERSISTENCY

void  CexmcEventAction::SaveEvent( const G4Event *  event,
                                   G4bool  edDigitizerHasTriggered,
                                   const CexmcEnergyDepositStore *  edStore,
                                   const CexmcTrackPointsStore *  tpStore,
                                   const CexmcProductionModelData &  pmData )
{
    CexmcRunManager *  runManager( static_cast< CexmcRunManager * >(
                                            G4RunManager::GetRunManager() ) );
    if ( ! runManager->ProjectIsSaved() )
        return;

    if ( runManager->GetEventDataVerboseLevel() == CexmcWriteNoEventData )
        return;

    if ( ! edDigitizerHasTriggered && runManager->GetEventDataVerboseLevel() !=
                                                CexmcWriteEventDataOnEveryTPT )
        return;

    boost::archive::binary_oarchive *  archive(
                                            runManager->GetEventsArchive() );
    if ( archive )
    {
        CexmcEventSObject  sObject = { event->GetEventID(),
            edDigitizerHasTriggered, edStore->monitorED,
            edStore->vetoCounterEDLeft, edStore->vetoCounterEDRight,
            edStore->calorimeterEDLeft, edStore->calorimeterEDRight,
            edStore->calorimeterEDLeftCollection,
            edStore->calorimeterEDRightCollection,
            tpStore->monitorTP, tpStore->targetTPBeamParticle,
            tpStore->targetTPOutputParticle, tpStore->targetTPNucleusParticle,
            tpStore->targetTPOutputParticleDecayProductParticle1,
            tpStore->targetTPOutputParticleDecayProductParticle2,
            tpStore->vetoCounterTPLeft, tpStore->vetoCounterTPRight,
            tpStore->calorimeterTPLeft, tpStore->calorimeterTPRight, pmData };
        archive->operator<<( sObject );
        const CexmcRun *  run( static_cast< const CexmcRun * >(
                                                runManager->GetCurrentRun() ) );
        CexmcRun *        theRun( const_cast< CexmcRun * >( run ) );
        theRun->IncrementNmbOfSavedEvents();
    }
}


void  CexmcEventAction::SaveEventFast( const G4Event *  event,
                                       G4bool  tpDigitizerHasTriggered,
                                       G4bool  edDigitizerHasTriggered,
                                       G4bool  edDigitizerMonitorHasTriggered,
                                       G4double  opCosThetaSCM )
{
    CexmcRunManager *  runManager( static_cast< CexmcRunManager * >(
                                            G4RunManager::GetRunManager() ) );
    if ( ! runManager->ProjectIsSaved() )
        return;

    if ( runManager->GetEventDataVerboseLevel() == CexmcWriteNoEventData )
        return;

    boost::archive::binary_oarchive *  archive(
                                        runManager->GetFastEventsArchive() );
    if ( archive )
    {
        if ( ! tpDigitizerHasTriggered )
            opCosThetaSCM = CexmcInvalidCosTheta;

        CexmcEventFastSObject  sObject = { event->GetEventID(), opCosThetaSCM,
                                           edDigitizerHasTriggered,
                                           edDigitizerMonitorHasTriggered };
        archive->operator<<( sObject );
        const CexmcRun *  run( static_cast< const CexmcRun * >(
                                                runManager->GetCurrentRun() ) );
        CexmcRun *        theRun( const_cast< CexmcRun * >( run ) );
        theRun->IncrementNmbOfSavedFastEvents();
    }
}

#endif


void  CexmcEventAction::EndOfEventAction( const G4Event *  event )
{
    G4DigiManager *                digiManager( G4DigiManager::GetDMpointer() );
    CexmcEnergyDepositDigitizer *  energyDepositDigitizer(
            static_cast< CexmcEnergyDepositDigitizer * >( digiManager->
                                FindDigitizerModule( CexmcEDDigitizerName ) ) );
    CexmcTrackPointsDigitizer *    trackPointsDigitizer(
            static_cast< CexmcTrackPointsDigitizer * >( digiManager->
                                FindDigitizerModule( CexmcTPDigitizerName ) ) );

    energyDepositDigitizer->Digitize();
    trackPointsDigitizer->Digitize();

    G4bool  edDigitizerMonitorHasTriggered(
                                energyDepositDigitizer->MonitorHasTriggered() );
    G4bool  edDigitizerHasTriggered( false );

    CexmcEventInfo *  eventInfo( static_cast< CexmcEventInfo * >(
                                                event->GetUserInformation() ) );
    if ( ! eventInfo || eventInfo->EdTriggerIsOk() )
        edDigitizerHasTriggered = energyDepositDigitizer->HasTriggered();

    G4bool  tpDigitizerHasTriggered( trackPointsDigitizer->HasTriggered() );
    G4bool  reconstructorHasBasicTrigger( false );
    G4bool  reconstructorHasFullTrigger( false );

    CexmcEnergyDepositStore *  edStore( MakeEnergyDepositStore(
                                                    energyDepositDigitizer ) );
    CexmcTrackPointsStore *    tpStore( MakeTrackPointsStore(
                                                    trackPointsDigitizer ) );

    try
    {
        CexmcProductionModel *  productionModel(
                                        physicsManager->GetProductionModel() );

        if ( ! productionModel )
            throw CexmcException( CexmcWeirdException );

        const CexmcAngularRangeList &     angularRanges(
                                productionModel->GetAngularRanges() );
        const CexmcAngularRangeList &     triggeredAngularRanges(
                                productionModel->GetTriggeredAngularRanges() );
        const CexmcProductionModelData &  pmData(
                                productionModel->GetProductionModelData() );

        if ( edDigitizerHasTriggered )
        {
            reconstructor->Reconstruct( edStore );
            reconstructorHasBasicTrigger = reconstructor->HasBasicTrigger();
            reconstructorHasFullTrigger = reconstructor->HasFullTrigger();
        }

        CexmcAngularRangeList  triggeredRecAngularRanges;

        if ( reconstructorHasBasicTrigger )
        {
            for ( CexmcAngularRangeList::const_iterator
                  k( angularRanges.begin() ); k != angularRanges.end(); ++k )
            {
                G4double  cosTheta( reconstructor->GetProductionModelData().
                                    outputParticleSCM.cosTheta() );
                if ( cosTheta <= k->top && cosTheta > k->bottom )
                    triggeredRecAngularRanges.push_back( CexmcAngularRange(
                                                k->top, k->bottom, k->index ) );
            }
        }

        CexmcAngularRange  angularGap( 0.0, 0.0, 0 );
        if ( triggeredRecAngularRanges.empty() )
        {
            CexmcAngularRangeList  angularGaps;
            GetAngularGaps( angularRanges, angularGaps );
            for ( CexmcAngularRangeList::const_iterator
                    k( angularGaps.begin() ); k != angularGaps.end(); ++k )
            {
                G4double  cosTheta( reconstructor->GetProductionModelData().
                                    outputParticleSCM.cosTheta() );
                if ( cosTheta <= k->top && cosTheta > k->bottom )
                {
                    angularGap = *k;
                    break;
                }
            }
        }

        UpdateRunHits( triggeredAngularRanges, triggeredRecAngularRanges,
                       tpDigitizerHasTriggered, edDigitizerHasTriggered,
                       edDigitizerMonitorHasTriggered,
                       reconstructorHasFullTrigger, angularGap );

        if ( verbose > 0 )
        {
            G4bool  printMessages( verbose > 3 ||
                        ( ( verbose == 1 ) && tpDigitizerHasTriggered ) ||
                        ( ( verbose == 2 ) && edDigitizerHasTriggered ) ||
                        ( ( verbose == 3 ) && ( tpDigitizerHasTriggered ||
                                                edDigitizerHasTriggered ) ) );
            if ( printMessages )
            {
                G4cout << "Event " << event->GetEventID() << G4endl;
                if ( tpDigitizerHasTriggered )
                {
                    PrintTrackPoints( tpStore );
                    PrintProductionModelData( triggeredAngularRanges, pmData );
                }
                if ( reconstructorHasBasicTrigger )
                    PrintReconstructedData( triggeredRecAngularRanges,
                                            angularGap );
                if ( edDigitizerHasTriggered )
                    PrintEnergyDeposit( edStore );
            }
        }

        if ( verboseDraw > 0 )
        {
            G4bool  drawTrajectories( verboseDraw > 3 ||
                        ( ( verboseDraw == 1 ) && tpDigitizerHasTriggered ) ||
                        ( ( verboseDraw == 2 ) && edDigitizerHasTriggered ) ||
                        ( ( verboseDraw == 3 ) && ( tpDigitizerHasTriggered ||
                                                edDigitizerHasTriggered ) ) );
            if ( drawTrajectories )
            {
                DrawTrajectories( event );
                if ( tpDigitizerHasTriggered )
                    DrawTrackPoints( tpStore );
                if ( reconstructorHasBasicTrigger )
                    DrawReconstructionData();
            }
        }

#ifdef CEXMC_USE_PERSISTENCY
        if ( edDigitizerHasTriggered || tpDigitizerHasTriggered )
        {
            SaveEventFast( event, tpDigitizerHasTriggered,
                           edDigitizerHasTriggered,
                           edDigitizerMonitorHasTriggered,
                           pmData.outputParticleSCM.cosTheta() );
            SaveEvent( event, edDigitizerHasTriggered, edStore, tpStore,
                       pmData );
        }
#endif

#ifdef CEXMC_USE_ROOT
        /* opKinEnergy will be used in several histos */
        if ( tpStore->targetTPOutputParticle.IsValid() )
        {
            opKinEnergy = CexmcGetKinEnergy(
                    tpStore->targetTPOutputParticle.momentumAmp,
                    tpStore->targetTPOutputParticle.particle->GetPDGMass() );
        }

        if ( edDigitizerHasTriggered )
            FillEDTHistos( edStore, triggeredAngularRanges );

        /* fill TPT histos only when the monitor has triggered because events
         * when it was missed have less value for us */
        if ( tpDigitizerHasTriggered && edDigitizerMonitorHasTriggered )
            FillTPTHistos( tpStore, pmData, triggeredAngularRanges );

        if ( reconstructorHasBasicTrigger )
            FillRTHistos( reconstructorHasFullTrigger, edStore, tpStore,
                          pmData, triggeredAngularRanges );
#endif

        G4Event *  theEvent( const_cast< G4Event * >( event ) );
        if ( eventInfo )
        {
            delete eventInfo;
            theEvent->SetUserInformation( NULL );
        }
        theEvent->SetUserInformation( new CexmcEventInfo(
                                                edDigitizerHasTriggered,
                                                tpDigitizerHasTriggered,
                                                reconstructorHasFullTrigger ) );
    }
    catch ( CexmcException &  e )
    {
        G4cout << e.what() << G4endl;
    }
    catch ( ... )
    {
        G4cout << "Unknown exception caught" << G4endl;
    }

    delete edStore;
    delete tpStore;
}

