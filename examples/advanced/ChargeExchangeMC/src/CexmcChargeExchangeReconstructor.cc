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
 *       Filename:  CexmcChargeExchangeReconstructor.cc
 *
 *    Description:  charge exchange reconstructor
 *
 *        Version:  1.0
 *        Created:  02.12.2009 15:17:13
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <cmath>
#include <G4ThreeVector.hh>
#include <G4LorentzVector.hh>
#include "CexmcChargeExchangeReconstructor.hh"
#include "CexmcChargeExchangeReconstructorMessenger.hh"
#include "CexmcEnergyDepositStore.hh"
#include "CexmcPrimaryGeneratorAction.hh"
#include "CexmcParticleGun.hh"
#include "CexmcProductionModel.hh"
#include "CexmcRunManager.hh"
#include "CexmcException.hh"


CexmcChargeExchangeReconstructor::CexmcChargeExchangeReconstructor(
                            const CexmcProductionModel *  productionModel ) :
    outputParticleMass( 0 ), nucleusOutputParticleMass( 0 ),
    useTableMass( false ), useMassCut( false ), massCutOPCenter( 0 ),
    massCutNOPCenter( 0 ), massCutOPWidth( 0 ), massCutNOPWidth( 0 ),
    massCutEllipseAngle( 0 ), useAbsorbedEnergyCut( false ),
    absorbedEnergyCutCLCenter( 0 ), absorbedEnergyCutCRCenter( 0 ),
    absorbedEnergyCutCLWidth( 0 ), absorbedEnergyCutCRWidth( 0 ),
    absorbedEnergyCutEllipseAngle( 0 ), expectedMomentumAmp( -1 ),
    edCollectionAlgorithm( CexmcCollectEDInAllCrystals ),
    hasMassCutTriggered( false ), hasAbsorbedEnergyCutTriggered( false ),
    beamParticleIsInitialized( false ), particleGun( NULL ), messenger( NULL )
{
    if ( ! productionModel )
        throw CexmcException( CexmcWeirdException );

    productionModelData.incidentParticle =
                                    productionModel->GetIncidentParticle();

    CexmcRunManager *  runManager( static_cast< CexmcRunManager * >(
                                            G4RunManager::GetRunManager() ) );
    const CexmcPrimaryGeneratorAction *  primaryGeneratorAction(
                        static_cast< const CexmcPrimaryGeneratorAction * >(
                                runManager->GetUserPrimaryGeneratorAction() ) );
    CexmcPrimaryGeneratorAction *  thePrimaryGeneratorAction(
                        const_cast< CexmcPrimaryGeneratorAction * >(
                                primaryGeneratorAction ) );
    particleGun = thePrimaryGeneratorAction->GetParticleGun();

    productionModelData.nucleusParticle =
                                    productionModel->GetNucleusParticle();
    productionModelData.outputParticle =
                                    productionModel->GetOutputParticle();
    productionModelData.nucleusOutputParticle =
                                    productionModel->GetNucleusOutputParticle();

    messenger = new CexmcChargeExchangeReconstructorMessenger( this );
}


CexmcChargeExchangeReconstructor::~CexmcChargeExchangeReconstructor()
{
    delete messenger;
}


void  CexmcChargeExchangeReconstructor::SetupBeamParticle( void )
{
    if ( *productionModelData.incidentParticle !=
                                        *particleGun->GetParticleDefinition() )
        throw CexmcException( CexmcBeamAndIncidentParticlesMismatch );

    beamParticleIsInitialized = true;
}


void  CexmcChargeExchangeReconstructor::Reconstruct(
                                    const CexmcEnergyDepositStore *  edStore )
{
    if ( ! beamParticleIsInitialized )
    {
        if ( *productionModelData.incidentParticle !=
                                        *particleGun->GetParticleDefinition() )
            throw CexmcException( CexmcBeamAndIncidentParticlesMismatch );

        beamParticleIsInitialized = true;
    }

    if ( edCollectionAlgorithm == CexmcCollectEDInAdjacentCrystals )
        collectEDInAdjacentCrystals = true;

    ReconstructEntryPoints( edStore );
    if ( hasBasicTrigger )
        ReconstructTargetPoint();
    if ( hasBasicTrigger )
        ReconstructAngle();

    G4ThreeVector  epLeft( calorimeterEPLeftWorldPosition -
                           targetEPWorldPosition );
    G4ThreeVector  epRight( calorimeterEPRightWorldPosition -
                            targetEPWorldPosition );

    G4double  cosTheAngle( std::cos( theAngle ) );
    G4double  calorimeterEDLeft( edStore->calorimeterEDLeft );
    G4double  calorimeterEDRight( edStore->calorimeterEDRight );

    if ( edCollectionAlgorithm == CexmcCollectEDInAdjacentCrystals )
    {
        calorimeterEDLeft = calorimeterEDLeftAdjacent;
        calorimeterEDRight = calorimeterEDRightAdjacent;
    }

    //G4double  cosOutputParticleLAB(
        //( calorimeterEDLeft * cosAngleLeft +
          //calorimeterEDRight * cosAngleRight ) /
          //std::sqrt( calorimeterEDLeft * calorimeterEDLeft +
                     //calorimeterEDRight * calorimeterEDRight +
                     //calorimeterEDLeft * calorimeterEDRight * cosTheAngle ) );

    outputParticleMass = std::sqrt( 2 * calorimeterEDLeft *
                                    calorimeterEDRight * ( 1 - cosTheAngle ) );

    G4ThreeVector    opdpLeftMomentum( epLeft );
    opdpLeftMomentum.setMag( calorimeterEDLeft );
    G4ThreeVector    opdpRightMomentum( epRight );
    opdpRightMomentum.setMag( calorimeterEDRight );
    G4ThreeVector    opMomentum( opdpLeftMomentum + opdpRightMomentum );

    /* opMass will be used only in calculation of output particle's total
     * energy, in other places outputParticleMass should be used instead */
    G4double         opMass( useTableMass ?
                             productionModelData.outputParticle->GetPDGMass() :
                             outputParticleMass );
    /* the formula below is equivalent to
     * calorimeterEDLeft + calorimeterEDRight if opMass = outputParticleMass */
    G4double         opEnergy( std::sqrt( opMomentum.mag2() +
                                          opMass * opMass ) );
    productionModelData.outputParticleLAB = G4LorentzVector( opMomentum,
                                                             opEnergy );

    G4ThreeVector  incidentParticleMomentum( particleGun->GetOrigDirection() );
    G4double       incidentParticleMomentumAmp( expectedMomentumAmp > 0 ?
                    expectedMomentumAmp : particleGun->GetOrigMomentumAmp() );
    incidentParticleMomentum *= incidentParticleMomentumAmp;

    G4double       incidentParticlePDGMass(
                        productionModelData.incidentParticle->GetPDGMass() );
    G4double       incidentParticlePDGMass2( incidentParticlePDGMass *
                                             incidentParticlePDGMass );
    G4double       incidentParticleEnergy(
        std::sqrt( incidentParticleMomentumAmp * incidentParticleMomentumAmp +
                   incidentParticlePDGMass2 ) );

    productionModelData.incidentParticleLAB = G4LorentzVector(
                        incidentParticleMomentum, incidentParticleEnergy );
    G4double       nucleusParticlePDGMass(
                        productionModelData.nucleusParticle->GetPDGMass() );
    productionModelData.nucleusParticleLAB = G4LorentzVector(
                        G4ThreeVector( 0, 0, 0 ), nucleusParticlePDGMass );

    G4LorentzVector  lVecSum( productionModelData.incidentParticleLAB +
                        productionModelData.nucleusParticleLAB );
    G4ThreeVector    boostVec( lVecSum.boostVector() );

    productionModelData.nucleusOutputParticleLAB =
            lVecSum - productionModelData.outputParticleLAB;

    productionModelData.incidentParticleSCM =
            productionModelData.incidentParticleLAB;
    productionModelData.nucleusParticleSCM =
            productionModelData.nucleusParticleLAB;
    productionModelData.outputParticleSCM =
            productionModelData.outputParticleLAB;
    productionModelData.nucleusOutputParticleSCM =
            productionModelData.nucleusOutputParticleLAB;

    productionModelData.incidentParticleSCM.boost( -boostVec );
    productionModelData.nucleusParticleSCM.boost( -boostVec );
    productionModelData.outputParticleSCM.boost( -boostVec );
    productionModelData.nucleusOutputParticleSCM.boost( -boostVec );

    G4ThreeVector  nopMomentum( incidentParticleMomentum - opMomentum );
    G4double       nopEnergy( incidentParticleEnergy + nucleusParticlePDGMass -
                              opEnergy );
    nucleusOutputParticleMass = std::sqrt( nopEnergy * nopEnergy -
                                           nopMomentum.mag2() );

    if ( useMassCut )
    {
        G4double  cosMassCutEllipseAngle( std::cos( massCutEllipseAngle ) );
        G4double  sinMassCutEllipseAngle( std::sin( massCutEllipseAngle ) );

        if ( massCutOPWidth <= 0. || massCutNOPWidth <= 0. )
        {
            hasMassCutTriggered = false;
        }
        else
        {
            G4double  massCutOPWidth2( massCutOPWidth * massCutOPWidth );
            G4double  massCutNOPWidth2( massCutNOPWidth * massCutNOPWidth );

            hasMassCutTriggered =
                std::pow( ( outputParticleMass - massCutOPCenter ) *
                              cosMassCutEllipseAngle +
                          ( nucleusOutputParticleMass - massCutNOPCenter ) *
                              sinMassCutEllipseAngle, 2 ) / massCutOPWidth2 +
                std::pow( - ( outputParticleMass - massCutOPCenter ) *
                              sinMassCutEllipseAngle +
                          ( nucleusOutputParticleMass - massCutNOPCenter ) *
                              cosMassCutEllipseAngle, 2 ) / massCutNOPWidth2 <
                1;
        }
    }

    if ( useAbsorbedEnergyCut )
    {
        G4double  cosAbsorbedEnergyCutEllipseAngle(
                                std::cos( absorbedEnergyCutEllipseAngle ) );
        G4double  sinAbsorbedEnergyCutEllipseAngle(
                                std::sin( absorbedEnergyCutEllipseAngle ) );

        if ( absorbedEnergyCutCLWidth <= 0. || absorbedEnergyCutCRWidth <= 0. )
        {
            hasAbsorbedEnergyCutTriggered = false;
        }
        else
        {
            G4double  absorbedEnergyCutCLWidth2(
                        absorbedEnergyCutCLWidth * absorbedEnergyCutCLWidth );
            G4double  absorbedEnergyCutCRWidth2(
                        absorbedEnergyCutCRWidth * absorbedEnergyCutCRWidth );

            hasAbsorbedEnergyCutTriggered =
                std::pow( ( calorimeterEDLeft - absorbedEnergyCutCLCenter ) *
                              cosAbsorbedEnergyCutEllipseAngle +
                          ( calorimeterEDRight - absorbedEnergyCutCRCenter ) *
                              sinAbsorbedEnergyCutEllipseAngle, 2 ) /
                absorbedEnergyCutCLWidth2 +
                std::pow( - ( calorimeterEDLeft - absorbedEnergyCutCLCenter ) *
                              sinAbsorbedEnergyCutEllipseAngle +
                          ( calorimeterEDRight - absorbedEnergyCutCRCenter ) *
                              cosAbsorbedEnergyCutEllipseAngle, 2 ) /
                absorbedEnergyCutCRWidth2 <
                1;
        }
    }

    hasBasicTrigger = true;
}


G4bool  CexmcChargeExchangeReconstructor::HasFullTrigger( void ) const
{
    if ( ! hasBasicTrigger )
        return false;
    if ( useMassCut && ! hasMassCutTriggered )
        return false;
    if ( useAbsorbedEnergyCut && ! hasAbsorbedEnergyCutTriggered )
        return false;

    return true;
}


void  CexmcChargeExchangeReconstructor::SetExpectedMomentumAmpDiff(
                                                            G4double  value )
{
    expectedMomentumAmp = particleGun->GetOrigMomentumAmp() + value;
}

