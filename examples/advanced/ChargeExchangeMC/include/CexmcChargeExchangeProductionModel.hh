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
 *       Filename:  CexmcChargeExchangeProductionModel.hh
 *
 *    Description:  charge exchange physics itself
 *
 *        Version:  1.0
 *        Created:  01.11.2009 00:30:46
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_CHARGE_EXCHANGE_PRODUCTION_MODEL_HH
#define CEXMC_CHARGE_EXCHANGE_PRODUCTION_MODEL_HH

#include <G4HadronicInteraction.hh>
#include <G4HadFinalState.hh>
#include <G4HadProjectile.hh>
#include <G4Nucleus.hh>
#include <G4Proton.hh>
#include <G4Neutron.hh>
#include "CexmcProductionModel.hh"
#include "CexmcGenbod.hh"
#include "CexmcReimplementedGenbod.hh"
#include "CexmcException.hh"


template  < typename  OutputParticle >
class  CexmcChargeExchangeProductionModel : public G4HadronicInteraction,
                                            public CexmcProductionModel
{
    public:
        CexmcChargeExchangeProductionModel();

        ~CexmcChargeExchangeProductionModel();

    public:
        G4HadFinalState *  ApplyYourself( const G4HadProjectile &  projectile,
                                          G4Nucleus &  targetNucleus );

    private:
        G4double                    nucleusParticleMass;

        CexmcPhaseSpaceGenerator *  phaseSpaceGenerator;
};


template  < typename  OutputParticle >
CexmcChargeExchangeProductionModel< OutputParticle >::
                                        CexmcChargeExchangeProductionModel() :
    G4HadronicInteraction( CexmcChargeExchangeInteractionName ),
    CexmcProductionModel( CexmcChargeExchangeProductionModelName ),
    nucleusParticleMass( 0 ), phaseSpaceGenerator( NULL )
{
    incidentParticle = G4PionMinus::Definition();
    nucleusParticle = G4Proton::Definition();
    outputParticle = OutputParticle::Definition();
    nucleusOutputParticle = G4Neutron::Definition();

    nucleusParticleMass = nucleusParticle->GetPDGMass();

    productionModelData.incidentParticle = incidentParticle;
    productionModelData.nucleusParticle = nucleusParticle;
    productionModelData.outputParticle = outputParticle;
    productionModelData.nucleusOutputParticle = nucleusOutputParticle;

    CexmcPhaseSpaceInVector   inVec;

    inVec.push_back( &productionModelData.incidentParticleSCM );
    inVec.push_back( &productionModelData.nucleusParticleSCM );

    CexmcPhaseSpaceOutVector  outVec;

    outVec.push_back( CexmcPhaseSpaceOutVectorElement(
                            &productionModelData.outputParticleSCM,
                            outputParticle->GetPDGMass() ) );
    outVec.push_back( CexmcPhaseSpaceOutVectorElement(
                            &productionModelData.nucleusOutputParticleSCM,
                            nucleusOutputParticle->GetPDGMass() ) );

#ifdef CEXMC_USE_GENBOD
    phaseSpaceGenerator = new CexmcGenbod;
#else
    phaseSpaceGenerator = new CexmcReimplementedGenbod;
#endif

    phaseSpaceGenerator->SetParticles( inVec, outVec );
}


template  < typename  OutputParticle >
CexmcChargeExchangeProductionModel< OutputParticle >::
                                        ~CexmcChargeExchangeProductionModel()
{
    delete phaseSpaceGenerator;
}


template  < typename  OutputParticle >
G4HadFinalState *  CexmcChargeExchangeProductionModel< OutputParticle >::
                            ApplyYourself( const G4HadProjectile &  projectile,
                                           G4Nucleus &  targetNucleus )
{
    theParticleChange.Clear();

    G4double           kinEnergy( projectile.GetKineticEnergy() );
    G4HadProjectile &  theProjectile( const_cast< G4HadProjectile & >(
                                                                projectile ) );
    const G4LorentzRotation &  projToLab(
                                    const_cast< const G4LorentzRotation & >(
                                            theProjectile.GetTrafoToLab() ) );
    productionModelData.incidentParticleLAB = projectile.Get4Momentum();
    productionModelData.incidentParticleLAB.transform( projToLab );
    productionModelData.nucleusParticleLAB.setPx( 0 );
    productionModelData.nucleusParticleLAB.setPy( 0 );
    productionModelData.nucleusParticleLAB.setPz( 0 );
    productionModelData.nucleusParticleLAB.setE( nucleusParticleMass );

    if ( fermiMotionIsOn )
    {
        G4ThreeVector  targetNucleusMomentum(
                                        targetNucleus.GetFermiMomentum() );
        G4double       targetNucleusEnergy(
                            std::sqrt( targetNucleusMomentum.mag2() +
                                nucleusParticleMass * nucleusParticleMass ) );
        productionModelData.nucleusParticleLAB = G4LorentzVector(
                                targetNucleusMomentum, targetNucleusEnergy );
    }
    productionModelData.nucleusParticleLAB.transform( projToLab );
    G4LorentzVector  lVecSum( productionModelData.incidentParticleLAB +
                              productionModelData.nucleusParticleLAB );
    G4ThreeVector    boostVec( lVecSum.boostVector() );

    productionModelData.incidentParticleSCM =
            productionModelData.incidentParticleLAB;
    productionModelData.nucleusParticleSCM =
            productionModelData.nucleusParticleLAB;

    productionModelData.incidentParticleSCM.boost( -boostVec );
    productionModelData.nucleusParticleSCM.boost( -boostVec );

    triggeredAngularRanges.clear();

    if ( ! phaseSpaceGenerator->CheckKinematics() )
    {
        theParticleChange.SetEnergyChange( kinEnergy );
        theParticleChange.SetMomentumChange(
                                    projectile.Get4Momentum().vect().unit());
        return &theParticleChange;
    }

    do
    {
        phaseSpaceGenerator->Generate();
        G4double  cosTheta( productionModelData.outputParticleSCM.cosTheta() );
        for ( CexmcAngularRangeList::iterator  k( angularRanges.begin() );
                                                k != angularRanges.end(); ++k )
        {
            if ( cosTheta <= k->top && cosTheta > k->bottom )
                triggeredAngularRanges.push_back( CexmcAngularRange(
                                                k->top, k->bottom, k->index ) );
        }
    } while ( triggeredAngularRanges.empty() );

    productionModelData.outputParticleLAB =
            productionModelData.outputParticleSCM;
    productionModelData.nucleusOutputParticleLAB =
            productionModelData.nucleusOutputParticleSCM;

    productionModelData.outputParticleLAB.boost( boostVec );
    productionModelData.nucleusOutputParticleLAB.boost( boostVec );

    theParticleChange.SetStatusChange( stopAndKill );
    theParticleChange.SetEnergyChange( 0.0 );

    G4DynamicParticle *  secOutParticle( new G4DynamicParticle(
                            outputParticle,
                            productionModelData.outputParticleLAB ) );
    theParticleChange.AddSecondary( secOutParticle );
    G4DynamicParticle *  secNeutron( new G4DynamicParticle(
                            nucleusOutputParticle,
                            productionModelData.nucleusOutputParticleLAB ) );
    theParticleChange.AddSecondary( secNeutron );

    /* projectile->GetDefinition() shall always be identical to incidentParticle
     * as far as CexmcHadronicProcess::IsApplicable() will check that only
     * incidentParticle is allowed. Here is mostly unnecessary assignment. */
    productionModelData.incidentParticle = projectile.GetDefinition();

    return &theParticleChange;
}


#endif

