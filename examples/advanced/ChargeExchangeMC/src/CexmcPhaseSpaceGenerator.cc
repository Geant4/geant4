/*
 * =============================================================================
 *
 *       Filename:  CexmcPhaseSpaceGenerator.cc
 *
 *    Description:  phase space generator interface
 *
 *        Version:  1.0
 *        Created:  08.09.2010 14:05:51
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#include "CexmcPhaseSpaceGenerator.hh"


CexmcPhaseSpaceGenerator::CexmcPhaseSpaceGenerator() :
    fermiEnergyDepIsOn( false ), totalEnergy( 0. ), totalMass( 0. )
{
}


CexmcPhaseSpaceGenerator::~CexmcPhaseSpaceGenerator()
{
}


void  CexmcPhaseSpaceGenerator::SetParticles(
                                    const CexmcPhaseSpaceInVector &  inVec_,
                                    const CexmcPhaseSpaceOutVector &  outVec_ )
{
    inVec = inVec_;
    outVec = outVec_;
    for ( CexmcPhaseSpaceOutVector::const_iterator  k( outVec.begin() );
                                                        k != outVec.end(); ++k )
    {
        totalMass += k->mass;
    }

    ParticleChangeHook();
}


void  CexmcPhaseSpaceGenerator::SetFermiEnergyDependence( G4bool  on )
{
    fermiEnergyDepIsOn = on;

    FermiEnergyDepStatusChangeHook();
}


G4bool  CexmcPhaseSpaceGenerator::CheckKinematics( void )
{
    totalEnergy = 0.;
    for ( CexmcPhaseSpaceInVector::const_iterator  k( inVec.begin() );
                                                        k != inVec.end(); ++k )
    {
        totalEnergy += ( *k )->e();
    }

    return totalEnergy > totalMass;
}


void  CexmcPhaseSpaceGenerator::ParticleChangeHook( void )
{
}


void  CexmcPhaseSpaceGenerator::FermiEnergyDepStatusChangeHook( void )
{
}

