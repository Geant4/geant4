/*
 * ============================================================================
 *
 *       Filename:  CexmcSimpleProductionModelDataStore.cc
 *
 *    Description:  serialization helper for production model data
 *
 *        Version:  1.0
 *        Created:  02.01.2010 14:56:53
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <G4ParticleTable.hh>
#include "CexmcSimpleProductionModelDataStore.hh"
#include "CexmcProductionModelData.hh"
#include "CexmcException.hh"


CexmcSimpleProductionModelDataStore::CexmcSimpleProductionModelDataStore()
{
}


CexmcSimpleProductionModelDataStore::CexmcSimpleProductionModelDataStore(
                                    const CexmcProductionModelData &  pmData )
{
    incidentParticleSCM = pmData.incidentParticleSCM;
    incidentParticleLAB = pmData.incidentParticleLAB;
    nucleusParticleSCM = pmData.nucleusParticleSCM;
    nucleusParticleLAB = pmData.nucleusParticleLAB;
    outputParticleSCM = pmData.outputParticleSCM;
    outputParticleLAB = pmData.outputParticleLAB;
    nucleusOutputParticleSCM = pmData.nucleusOutputParticleSCM;
    nucleusOutputParticleLAB = pmData.nucleusOutputParticleLAB;
    incidentParticle = pmData.incidentParticle->GetPDGEncoding();
    nucleusParticle = pmData.nucleusParticle->GetPDGEncoding();
    outputParticle = pmData.outputParticle->GetPDGEncoding();
    nucleusOutputParticle = pmData.nucleusOutputParticle->GetPDGEncoding();
}


CexmcSimpleProductionModelDataStore::operator CexmcProductionModelData() const
{
    G4ParticleDefinition *  incidentParticleDefinition(
                    G4ParticleTable::GetParticleTable()->FindParticle(
                                                    incidentParticle ) );
    if ( ! incidentParticleDefinition )
        throw CexmcException( CexmcWeirdException );
    G4ParticleDefinition *  nucleusParticleDefinition(
                    G4ParticleTable::GetParticleTable()->FindParticle(
                                                    nucleusParticle ) );
    if ( ! nucleusParticleDefinition )
        throw CexmcException( CexmcWeirdException );
    G4ParticleDefinition *  outputParticleDefinition(
                    G4ParticleTable::GetParticleTable()->FindParticle(
                                                    outputParticle ) );
    if ( ! outputParticleDefinition )
        throw CexmcException( CexmcWeirdException );
    G4ParticleDefinition *  nucleusOutputParticleDefinition(
                    G4ParticleTable::GetParticleTable()->FindParticle(
                                                    nucleusOutputParticle ) );
    if ( ! nucleusOutputParticleDefinition )
        throw CexmcException( CexmcWeirdException );

    return CexmcProductionModelData( incidentParticleSCM, incidentParticleLAB,
                    nucleusParticleSCM, nucleusParticleLAB,
                    outputParticleSCM, outputParticleLAB,
                    nucleusOutputParticleSCM, nucleusOutputParticleLAB,
                    incidentParticleDefinition, nucleusParticleDefinition,
                    outputParticleDefinition, nucleusOutputParticleDefinition );
}

