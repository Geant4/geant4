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

#ifdef CEXMC_USE_PERSISTENCY

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

#endif

