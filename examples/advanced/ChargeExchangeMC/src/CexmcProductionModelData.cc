/*
 * ============================================================================
 *
 *       Filename:  CexmcProductionModelData.cc
 *
 *    Description:  SCM/LAB lorentz vector of the particles in reaction
 *
 *        Version:  1.0
 *        Created:  28.12.2009 22:59:10
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <iostream>
#include "CexmcProductionModelData.hh"


CexmcProductionModelData::CexmcProductionModelData()
{
}


CexmcProductionModelData::CexmcProductionModelData(
                      const G4LorentzVector &  incidentParticleSCM,
                      const G4LorentzVector &  incidentParticleLAB,
                      const G4LorentzVector &  nucleusParticleSCM,
                      const G4LorentzVector &  nucleusParticleLAB,
                      const G4LorentzVector &  outputParticleSCM,
                      const G4LorentzVector &  outputParticleLAB,
                      const G4LorentzVector &  nucleusOutputParticleSCM,
                      const G4LorentzVector &  nucleusOutputParticleLAB,
                      const G4ParticleDefinition *  incidentParticle,
                      const G4ParticleDefinition *  nucleusParticle,
                      const G4ParticleDefinition *  outputParticle,
                      const G4ParticleDefinition *  nucleusOutputParticle ) :
    incidentParticleSCM( incidentParticleSCM ),
    incidentParticleLAB( incidentParticleLAB ),
    nucleusParticleSCM( nucleusParticleSCM ),
    nucleusParticleLAB( nucleusParticleLAB ),
    outputParticleSCM( outputParticleSCM ),
    outputParticleLAB( outputParticleLAB ),
    nucleusOutputParticleSCM( nucleusOutputParticleSCM ),
    nucleusOutputParticleLAB( nucleusOutputParticleLAB ),
    incidentParticle( incidentParticle ), nucleusParticle( nucleusParticle ),
    outputParticle( outputParticle ),
    nucleusOutputParticle( nucleusOutputParticle )
{
}


std::ostream &  operator<<( std::ostream &  out,
                            const CexmcProductionModelData &  data )
{
    std::ostream::fmtflags  savedFlags( out.flags() );
    std::streamsize         prec( out.precision() );

    out.precision( 4 );
    out.flags( std::ios::fixed );

    out << std::endl;
    out << "       Incident particle       (LAB) : " <<
           data.incidentParticle->GetParticleName() << " " <<
           G4BestUnit( data.incidentParticleLAB, "Energy" ) << " -- " <<
           G4BestUnit( data.incidentParticleLAB.e(), "Energy" ) << std::endl;
    out << "                               (SCM) : " <<
           data.incidentParticle->GetParticleName() << " " <<
           G4BestUnit( data.incidentParticleSCM, "Energy" ) << " -- " <<
           G4BestUnit( data.incidentParticleSCM.e(), "Energy" ) << std::endl;
    out << "       Nucleus particle        (LAB) : " <<
           data.nucleusParticle->GetParticleName() << " " <<
           G4BestUnit( data.nucleusParticleLAB, "Energy" ) << " -- " <<
           G4BestUnit( data.nucleusParticleLAB.e(), "Energy" ) << std::endl;
    out << "                               (SCM) : " <<
           data.nucleusParticle->GetParticleName() << " " <<
           G4BestUnit( data.nucleusParticleSCM, "Energy" ) << " -- " <<
           G4BestUnit( data.nucleusParticleSCM.e(), "Energy" ) << std::endl;
    out << "       Output particle         (LAB) : " <<
           data.outputParticle->GetParticleName() << " " <<
           G4BestUnit( data.outputParticleLAB, "Energy" ) << " -- " <<
           G4BestUnit( data.outputParticleLAB.e(), "Energy" ) << std::endl;
    out << "                               (SCM) : " <<
           data.outputParticle->GetParticleName() << " " <<
           G4BestUnit( data.outputParticleSCM, "Energy" ) << " -- " <<
           G4BestUnit( data.outputParticleSCM.e(), "Energy" ) << std::endl;
    out << "       Nucleus output particle (LAB) : " <<
           data.nucleusOutputParticle->GetParticleName() << " " <<
           G4BestUnit( data.nucleusOutputParticleLAB, "Energy" ) << " -- " <<
           G4BestUnit( data.nucleusOutputParticleLAB.e(), "Energy" ) <<
           std::endl;
    out << "                               (SCM) : " <<
           data.nucleusOutputParticle->GetParticleName() << " " <<
           G4BestUnit( data.nucleusOutputParticleSCM, "Energy" ) << " -- " <<
           G4BestUnit( data.nucleusOutputParticleSCM.e(), "Energy" ) <<
           std::endl;

    out.precision( prec );
    out.flags( savedFlags );

    return out;
}

