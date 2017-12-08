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
                      const G4LorentzVector &  incidentParticleSCM_,
                      const G4LorentzVector &  incidentParticleLAB_,
                      const G4LorentzVector &  nucleusParticleSCM_,
                      const G4LorentzVector &  nucleusParticleLAB_,
                      const G4LorentzVector &  outputParticleSCM_,
                      const G4LorentzVector &  outputParticleLAB_,
                      const G4LorentzVector &  nucleusOutputParticleSCM_,
                      const G4LorentzVector &  nucleusOutputParticleLAB_,
                      const G4ParticleDefinition *  incidentParticle_,
                      const G4ParticleDefinition *  nucleusParticle_,
                      const G4ParticleDefinition *  outputParticle_,
                      const G4ParticleDefinition *  nucleusOutputParticle_ ) :
    incidentParticleSCM( incidentParticleSCM_ ),
    incidentParticleLAB( incidentParticleLAB_ ),
    nucleusParticleSCM( nucleusParticleSCM_ ),
    nucleusParticleLAB( nucleusParticleLAB_ ),
    outputParticleSCM( outputParticleSCM_ ),
    outputParticleLAB( outputParticleLAB_ ),
    nucleusOutputParticleSCM( nucleusOutputParticleSCM_ ),
    nucleusOutputParticleLAB( nucleusOutputParticleLAB_ ),
    incidentParticle( incidentParticle_ ), nucleusParticle( nucleusParticle_ ),
    outputParticle( outputParticle_ ),
    nucleusOutputParticle( nucleusOutputParticle_ )
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

