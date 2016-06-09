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
 *       Filename:  CexmcSimpleProductionModelDataStore.hh
 *
 *    Description:  serialization helper for production model data
 *
 *        Version:  1.0
 *        Created:  02.01.2010 14:37:16
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_SIMPLE_PRODUCTION_MODEL_DATA_STORE_HH
#define CEXMC_SIMPLE_PRODUCTION_MODEL_DATA_STORE_HH

#ifdef CEXMC_USE_PERSISTENCY

#include <boost/serialization/access.hpp>
#include "CexmcSimpleLorentzVectorStore.hh"

class  CexmcProductionModelData;


class  CexmcSimpleProductionModelDataStore
{
    friend class  boost::serialization::access;
#ifdef CEXMC_USE_CUSTOM_FILTER
    friend class  CexmcASTEval;
#endif

    public:
        CexmcSimpleProductionModelDataStore();

        CexmcSimpleProductionModelDataStore(
                                    const CexmcProductionModelData &  pmData );

    public:
        operator CexmcProductionModelData() const;

    private:
        template  < typename  Archive >
        void  serialize( Archive &  archive, const unsigned int  version );

    private:
        CexmcSimpleLorentzVectorStore  incidentParticleSCM;

        CexmcSimpleLorentzVectorStore  incidentParticleLAB;

        CexmcSimpleLorentzVectorStore  nucleusParticleSCM;

        CexmcSimpleLorentzVectorStore  nucleusParticleLAB;

        CexmcSimpleLorentzVectorStore  outputParticleSCM;

        CexmcSimpleLorentzVectorStore  outputParticleLAB;

        CexmcSimpleLorentzVectorStore  nucleusOutputParticleSCM;

        CexmcSimpleLorentzVectorStore  nucleusOutputParticleLAB;

        G4int                          incidentParticle;

        G4int                          nucleusParticle;

        G4int                          outputParticle;

        G4int                          nucleusOutputParticle;
};


template  < typename  Archive >
void  CexmcSimpleProductionModelDataStore::serialize( Archive &  archive,
                                                      const unsigned int )
{
    archive & incidentParticleSCM;
    archive & incidentParticleLAB;
    archive & nucleusParticleSCM;
    archive & nucleusParticleLAB;
    archive & outputParticleSCM;
    archive & outputParticleLAB;
    archive & nucleusOutputParticleSCM;
    archive & nucleusOutputParticleLAB;
    archive & incidentParticle;
    archive & nucleusParticle;
    archive & outputParticle;
    archive & nucleusOutputParticle;
}

#endif

#endif

