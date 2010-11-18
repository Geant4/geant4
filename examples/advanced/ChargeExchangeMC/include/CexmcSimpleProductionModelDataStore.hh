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

