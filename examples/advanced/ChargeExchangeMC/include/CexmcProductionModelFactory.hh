/*
 * =============================================================================
 *
 *       Filename:  CexmcProductionModelFactory.hh
 *
 *    Description:  production model factory
 *
 *        Version:  1.0
 *        Created:  03.11.2009 23:20:35
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_PRODUCTION_MODEL_FACTORY_HH
#define CEXMC_PRODUCTION_MODEL_FACTORY_HH

#include <G4VUserPhysicsList.hh>
#include <G4PionZero.hh>
#include <G4Eta.hh>
#include "CexmcPhysicsList.hh"
#include "CexmcCommon.hh"


namespace  CexmcPrivate
{
    template  < typename  BasePhysics >
    class  CexmcBasePhysicsInstance
    {
        public:
            static const CexmcBasePhysicsUsed  value = CexmcNoBasePhysics;
    };

#ifdef CEXMC_USE_QGSP_BIC_EMY
    template<>
    class CexmcBasePhysicsInstance< QGSP_BIC_EMY >
    {
        public:
            static const CexmcBasePhysicsUsed  value = Cexmc_QGSP_BIC_EMY;
    };
#else
    template<>
    class CexmcBasePhysicsInstance< QGSP_BERT >
    {
        public:
            static const CexmcBasePhysicsUsed  value = Cexmc_QGSP_BERT;
    };
#endif
}


template  < typename  BasePhysics,
            template  < typename > class  StudiedPhysics,
            template  < typename > class  ProductionModel >
class  CexmcProductionModelFactory
{
    public:
        static G4VUserPhysicsList *    Create(
                                CexmcProductionModelType  productionModelType );

        static CexmcBasePhysicsUsed    GetBasePhysics( void );

    private:
        CexmcProductionModelFactory();
};


template  < typename  BasePhysics,
            template  < typename > class  StudiedPhysics,
            template  < typename > class  ProductionModel >
G4VUserPhysicsList *  CexmcProductionModelFactory<
                                BasePhysics, StudiedPhysics, ProductionModel >::
                    Create( CexmcProductionModelType  productionModelType )
{
    switch ( productionModelType )
    {
    case CexmcPionZeroProduction :
        return new CexmcPhysicsList< BasePhysics, StudiedPhysics,
                                                ProductionModel< G4PionZero > >;
    case CexmcEtaProduction :
        return new CexmcPhysicsList< BasePhysics, StudiedPhysics,
                                                ProductionModel< G4Eta > >;
    default :
        return NULL;
    }
}


template  < typename  BasePhysics,
            template  < typename > class  StudiedPhysics,
            template  < typename > class  ProductionModel >
CexmcBasePhysicsUsed  CexmcProductionModelFactory<
                                BasePhysics, StudiedPhysics, ProductionModel >::
            GetBasePhysics( void )
{
    return CexmcPrivate::CexmcBasePhysicsInstance< BasePhysics >::value;
}


#endif

