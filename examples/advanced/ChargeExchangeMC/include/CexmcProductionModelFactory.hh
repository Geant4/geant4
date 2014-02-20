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
    template  < typename >
    struct  CexmcBasePhysicsInstance
    {
        static const CexmcBasePhysicsUsed  value = CexmcNoBasePhysics;
    };

#ifdef CEXMC_USE_QGSP_BIC_EMY
    template  <>
    struct  CexmcBasePhysicsInstance< QGSP_BIC_EMY >
    {
        static const CexmcBasePhysicsUsed  value = Cexmc_QGSP_BIC_EMY;
    };
#else
#ifdef CEXMC_USE_QGSP_BERT
    template  <>
    struct  CexmcBasePhysicsInstance< QGSP_BERT >
    {
        static const CexmcBasePhysicsUsed  value = Cexmc_QGSP_BERT;
    };
#else
    template  <>
    struct  CexmcBasePhysicsInstance< FTFP_BERT >
    {
        static const CexmcBasePhysicsUsed  value = Cexmc_FTFP_BERT;
    };
#endif
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

