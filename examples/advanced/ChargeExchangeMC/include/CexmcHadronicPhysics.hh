/*
 * =============================================================================
 *
 *       Filename:  CexmcHadronicPhysics.hh
 *
 *    Description:  hadronic physics with adjustable production model
 *
 *        Version:  1.0
 *        Created:  28.10.2009 17:16:44
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_HADRONIC_PHYSICS_HH
#define CEXMC_HADRONIC_PHYSICS_HH

#include "CexmcStudiedPhysics.hh"
#include "CexmcHadronicProcess.hh"


template  < typename  ProductionModel >
class  CexmcHadronicPhysics : public CexmcStudiedPhysics< CexmcHadronicProcess >
{
    public:
        explicit CexmcHadronicPhysics( CexmcPhysicsManager *  physicsManager );

        ~CexmcHadronicPhysics();

    private:
        void  ApplyInteractionModel( G4VProcess *  process );
};


template  < typename  ProductionModel >
CexmcHadronicPhysics< ProductionModel >::CexmcHadronicPhysics(
                                    CexmcPhysicsManager *  physicsManager ) :
    CexmcStudiedPhysics< CexmcHadronicProcess >( physicsManager )
{
    productionModel = new ProductionModel;
}


template  < typename  ProductionModel >
CexmcHadronicPhysics< ProductionModel >::~CexmcHadronicPhysics()
{
    delete productionModel;
}


template  < typename  ProductionModel >
void  CexmcHadronicPhysics< ProductionModel >::ApplyInteractionModel(
                                                        G4VProcess *  process )
{
    CexmcHadronicProcess *  theProcess( static_cast< CexmcHadronicProcess * >(
                                                                    process ) );
    theProcess->RegisterProductionModel( productionModel );
}


#endif

