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
                                    CexmcPhysicsManager *  physicsManager_ ) :
    CexmcStudiedPhysics< CexmcHadronicProcess >( physicsManager_ )
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

