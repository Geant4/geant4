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
 *       Filename:  CexmcEnergyDepositStore.hh
 *
 *    Description:  store energy deposit data and const references to
 *                  energy deposit collections in calorimeters
 *
 *        Version:  1.0
 *        Created:  25.11.2009 15:32:51
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_ENERGY_DEPOSIT_STORE_HH
#define CEXMC_ENERGY_DEPOSIT_STORE_HH

#include <G4Allocator.hh>
#include "CexmcCommon.hh"


struct  CexmcEnergyDepositStore
{
    CexmcEnergyDepositStore( G4double  monitorED_,
                             G4double  vetoCounterEDLeft_,
                             G4double  vetoCounterEDRight_,
                             G4double  calorimeterEDLeft_,
                             G4double  calorimeterEDRight_,
                             G4int     calorimeterEDLeftMaxX_,
                             G4int     calorimeterEDLeftMaxY_,
                             G4int     calorimeterEDRightMaxX_,
                             G4int     calorimeterEDRightMaxY_,
                             const CexmcEnergyDepositCalorimeterCollection &
                                       calorimeterEDLeftCollection_,
                             const CexmcEnergyDepositCalorimeterCollection &
                                       calorimeterEDRightCollection_ ) :
        monitorED( monitorED_ ), vetoCounterEDLeft( vetoCounterEDLeft_ ),
        vetoCounterEDRight( vetoCounterEDRight_ ),
        calorimeterEDLeft( calorimeterEDLeft_ ),
        calorimeterEDRight( calorimeterEDRight_ ),
        calorimeterEDLeftMaxX( calorimeterEDLeftMaxX_ ),
        calorimeterEDLeftMaxY( calorimeterEDLeftMaxY_ ),
        calorimeterEDRightMaxX( calorimeterEDRightMaxX_ ),
        calorimeterEDRightMaxY( calorimeterEDRightMaxY_ ),
        calorimeterEDLeftCollection( calorimeterEDLeftCollection_ ),
        calorimeterEDRightCollection( calorimeterEDRightCollection_ )
    {}

    void *  operator new( size_t  size );

    void    operator delete( void *  obj );

    G4double  monitorED;

    G4double  vetoCounterEDLeft;

    G4double  vetoCounterEDRight;

    G4double  calorimeterEDLeft;

    G4double  calorimeterEDRight;

    G4int     calorimeterEDLeftMaxX;

    G4int     calorimeterEDLeftMaxY;

    G4int     calorimeterEDRightMaxX;

    G4int     calorimeterEDRightMaxY;

    const CexmcEnergyDepositCalorimeterCollection &
              calorimeterEDLeftCollection;

    const CexmcEnergyDepositCalorimeterCollection &
              calorimeterEDRightCollection;
};


extern G4Allocator< CexmcEnergyDepositStore >  energyDepositStoreAllocator;


inline void *  CexmcEnergyDepositStore::operator new( size_t )
{
    return energyDepositStoreAllocator.MallocSingle();
}


inline void  CexmcEnergyDepositStore::operator delete( void *  obj )
{
    energyDepositStoreAllocator.FreeSingle(
                        reinterpret_cast< CexmcEnergyDepositStore * >( obj ) );
}


#endif

