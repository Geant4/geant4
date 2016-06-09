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
 *       Filename:  CexmcSimpleDecayTableStore.hh
 *
 *    Description:  decay table serialization helper
 *
 *        Version:  1.0
 *        Created:  24.12.2009 14:17:40
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_SIMPLE_DECAY_TABLE_STORE_HH
#define CEXMC_SIMPLE_DECAY_TABLE_STORE_HH

#ifdef CEXMC_USE_PERSISTENCY

#include <boost/serialization/access.hpp>
#include <boost/serialization/map.hpp>
#include <G4Types.hh>

class  G4DecayTable;


typedef std::map< G4int, G4double >  CexmcDecayBranchesStore;


class  CexmcSimpleDecayTableStore
{
    friend class  boost::serialization::access;

    public:
        CexmcSimpleDecayTableStore();

        CexmcSimpleDecayTableStore( const G4DecayTable *  decayTable );

    public:
        const CexmcDecayBranchesStore &  GetDecayBranches( void ) const;

    private:
        template  < typename  Archive >
        void  serialize( Archive &  archive, const unsigned int  version );

    private:
        CexmcDecayBranchesStore  decayBranches;
};


inline const CexmcDecayBranchesStore &
                    CexmcSimpleDecayTableStore::GetDecayBranches( void ) const
{
    return decayBranches;
}


template  < typename  Archive >
void  CexmcSimpleDecayTableStore::serialize( Archive &  archive,
                                             const unsigned int )
{
    archive & decayBranches;
}

#endif

#endif

