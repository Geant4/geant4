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

