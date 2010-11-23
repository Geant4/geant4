/*
 * ============================================================================
 *
 *       Filename:  CexmcSimpleDecayTableStore.cc
 *
 *    Description:  decay table serialization helper
 *
 *        Version:  1.0
 *        Created:  24.12.2009 15:48:13
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#ifdef CEXMC_USE_PERSISTENCY

#include <G4DecayTable.hh>
#include "CexmcSimpleDecayTableStore.hh"


CexmcSimpleDecayTableStore::CexmcSimpleDecayTableStore()
{
}


CexmcSimpleDecayTableStore::CexmcSimpleDecayTableStore(
                                            const G4DecayTable *  decayTable )
{
    for ( G4int  i( 0 ); i < decayTable->entries(); ++i )
    {
        decayBranches.insert(
                std::pair< G4int, G4double >( i,
                            decayTable->GetDecayChannel( i )->GetBR() ) );
    }
}

#endif

