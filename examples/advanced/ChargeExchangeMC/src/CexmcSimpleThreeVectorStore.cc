/*
 * ============================================================================
 *
 *       Filename:  CexmcSimpleThreeVectorStore.cc
 *
 *    Description:  G4ThreeVector serialization helper
 *
 *        Version:  1.0
 *        Created:  24.12.2009 22:55:10
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#ifdef CEXMC_USE_PERSISTENCY

#include "CexmcSimpleThreeVectorStore.hh"


CexmcSimpleThreeVectorStore::CexmcSimpleThreeVectorStore()
{
}


CexmcSimpleThreeVectorStore::CexmcSimpleThreeVectorStore(
                                            const G4ThreeVector &  threeVector )
{
    x = threeVector.x();
    y = threeVector.y();
    z = threeVector.z();
}


CexmcSimpleThreeVectorStore::operator G4ThreeVector() const
{
    return G4ThreeVector( x, y, z );
}

#endif

