/*
 * ============================================================================
 *
 *       Filename:  CexmcSimpleLorentzVectorStore.cc
 *
 *    Description:  G4LorentzVector serialization helper
 *
 *        Version:  1.0
 *        Created:  02.01.2010 14:12:56
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#ifdef CEXMC_USE_PERSISTENCY

#include "CexmcSimpleLorentzVectorStore.hh"


CexmcSimpleLorentzVectorStore::CexmcSimpleLorentzVectorStore()
{
}


CexmcSimpleLorentzVectorStore::CexmcSimpleLorentzVectorStore(
                                        const G4LorentzVector &  lorentzVector )
{
    px = lorentzVector.px();
    py = lorentzVector.py();
    pz = lorentzVector.pz();
    e = lorentzVector.e();
}


CexmcSimpleLorentzVectorStore::operator G4LorentzVector() const
{
    return G4LorentzVector( px, py, pz, e );
}

#endif

