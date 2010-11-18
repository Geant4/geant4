/*
 * =============================================================================
 *
 *       Filename:  CexmcSimpleLorentzVectorStore.hh
 *
 *    Description:  G4LorentzVector serialization helper
 *
 *        Version:  1.0
 *        Created:  02.01.2010 14:08:42
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_SIMPLE_LORENTZ_VECTOR_STORE_HH
#define CEXMC_SIMPLE_LORENTZ_VECTOR_STORE_HH

#include <boost/serialization/access.hpp>
#include <G4LorentzVector.hh>


class  CexmcSimpleLorentzVectorStore
{
    friend class boost::serialization::access;
#ifdef CEXMC_USE_CUSTOM_FILTER
    friend class  CexmcASTEval;
#endif

    public:
        CexmcSimpleLorentzVectorStore();

        CexmcSimpleLorentzVectorStore( const G4LorentzVector &  lorentzVector );

    public:
        operator G4LorentzVector() const;

    private:
        template  < typename  Archive >
        void  serialize( Archive &  archive, const unsigned int  version );

    private:
        G4double  px;

        G4double  py;

        G4double  pz;

        G4double  e;
};


template  < typename  Archive >
void  CexmcSimpleLorentzVectorStore::serialize( Archive &  archive,
                                                const unsigned int )
{
    archive & px;
    archive & py;
    archive & pz;
    archive & e;
}


#endif

