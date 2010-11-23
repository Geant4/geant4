/*
 * =============================================================================
 *
 *       Filename:  CexmcSimpleThreeVectorStore.hh
 *
 *    Description:  G4ThreeVector serialization helper
 *
 *        Version:  1.0
 *        Created:  24.12.2009 22:45:36
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_SIMPLE_THREE_VECTOR_STORE_HH
#define CEXMC_SIMPLE_THREE_VECTOR_STORE_HH

#ifdef CEXMC_USE_PERSISTENCY

#include <boost/serialization/access.hpp>
#include <G4ThreeVector.hh>


class  CexmcSimpleThreeVectorStore
{
    friend class  boost::serialization::access;
#ifdef CEXMC_USE_CUSTOM_FILTER
    friend class  CexmcASTEval;
#endif

    public:
        CexmcSimpleThreeVectorStore();

        CexmcSimpleThreeVectorStore( const G4ThreeVector &  threeVector );

    public:
        operator G4ThreeVector() const;

    private:
        template  < typename  Archive >
        void  serialize( Archive &  archive, const unsigned int  version );

    private:
        G4double  x;

        G4double  y;

        G4double  z;
};


template  < typename  Archive >
void  CexmcSimpleThreeVectorStore::serialize( Archive &  archive,
                                              const unsigned int )
{
    archive & x;
    archive & y;
    archive & z;
}

#endif

#endif

