/*
 * =============================================================================
 *
 *       Filename:  CexmcASTEval.hh
 *
 *    Description:  abstract syntax tree for custom filter eval
 *
 *        Version:  1.0
 *        Created:  17.07.2010 15:43:09
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_AST_EVAL_HH
#define CEXMC_AST_EVAL_HH

#ifdef CEXMC_USE_CUSTOM_FILTER

#include <map>
#include <string>
#include <boost/variant/variant.hpp>
#include "CexmcAST.hh"
#include "CexmcEventSObject.hh"
#include "CexmcEventFastSObject.hh"
#include "CexmcException.hh"
#include "CexmcCommon.hh"


class  CexmcASTEval : public CexmcAST::BasicEval
{
    private:
        typedef boost::variant< const CexmcEnergyDepositCalorimeterCollection *,
                                const bool * >     VarAddr;

        typedef std::map< std::string, VarAddr >   VarAddrMap;

        typedef std::pair< std::string, VarAddr >  VarAddrMapData;

    public:
        explicit CexmcASTEval(
                        const CexmcEventFastSObject *  evFastSObject = NULL,
                        const CexmcEventSObject *  evSObject = NULL );

    public:
        void  SetAddressedData(
                        const CexmcEventFastSObject *  evFastSObject_ = NULL,
                        const CexmcEventSObject *  evSObject_ = NULL );

        void  BindAddresses( CexmcAST::Subtree &  ast );

        void  ResetAddressBinding( CexmcAST::Subtree &  ast );

    private:
        ScalarValueType  GetFunScalarValue( const CexmcAST::Subtree &  ast )
                                                                        const;

        ScalarValueType  GetVarScalarValue( const CexmcAST::Variable &  var )
                                                                        const;

        void             GetEDCollectionValue( const CexmcAST::Node &  node,
                    CexmcEnergyDepositCalorimeterCollection &  edCol ) const;

    private:
        const G4double *  GetThreeVectorElementAddrByIndex(
                                    const CexmcSimpleThreeVectorStore &  vect,
                                    G4int  index ) const;

        const G4double *  GetLorentzVectorElementAddrByIndex(
                                    const CexmcSimpleLorentzVectorStore &  vect,
                                    G4int  index ) const;

    private:
        const CexmcEventFastSObject *  evFastSObject;

        const CexmcEventSObject *      evSObject;

    private:
        VarAddrMap                     varAddrMap;

    private:
        static const G4double          constants[];
};


inline void  CexmcASTEval::SetAddressedData(
                                const CexmcEventFastSObject *  evFastSObject_,
                                const CexmcEventSObject *  evSObject_ )
{
    varAddrMap.clear();
    evFastSObject = evFastSObject_;
    evSObject = evSObject_;
}


inline const G4double *  CexmcASTEval::GetThreeVectorElementAddrByIndex(
                                const CexmcSimpleThreeVectorStore &  vect,
                                G4int  index ) const
{
    switch ( index )
    {
    case 1 :
        return &vect.x;
    case 2 :
        return &vect.y;
    case 3 :
        return &vect.z;
    default :
        throw CexmcException( CexmcCFUnexpectedVectorIndex );
        return NULL;
    }
}


inline const G4double *  CexmcASTEval::GetLorentzVectorElementAddrByIndex(
                                const CexmcSimpleLorentzVectorStore &  vect,
                                G4int  index ) const
{
    switch ( index )
    {
    case 1 :
        return &vect.px;
    case 2 :
        return &vect.py;
    case 3 :
        return &vect.pz;
    case 4 :
        return &vect.e;
    default :
        throw CexmcException( CexmcCFUnexpectedVectorIndex );
        return NULL;
    }
}

#endif

#endif

