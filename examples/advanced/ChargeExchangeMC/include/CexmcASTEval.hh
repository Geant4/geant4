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
                                const bool * >    VarAddr;

        typedef std::map< std::string, VarAddr >  VarAddrMap;

        typedef VarAddrMap::value_type            VarAddrMapData;

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

