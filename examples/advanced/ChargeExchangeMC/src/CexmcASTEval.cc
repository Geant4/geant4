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
 *       Filename:  CexmcASTEval.cc
 *
 *    Description:  abstract syntax tree for custom filter eval
 *
 *        Version:  1.0
 *        Created:  17.07.2010 15:46:01
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifdef CEXMC_USE_CUSTOM_FILTER

#include <numeric>
#include <boost/variant/get.hpp>
#include <G4SystemOfUnits.hh>
#include "CexmcASTEval.hh"


namespace
{
    const std::string  CexmcCFVarEvent( "event" );
    const std::string  CexmcCFVarOpCosThetaSCM( "op_cosTh_SCM" );
    const std::string  CexmcCFVarEDT( "edt" );
    const std::string  CexmcCFVarTPT( "tpt" );
    const std::string  CexmcCFVarMon( "mon" );
    const std::string  CexmcCFVarMonED( "monED" );
    const std::string  CexmcCFVarVclED( "vclED" );
    const std::string  CexmcCFVarVcrED( "vcrED" );
    const std::string  CexmcCFVarClED( "clED" );
    const std::string  CexmcCFVarCrED( "crED" );
    const std::string  CexmcCFVarClEDCol( "clEDcol" );
    const std::string  CexmcCFVarCrEDCol( "crEDcol" );
    const std::string  CexmcCFVarBpMonPosL( "bp_mon_posl" );
    const std::string  CexmcCFVarBpMonPosW( "bp_mon_posw" );
    const std::string  CexmcCFVarBpMonDirL( "bp_mon_dirl" );
    const std::string  CexmcCFVarBpMonDirW( "bp_mon_dirw" );
    const std::string  CexmcCFVarBpMonMom( "bp_mon_mom" );
    const std::string  CexmcCFVarBpMonTid( "bp_mon_tid" );
    const std::string  CexmcCFVarBpTgtPosL( "bp_tgt_posl" );
    const std::string  CexmcCFVarBpTgtPosW( "bp_tgt_posw" );
    const std::string  CexmcCFVarBpTgtDirL( "bp_tgt_dirl" );
    const std::string  CexmcCFVarBpTgtDirW( "bp_tgt_dirw" );
    const std::string  CexmcCFVarBpTgtMom( "bp_tgt_mom" );
    const std::string  CexmcCFVarBpTgtTid( "bp_tgt_tid" );
    const std::string  CexmcCFVarOpTgtPosL( "op_tgt_posl" );
    const std::string  CexmcCFVarOpTgtPosW( "op_tgt_posw" );
    const std::string  CexmcCFVarOpTgtDirL( "op_tgt_dirl" );
    const std::string  CexmcCFVarOpTgtDirW( "op_tgt_dirw" );
    const std::string  CexmcCFVarOpTgtMom( "op_tgt_mom" );
    const std::string  CexmcCFVarOpTgtTid( "op_tgt_tid" );
    const std::string  CexmcCFVarNpTgtPosL( "np_tgt_posl" );
    const std::string  CexmcCFVarNpTgtPosW( "np_tgt_posw" );
    const std::string  CexmcCFVarNpTgtDirL( "np_tgt_dirl" );
    const std::string  CexmcCFVarNpTgtDirW( "np_tgt_dirw" );
    const std::string  CexmcCFVarNpTgtMom( "np_tgt_mom" );
    const std::string  CexmcCFVarNpTgtTid( "np_tgt_tid" );
    const std::string  CexmcCFVarOpdp1TgtPosL( "opdp1_tgt_posl" );
    const std::string  CexmcCFVarOpdp1TgtPosW( "opdp1_tgt_posw" );
    const std::string  CexmcCFVarOpdp1TgtDirL( "opdp1_tgt_dirl" );
    const std::string  CexmcCFVarOpdp1TgtDirW( "opdp1_tgt_dirw" );
    const std::string  CexmcCFVarOpdp1TgtMom( "opdp1_tgt_mom" );
    const std::string  CexmcCFVarOpdp1TgtTid( "opdp1_tgt_tid" );
    const std::string  CexmcCFVarOpdp2TgtPosL( "opdp2_tgt_posl" );
    const std::string  CexmcCFVarOpdp2TgtPosW( "opdp2_tgt_posw" );
    const std::string  CexmcCFVarOpdp2TgtDirL( "opdp2_tgt_dirl" );
    const std::string  CexmcCFVarOpdp2TgtDirW( "opdp2_tgt_dirw" );
    const std::string  CexmcCFVarOpdp2TgtMom( "opdp2_tgt_mom" );
    const std::string  CexmcCFVarOpdp2TgtTid( "opdp2_tgt_tid" );
    const std::string  CexmcCFVarOpdpVclPosL( "opdp_vcl_posl" );
    const std::string  CexmcCFVarOpdpVclPosW( "opdp_vcl_posw" );
    const std::string  CexmcCFVarOpdpVclDirL( "opdp_vcl_dirl" );
    const std::string  CexmcCFVarOpdpVclDirW( "opdp_vcl_dirw" );
    const std::string  CexmcCFVarOpdpVclMom( "opdp_vcl_mom" );
    const std::string  CexmcCFVarOpdpVclTid( "opdp_vcl_tid" );
    const std::string  CexmcCFVarOpdpVcrPosL( "opdp_vcr_posl" );
    const std::string  CexmcCFVarOpdpVcrPosW( "opdp_vcr_posw" );
    const std::string  CexmcCFVarOpdpVcrDirL( "opdp_vcr_dirl" );
    const std::string  CexmcCFVarOpdpVcrDirW( "opdp_vcr_dirw" );
    const std::string  CexmcCFVarOpdpVcrMom( "opdp_vcr_mom" );
    const std::string  CexmcCFVarOpdpVcrTid( "opdp_vcr_tid" );
    const std::string  CexmcCFVarOpdpClPosL( "opdp_cl_posl" );
    const std::string  CexmcCFVarOpdpClPosW( "opdp_cl_posw" );
    const std::string  CexmcCFVarOpdpClDirL( "opdp_cl_dirl" );
    const std::string  CexmcCFVarOpdpClDirW( "opdp_cl_dirw" );
    const std::string  CexmcCFVarOpdpClMom( "opdp_cl_mom" );
    const std::string  CexmcCFVarOpdpClTid( "opdp_cl_tid" );
    const std::string  CexmcCFVarOpdpCrPosL( "opdp_cr_posl" );
    const std::string  CexmcCFVarOpdpCrPosW( "opdp_cr_posw" );
    const std::string  CexmcCFVarOpdpCrDirL( "opdp_cr_dirl" );
    const std::string  CexmcCFVarOpdpCrDirW( "opdp_cr_dirw" );
    const std::string  CexmcCFVarOpdpCrMom( "opdp_cr_mom" );
    const std::string  CexmcCFVarOpdpCrTid( "opdp_cr_tid" );
    const std::string  CexmcCFVarIpSCM( "ipSCM" );
    const std::string  CexmcCFVarIpLAB( "ipLAB" );
    const std::string  CexmcCFVarNpSCM( "npSCM" );
    const std::string  CexmcCFVarNpLAB( "npLAB" );
    const std::string  CexmcCFVarOpSCM( "opSCM" );
    const std::string  CexmcCFVarOpLAB( "opLAB" );
    const std::string  CexmcCFVarNopSCM( "nopSCM" );
    const std::string  CexmcCFVarNopLAB( "nopLAB" );
    const std::string  CexmcCFVarIpId( "ipId" );
    const std::string  CexmcCFVarNpId( "npId" );
    const std::string  CexmcCFVarOpId( "opId" );
    const std::string  CexmcCFVarNopId( "nopId" );
    const std::string  CexmcCFVarConst_eV( "eV" );
    const std::string  CexmcCFVarConst_keV( "keV" );
    const std::string  CexmcCFVarConst_MeV( "MeV" );
    const std::string  CexmcCFVarConst_GeV( "GeV" );
    const std::string  CexmcCFVarConst_mm( "mm" );
    const std::string  CexmcCFVarConst_cm( "cm" );
    const std::string  CexmcCFVarConst_m( "m" );
}


const G4double  CexmcASTEval::constants[] = { eV, keV, MeV, GeV, mm, cm, m };


CexmcASTEval::CexmcASTEval( const CexmcEventFastSObject *  evFastSObject_,
                            const CexmcEventSObject *  evSObject_ ) :
    evFastSObject( evFastSObject_ ), evSObject( evSObject_ )
{
}


CexmcAST::BasicEval::ScalarValueType  CexmcASTEval::GetFunScalarValue(
                                        const CexmcAST::Subtree &  ast ) const
{
    const CexmcAST::Function &  fun( boost::get< CexmcAST::Function >(
                                                                ast.type ) );

    if ( fun == "Sum" )
    {
        CexmcEnergyDepositCalorimeterCollection  edCol;
        GetEDCollectionValue( ast.children[ 0 ], edCol );

        G4double  result( 0. );

        for ( CexmcEnergyDepositCalorimeterCollection::iterator
                                    k( edCol.begin() ); k != edCol.end(); ++k )
        {
            result += std::accumulate( k->begin(), k->end(), G4double( 0. ) );
        }

        return result;
    }

    bool             evalResult( false );
    ScalarValueType  result( GetBasicFunScalarValue( ast, evalResult ) );

    if ( evalResult )
        return result;

    throw CexmcException( CexmcCFUnexpectedFunction );

    return 0;
}


CexmcAST::BasicEval::ScalarValueType  CexmcASTEval::GetVarScalarValue(
                                        const CexmcAST::Variable &  var ) const
{
    if ( evFastSObject == NULL || evSObject == NULL )
        throw CexmcException( CexmcCFUninitialized );

    /* Variables with initialized address */

    /* bound to CexmcAST::Variable:addr */

    const double * const *  addr( boost::get< const double * >( &var.addr ) );

    if ( addr )
    {
        if ( *addr )
            return **addr;
    }
    else
    {
        const int * const &  addr_( boost::get< const int * >( var.addr ) );

        if ( addr_ )
            return *addr_;
    }

    /* found in varAddrMap */

    VarAddrMap::const_iterator  found( varAddrMap.find( var.name ) );

    if ( found != varAddrMap.end() )
    {
        const CexmcEnergyDepositCalorimeterCollection * const *  addr_(
                boost::get< const CexmcEnergyDepositCalorimeterCollection * >(
                                                            &found->second ) );
        if ( addr_ )
        {
            if ( *addr_ )
            {
                if ( ( *addr_ )->size() == 0 )
                    throw CexmcException( CexmcCFUninitializedVector );
                if ( var.index1 == 0 || var.index2 == 0 )
                    throw CexmcException( CexmcCFUnexpectedVectorIndex );
                return ( *addr_ )->at( var.index1 - 1 ).at( var.index2 - 1 );
            }
        }
        else
        {
            const bool * const &  addr__( boost::get< const bool * >(
                                                            found->second ) );
            if ( addr__ )
                return int( *addr__ );
        }
    }

    /* Variables without address */

    if ( var.name == CexmcCFVarTPT )
    {
        return int( evSObject->targetTPOutputParticle.trackId !=
                    CexmcInvalidTrackId );
    } 

    throw CexmcException( CexmcCFUnexpectedVariable );

    return 0;
}


void  CexmcASTEval::GetEDCollectionValue( const CexmcAST::Node &  node,
                        CexmcEnergyDepositCalorimeterCollection &  edCol ) const
{
    if ( evSObject == NULL )
        throw CexmcException( CexmcCFUninitialized );

    const CexmcAST::Subtree *  ast( boost::get< CexmcAST::Subtree >( &node ) );

    if ( ast )
    {
        const CexmcAST::Function &  fun( boost::get< CexmcAST::Function >(
                                                                ast->type ) );

        if ( fun == "Inner" )
        {
            GetEDCollectionValue( ast->children[ 0 ], edCol );
            edCol.pop_back();
            edCol.erase( edCol.begin() );
            for ( CexmcEnergyDepositCalorimeterCollection::iterator
                            k( edCol.begin() ); k != edCol.end(); ++k )
            {
                k->pop_back();
                k->erase( k->begin() );
            }
            return;
        }
        if ( fun == "Outer" )
        {
            GetEDCollectionValue( ast->children[ 0 ], edCol );
            if ( edCol.size() < 3 )
                return;
            for ( CexmcEnergyDepositCalorimeterCollection::iterator
                            k( edCol.begin() + 1 ); k != edCol.end() - 1; ++k )
            {
                if ( k->size() < 3 )
                    continue;
                k->erase( k->begin() + 1, k->end() - 1 );
            }
            return;
        }
    }
    else
    {
        const CexmcAST::Leaf &      leaf( boost::get< CexmcAST::Leaf >(
                                                                    node ) );
        const CexmcAST::Variable &  var( boost::get< CexmcAST::Variable >(
                                                                    leaf ) );

        if ( var.index1 != 0 || var.index2 != 0 )
            throw CexmcException( CexmcCFUnexpectedVariableUsage );

        VarAddrMap::const_iterator  found( varAddrMap.find( var.name ) );

        if ( found == varAddrMap.end() )
            throw CexmcException( CexmcCFUnexpectedVariable );

        const CexmcEnergyDepositCalorimeterCollection * const *  addr(
                boost::get< const CexmcEnergyDepositCalorimeterCollection * >(
                                                            &found->second ) );
        if ( ! addr )
        {
            throw CexmcException( CexmcCFUnexpectedVariableUsage );
        }
        else
        {
            if ( *addr )
                edCol = **addr;
            return;
        }
    }
}


void  CexmcASTEval::BindAddresses( CexmcAST::Subtree &  ast )
{
    if ( evFastSObject == NULL || evSObject == NULL )
        return;

    for ( std::vector< CexmcAST::Node >::iterator  k( ast.children.begin() );
                                                  k != ast.children.end(); ++k )
    {
        CexmcAST::Subtree *  subtree( boost::get< CexmcAST::Subtree >( &*k ) );

        if ( subtree )
        {
            BindAddresses( *subtree );
        }
        else
        {
            CexmcAST::Leaf &      leaf( boost::get< CexmcAST::Leaf >( *k ) );
            CexmcAST::Variable *  var( boost::get< CexmcAST::Variable >(
                                                                    &leaf ) );
            if ( ! var )
                continue;

            const int * const *   intVarAddr(
                                    boost::get< const int * >( &var->addr ) );
            if ( intVarAddr )
            {
                if ( *intVarAddr )
                    continue;
            }
            else
            {
                const double * const &  doubleVarAddr(
                                    boost::get< const double * >( var->addr ) );
                if ( doubleVarAddr )
                    continue;
            }

            VarAddrMap::const_iterator  found( varAddrMap.find( var->name ) );

            if ( found != varAddrMap.end() )
                continue;

            do
            {
                if ( var->name == CexmcCFVarEvent )
                {
                    var->addr = &evFastSObject->eventId;
                    break;
                }
                if ( var->name == CexmcCFVarOpCosThetaSCM )
                {
                    var->addr = &evFastSObject->opCosThetaSCM;
                    break;
                }
                if ( var->name == CexmcCFVarEDT )
                {
                    varAddrMap.insert( VarAddrMapData( var->name,
                                   &evFastSObject->edDigitizerHasTriggered ) );
                    break;
                } 
                if ( var->name == CexmcCFVarMon )
                {
                    varAddrMap.insert( VarAddrMapData( var->name,
                           &evFastSObject->edDigitizerMonitorHasTriggered ) );
                    break;
                }
                if ( var->name == CexmcCFVarMonED )
                {
                    var->addr = &evSObject->monitorED;
                    break;
                }
                if ( var->name == CexmcCFVarVclED )
                {
                    var->addr = &evSObject->vetoCounterEDLeft;
                    break;
                }
                if ( var->name == CexmcCFVarVcrED )
                {
                    var->addr = &evSObject->vetoCounterEDRight;
                    break;
                }
                if ( var->name == CexmcCFVarClED )
                {
                    var->addr = &evSObject->calorimeterEDLeft;
                    break;
                }
                if ( var->name == CexmcCFVarCrED )
                {
                    var->addr = &evSObject->calorimeterEDRight;
                    break;
                }
                if ( var->name == CexmcCFVarClEDCol )
                {
                    varAddrMap.insert( VarAddrMapData( var->name,
                                &evSObject->calorimeterEDLeftCollection ) );
                    break;
                }
                if ( var->name == CexmcCFVarCrEDCol )
                {
                    varAddrMap.insert( VarAddrMapData( var->name,
                                &evSObject->calorimeterEDRightCollection ) );
                    break;
                }
                if ( var->name == CexmcCFVarBpMonPosL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->monitorTP.positionLocal, var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarBpMonPosW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->monitorTP.positionWorld, var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarBpMonDirL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->monitorTP.directionLocal, var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarBpMonDirW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->monitorTP.directionWorld, var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarBpMonMom )
                {
                    var->addr = &evSObject->monitorTP.momentumAmp;
                    break;
                }
                if ( var->name == CexmcCFVarBpMonTid )
                {
                    var->addr = &evSObject->monitorTP.trackId;
                    break;
                }
                if ( var->name == CexmcCFVarBpTgtPosL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->targetTPBeamParticle.positionLocal,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarBpTgtPosW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->targetTPBeamParticle.positionWorld,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarBpTgtDirL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->targetTPBeamParticle.directionLocal,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarBpTgtDirW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->targetTPBeamParticle.directionWorld,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarBpTgtMom )
                {
                    var->addr = &evSObject->targetTPBeamParticle.momentumAmp;
                    break;
                }
                if ( var->name == CexmcCFVarBpTgtTid )
                {
                    var->addr = &evSObject->targetTPBeamParticle.trackId;
                    break;
                }
                if ( var->name == CexmcCFVarOpTgtPosL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->targetTPOutputParticle.positionLocal,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpTgtPosW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->targetTPOutputParticle.positionWorld,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpTgtDirL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->targetTPOutputParticle.directionLocal,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpTgtDirW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->targetTPOutputParticle.directionWorld,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpTgtMom )
                {
                    var->addr = &evSObject->targetTPOutputParticle.momentumAmp;
                    break;
                }
                if ( var->name == CexmcCFVarOpTgtTid )
                {
                    var->addr = &evSObject->targetTPOutputParticle.trackId;
                    break;
                }
                if ( var->name == CexmcCFVarNpTgtPosL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->targetTPNucleusParticle.positionLocal,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarNpTgtPosW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->targetTPNucleusParticle.positionWorld,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarNpTgtDirL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->targetTPNucleusParticle.directionLocal,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarNpTgtDirW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->targetTPNucleusParticle.directionWorld,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarNpTgtMom )
                {
                    var->addr = &evSObject->targetTPNucleusParticle.momentumAmp;
                    break;
                }
                if ( var->name == CexmcCFVarNpTgtTid )
                {
                    var->addr = &evSObject->targetTPNucleusParticle.trackId;
                    break;
                }
                if ( var->name == CexmcCFVarOpdp1TgtPosL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                        evSObject->targetTPOutputParticleDecayProductParticle1.
                                                                positionLocal,
                        var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdp1TgtPosW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                        evSObject->targetTPOutputParticleDecayProductParticle1.
                                                                positionWorld,
                        var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdp1TgtDirL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                        evSObject->targetTPOutputParticleDecayProductParticle1.
                                                                directionLocal,
                        var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdp1TgtDirW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                        evSObject->targetTPOutputParticleDecayProductParticle1.
                                                                directionWorld,
                        var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdp1TgtMom )
                {
                    var->addr = &evSObject->
                        targetTPOutputParticleDecayProductParticle1.momentumAmp;
                    break;
                }
                if ( var->name == CexmcCFVarOpdp1TgtTid )
                {
                    var->addr = &evSObject->
                        targetTPOutputParticleDecayProductParticle1.trackId;
                    break;
                }
                if ( var->name == CexmcCFVarOpdp2TgtPosL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                        evSObject->targetTPOutputParticleDecayProductParticle2.
                                                                positionLocal,
                        var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdp2TgtPosW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                        evSObject->targetTPOutputParticleDecayProductParticle2.
                                                                positionWorld,
                        var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdp2TgtDirL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                        evSObject->targetTPOutputParticleDecayProductParticle2.
                                                                directionLocal,
                        var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdp2TgtDirW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                        evSObject->targetTPOutputParticleDecayProductParticle2.
                                                                directionWorld,
                        var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdp2TgtMom )
                {
                    var->addr = &evSObject->
                        targetTPOutputParticleDecayProductParticle2.momentumAmp;
                    break;
                }
                if ( var->name == CexmcCFVarOpdp2TgtTid )
                {
                    var->addr = &evSObject->
                        targetTPOutputParticleDecayProductParticle2.trackId;
                    break;
                }
                if ( var->name == CexmcCFVarOpdpVclPosL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->vetoCounterTPLeft.positionLocal,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdpVclPosW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->vetoCounterTPLeft.positionWorld,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdpVclDirL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->vetoCounterTPLeft.directionLocal,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdpVclDirW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->vetoCounterTPLeft.directionWorld,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdpVclMom )
                {
                    var->addr = &evSObject->vetoCounterTPLeft.momentumAmp;
                    break;
                }
                if ( var->name == CexmcCFVarOpdpVclTid )
                {
                    var->addr = &evSObject->vetoCounterTPLeft.trackId;
                    break;
                }
                if ( var->name == CexmcCFVarOpdpVcrPosL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->vetoCounterTPRight.positionLocal,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdpVcrPosW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->vetoCounterTPRight.positionWorld,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdpVcrDirL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->vetoCounterTPRight.directionLocal,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdpVcrDirW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->vetoCounterTPRight.directionWorld,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdpVcrMom )
                {
                    var->addr = &evSObject->vetoCounterTPRight.momentumAmp;
                    break;
                }
                if ( var->name == CexmcCFVarOpdpVcrTid )
                {
                    var->addr = &evSObject->vetoCounterTPRight.trackId;
                    break;
                }
                if ( var->name == CexmcCFVarOpdpClPosL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->calorimeterTPLeft.positionLocal,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdpClPosW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->calorimeterTPLeft.positionWorld,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdpClDirL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->calorimeterTPLeft.directionLocal,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdpClDirW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->calorimeterTPLeft.directionWorld,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdpClMom )
                {
                    var->addr = &evSObject->calorimeterTPLeft.momentumAmp;
                    break;
                }
                if ( var->name == CexmcCFVarOpdpClTid )
                {
                    var->addr = &evSObject->calorimeterTPLeft.trackId;
                    break;
                }
                if ( var->name == CexmcCFVarOpdpCrPosL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->calorimeterTPRight.positionLocal,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdpCrPosW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->calorimeterTPRight.positionWorld,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdpCrDirL )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->calorimeterTPRight.directionLocal,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdpCrDirW )
                {
                    var->addr = GetThreeVectorElementAddrByIndex(
                            evSObject->calorimeterTPRight.directionWorld,
                            var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpdpCrMom )
                {
                    var->addr = &evSObject->calorimeterTPRight.momentumAmp;
                    break;
                }
                if ( var->name == CexmcCFVarOpdpCrTid )
                {
                    var->addr = &evSObject->calorimeterTPRight.trackId;
                    break;
                }
                if ( var->name == CexmcCFVarIpSCM )
                {
                    var->addr = GetLorentzVectorElementAddrByIndex(
                        evSObject->productionModelData.incidentParticleSCM,
                        var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarIpLAB )
                {
                    var->addr = GetLorentzVectorElementAddrByIndex(
                        evSObject->productionModelData.incidentParticleLAB,
                        var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarNpSCM )
                {
                    var->addr = GetLorentzVectorElementAddrByIndex(
                        evSObject->productionModelData.nucleusParticleSCM,
                        var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarNpLAB )
                {
                    var->addr = GetLorentzVectorElementAddrByIndex(
                        evSObject->productionModelData.nucleusParticleLAB,
                        var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpSCM )
                {
                    var->addr = GetLorentzVectorElementAddrByIndex(
                        evSObject->productionModelData.outputParticleSCM,
                        var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarOpLAB )
                {
                    var->addr = GetLorentzVectorElementAddrByIndex(
                        evSObject->productionModelData.outputParticleLAB,
                        var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarNopSCM )
                {
                    var->addr = GetLorentzVectorElementAddrByIndex(
                        evSObject->productionModelData.nucleusOutputParticleSCM,
                        var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarNopLAB )
                {
                    var->addr = GetLorentzVectorElementAddrByIndex(
                        evSObject->productionModelData.nucleusOutputParticleLAB,
                        var->index1 );
                    break;
                }
                if ( var->name == CexmcCFVarIpId )
                {
                    var->addr =
                        &evSObject->productionModelData.incidentParticle;
                    break;
                }
                if ( var->name == CexmcCFVarNpId )
                {
                    var->addr = &evSObject->productionModelData.nucleusParticle;
                    break;
                }
                if ( var->name == CexmcCFVarOpId )
                {
                    var->addr = &evSObject->productionModelData.outputParticle;
                    break;
                }
                if ( var->name == CexmcCFVarNopId )
                {
                    var->addr =
                        &evSObject->productionModelData.nucleusOutputParticle;
                    break;
                }
                if ( var->name == CexmcCFVarConst_eV )
                {
                    var->addr = &constants[ 0 ];
                    break;
                }
                if ( var->name == CexmcCFVarConst_keV )
                {
                    var->addr = &constants[ 1 ];
                    break;
                }
                if ( var->name == CexmcCFVarConst_MeV )
                {
                    var->addr = &constants[ 2 ];
                    break;
                }
                if ( var->name == CexmcCFVarConst_GeV )
                {
                    var->addr = &constants[ 3 ];
                    break;
                }
                if ( var->name == CexmcCFVarConst_mm )
                {
                    var->addr = &constants[ 4 ];
                    break;
                }
                if ( var->name == CexmcCFVarConst_cm )
                {
                    var->addr = &constants[ 5 ];
                    break;
                }
                if ( var->name == CexmcCFVarConst_m )
                {
                    var->addr = &constants[ 6 ];
                    break;
                }
            } while ( false );
        }
    }
}


void  CexmcASTEval::ResetAddressBinding( CexmcAST::Subtree &  ast )
{
    for ( std::vector< CexmcAST::Node >::iterator  k( ast.children.begin() );
                                                  k != ast.children.end(); ++k )
    {
        CexmcAST::Subtree *  subtree( boost::get< CexmcAST::Subtree >( &*k ) );

        if ( subtree )
        {
            ResetAddressBinding( *subtree );
        }
        else
        {
            CexmcAST::Leaf &      leaf( boost::get< CexmcAST::Leaf >( *k ) );
            CexmcAST::Variable *  var( boost::get< CexmcAST::Variable >(
                                                                    &leaf ) );
            if ( var )
                var->addr = ( const int * ) NULL;
        }
    }
}

#endif

