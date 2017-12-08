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
 *       Filename:  CexmcCustomFilter.cc
 *
 *    Description:  custom filter grammar and compiler
 *
 *        Version:  1.0
 *        Created:  17.07.2010 15:37:38
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifdef CEXMC_USE_CUSTOM_FILTER

#include "CexmcCustomFilter.hh"


namespace  CexmcCustomFilter
{
    void  Compiler::operator()( ParseResult &  parseResult, Action  value )
                                                                        const
    {
        parseResult.action = value;
    }


    void  Compiler::operator()( ParseResult &  parseResult, Subtree &  value )
                                                                        const
    {
        parseResult.expression = value;
    }


    void  Compiler::operator()( Subtree &  ast, Node &  node ) const
    {
        try
        {
            ast = boost::get< Subtree >( node );
        }
        catch ( const boost::bad_get & )
        {
            ast.type = Operator( Top );
            ast.children.push_back( node );
        }
    }


    void  Compiler::operator()( Node &  self, Node &  left, Node &  right,
                                Operator  value ) const
    {
        Subtree &  ast( boost::get< Subtree >( self ) );

        ast.children.push_back( left );
        ast.type = value;

        Subtree *  astRight( boost::get< Subtree >( &right ) );

        if ( ! astRight )
        {
            ast.children.push_back( right );
            return;
        }

        bool        haveSamePriorities( false );
        Operator *  rightOp( boost::get< Operator >( &astRight->type ) );

        if ( rightOp )
            haveSamePriorities = value.priority == rightOp->priority;

        if ( value.hasRLAssoc || ! haveSamePriorities )
        {
            ast.children.push_back( right );
            return;
        }

        Subtree *  astDeepestRight( astRight );

        /* propagate left binary operators with LR associativity (i.e. all in
         * our grammar) deep into the AST until any operator with a different
         * priority (which includes operators in parentheses that have priority
         * 0) or a unary operator or a function occured */
        while ( true )
        {
            Subtree *   candidate = boost::get< Subtree >(
                                            &astDeepestRight->children[ 0 ] );
            if ( ! candidate )
                break;

            if ( candidate->children.size() < 2 )
                break;

            bool        haveSamePriorities_( false );
            Operator *  candidateOp( boost::get< Operator >(
                                                        &candidate->type ) );

            if ( candidateOp )
                haveSamePriorities_ = value.priority == candidateOp->priority;

            /* FIXME: what to do if candidate has RL association? Our grammar is
             * not a subject of this issue; probably no grammar is a subject */
            if ( ! haveSamePriorities_ )
                break;

            astDeepestRight = candidate;
        }

        Subtree    astResult;
        astResult.children.push_back( ast.children[ 0 ] );
        astResult.children.push_back( astDeepestRight->children[ 0 ] );
        astResult.type = value;
        astDeepestRight->children[ 0 ] = astResult;
        self = right;
    }


    void  Compiler::operator()( Node &  self, Node &  child, Operator  value )
                                                                        const
    {
        Subtree &  ast( boost::get< Subtree >( self ) );
        ast.children.push_back( child );
        ast.type = value;
    }


    void  Compiler::operator()( Node &  self, Node &  primary ) const
    {
        self = primary;

        Subtree *  ast( boost::get< Subtree >( &self ) );

        if ( ! ast )
            return;

        Operator *  op( boost::get< Operator >( &ast->type ) );

        if ( op )
            op->priority = 0;
    }


    void  Compiler::operator()( Node &  self, Node &  child,
                                std::string &  value ) const
    {
        Subtree &  ast( boost::get< Subtree >( self ) );

        ast.children.push_back( child );
        ast.type = value;
    }


    void  Compiler::operator()( Leaf &  self, std::string &  name ) const
    {
        Variable &  variable( boost::get< Variable >( self ) );
        variable.name = name;
    }


    void  Compiler::operator()( Leaf &  self, int  value, size_t  index ) const
    {
        Variable &  variable( boost::get< Variable >( self ) );
        switch ( index )
        {
        case 0 :
            variable.index1 = value;
            break;
        case 1 :
            variable.index2 = value;
            break;
        default :
            break;
        }
    }
}

#endif

