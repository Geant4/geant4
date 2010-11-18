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
        catch( const boost::bad_get & )
        {
            ast.type = Top;
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

        bool        hasLessPriorityThanRightOp( false );
        Operator *  rightOp( boost::get< Operator >( &astRight->type ) );

        if ( rightOp )
        {
            hasLessPriorityThanRightOp = OperatorPriority::IsLessThan(
                                                            value, *rightOp );
        }

        if ( astRight->hasRLAssoc || hasLessPriorityThanRightOp )
        {
            ast.children.push_back( right );
            return;
        }
        ast.children.push_back( astRight->children[ 0 ] );
        Subtree    astResult;
        astResult.children.push_back( self );
        astResult.children.push_back( astRight->children[ 1 ] );
        astResult.type = astRight->type;
        self = astResult;
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

        if ( ast )
            ast->hasRLAssoc = true;
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

