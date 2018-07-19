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
 *       Filename:  CexmcAST.hh
 *
 *    Description:  abstract syntax tree for custom filter scripting language
 *
 *        Version:  1.0
 *        Created:  17.07.2010 14:39:37
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_AST_HH
#define CEXMC_AST_HH

#ifdef CEXMC_USE_CUSTOM_FILTER

#include <vector>
#include <boost/variant/recursive_variant.hpp>


namespace  CexmcAST
{
    using boost::variant;
    using boost::recursive_wrapper;

    enum  OperatorType
    {
        Uninitialized,
        Top,
        UMinus,
        Not,
        Mult,
        Div,
        Plus,
        Minus,
        Less,
        LessEq,
        More,
        MoreEq,
        Eq,
        NotEq,
        And,
        Or
    };


    struct  Operator
    {
        Operator( OperatorType  type_ = Uninitialized, int  priority_ = 0,
                  bool  hasRLAssoc_ = false ) :
            type( type_ ), priority( priority_ ), hasRLAssoc( hasRLAssoc_ )
        {}

        OperatorType  type;

        int           priority;

        bool          hasRLAssoc;
    };


    struct  Variable
    {
        Variable() : index1( 0 ), index2( 0 ), addr( ( const int * ) NULL )
        {}

        std::string                             name;

        int                                     index1;

        int                                     index2;

        variant< const int *, const double * >  addr;
    };


    struct  Subtree;

    typedef std::string                    Function;

    typedef variant< int, double >         Constant;

    typedef variant< Variable, Constant >  Leaf;

    typedef recursive_wrapper< Subtree >   Tree;

    typedef variant< Tree, Leaf >          Node;

    typedef variant< Operator, Function >  NodeType;


    struct  Subtree
    {
        Subtree() : type( Operator( Uninitialized ) )
        {}

        void  Print( int  level = 0 ) const;

        void  PrintLeaf( const Leaf *  leaf, int  level = 0 ) const;

        std::vector< Node >  children;

        NodeType             type;

        static const int     printIndent = 4;
    };


    class  BasicEval
    {
        protected:
            typedef variant< int, double >  ScalarValueType;

        protected:
            virtual ~BasicEval();

        public:
            bool  operator()( const Subtree &  ast ) const;

        protected:
            ScalarValueType          GetScalarValue( const Node &  node ) const;

            virtual ScalarValueType  GetFunScalarValue( const Subtree &  ast )
                                                                        const;

            virtual ScalarValueType  GetVarScalarValue( const Variable &  var )
                                                                        const;

            ScalarValueType          GetBasicFunScalarValue(
                                        const Subtree &  ast, bool &  result )
                                                                        const;
    };
}

#endif

#endif

