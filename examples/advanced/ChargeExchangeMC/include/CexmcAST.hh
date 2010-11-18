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

    enum  Operator
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


    struct  OperatorPriority
    {
        static bool  IsLessThan( Operator  left, Operator  right )
        {
            return priority[ left ] < priority[ right ];
        }

        static const int  priority[];
            
    };


    struct  Variable
    {
        Variable() : index1 ( 0 ), index2( 0 ), addr( ( const int * ) NULL )
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
        Subtree() : type ( Uninitialized ), hasRLAssoc( false )
        {}

        void  Print( int  level = 0 ) const;

        void  PrintLeaf( const Leaf *  leaf, int  level = 0 ) const;

        std::vector< Node >  children;

        NodeType             type;

        bool                 hasRLAssoc;

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

