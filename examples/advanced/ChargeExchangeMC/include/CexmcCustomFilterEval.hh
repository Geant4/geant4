/*
 * =============================================================================
 *
 *       Filename:  CexmcCustomFilterEval.hh
 *
 *    Description:  custom filter eval
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

#ifndef CEXMC_CUSTOM_FILTER_EVAL_HH
#define CEXMC_CUSTOM_FILTER_EVAL_HH

#ifdef CEXMC_USE_CUSTOM_FILTER

#include <vector>
#include <string>
#include "CexmcASTEval.hh"
#include "CexmcCustomFilter.hh"

class  CexmcEventFastSObject;
class  CexmcEventSObject;


class  CexmcCustomFilterEval : public CexmcAST::BasicEval
{
    private:
        typedef std::vector< CexmcCustomFilter::ParseResult >
                                                        ParseResultVector;

    public:
        explicit CexmcCustomFilterEval( const G4String &  sourceFileName,
                            const CexmcEventFastSObject *  evFastSObject = NULL,
                            const CexmcEventSObject *  evSObject = NULL );

    public:
        void  SetAddressedData( const CexmcEventFastSObject *  evFastSObject,
                                const CexmcEventSObject *  evSObject );

        bool  EvalTPT( void ) const;

        bool  EvalEDT( void ) const;

    private:
        CexmcASTEval       astEval;

        ParseResultVector  parseResultTPT;

        ParseResultVector  parseResultEDT;

        CexmcCustomFilter::Grammar< std::string::const_iterator >  grammar;
};

#endif

#endif

