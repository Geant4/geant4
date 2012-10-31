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
 *       Filename:  CexmcCustomFilterEval.cc
 *
 *    Description:  custom filter eval
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

#include <fstream>
#include <string>
#include "CexmcCustomFilterEval.hh"
#include "CexmcException.hh"


CexmcCustomFilterEval::CexmcCustomFilterEval( const G4String &  sourceFileName,
                                  const CexmcEventFastSObject *  evFastSObject,
                                  const CexmcEventSObject *  evSObject ) :
    astEval( evFastSObject, evSObject )
{
    std::string     command;
    std::ifstream   sourceFile( sourceFileName );

    if ( ! sourceFile )
        throw CexmcException( CexmcCFBadSource );

    bool            commandIsPending( false );

    while ( ! sourceFile.eof() )
    {
        std::string  line;
        std::getline( sourceFile, line );

        size_t  commentStartPos( line.find_first_of( '#' ) );
        if ( commentStartPos != std::string::npos )
            line.erase( commentStartPos );

        if ( line.empty() ||
             line.find_first_not_of( " \t" ) == std::string::npos )
        {
            if ( commandIsPending )
            {
                sourceFile.close();
                throw CexmcException( CexmcCFParseError );
            }
            continue;
        }

        command += line;

        size_t  length( command.length() );

        if ( command[ length - 1 ] == '\\' )
        {
            command.erase( length - 1 );
            commandIsPending = true;
            continue;
        }

        CexmcCustomFilter::ParseResult  curParseResult;

        std::string::const_iterator   begin( command.begin() );
        std::string::const_iterator   end( command.end() );

        try
        {
            if ( ! CexmcCustomFilter::phrase_parse( begin, end, grammar,
                                CexmcCustomFilter::space, curParseResult ) ||
                 begin != end )
            {
                throw CexmcException( CexmcCFParseError );
            }
        }
        catch ( ... )
        {
            sourceFile.close();
            throw;
        }

#ifdef CEXMC_DEBUG_CF
        G4cout << "Parsed expression AST:" << G4endl;
        curParseResult.expression.Print();
#endif

        switch ( curParseResult.action )
        {
        case CexmcCustomFilter::KeepTPT :
        case CexmcCustomFilter::DeleteTPT :
            parseResultTPT.push_back( curParseResult );
            break;
        case CexmcCustomFilter::KeepEDT :
        case CexmcCustomFilter::DeleteEDT :
            parseResultEDT.push_back( curParseResult );
            break;
        default :
            break;
        }

        command = "";
        commandIsPending = false;
    }

    sourceFile.close();

    if ( commandIsPending )
        throw CexmcException( CexmcCFParseError );
}


void  CexmcCustomFilterEval::SetAddressedData(
                                const CexmcEventFastSObject *  evFastSObject,
                                const CexmcEventSObject *  evSObject )
{
    astEval.SetAddressedData( evFastSObject, evSObject );

    for ( ParseResultVector::iterator  k( parseResultTPT.begin() );
          k != parseResultTPT.end(); ++k )
    {
        if ( evFastSObject == NULL || evSObject == NULL )
            astEval.ResetAddressBinding( k->expression );
        else
            astEval.BindAddresses( k->expression );
    }
    for ( ParseResultVector::iterator  k( parseResultEDT.begin() );
          k != parseResultEDT.end(); ++k )
    {
        if ( evFastSObject == NULL || evSObject == NULL )
            astEval.ResetAddressBinding( k->expression );
        else
            astEval.BindAddresses( k->expression );
    }
}


bool  CexmcCustomFilterEval::EvalTPT( void ) const
{
    for ( ParseResultVector::const_iterator  k( parseResultTPT.begin() );
          k != parseResultTPT.end(); ++k )
    {
        if ( astEval( k->expression ) )
            return k->action == CexmcCustomFilter::KeepTPT;
    }

    return true;
}


bool  CexmcCustomFilterEval::EvalEDT( void ) const
{
    for ( ParseResultVector::const_iterator  k( parseResultEDT.begin() );
          k != parseResultEDT.end(); ++k )
    {
        if ( astEval( k->expression ) )
            return k->action == CexmcCustomFilter::KeepEDT;
    }

    return true;
}


#endif

