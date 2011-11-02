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


class  CexmcCustomFilterEval
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

