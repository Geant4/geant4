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
 *       Filename:  CexmcSimpleRangeWithValue.cc
 *
 *    Description:  auxiliary functions for simple range instances
 *
 *        Version:  1.0
 *        Created:  17.02.2010 22:46:01
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#include <iostream>
#include <iomanip>
#include <G4UnitsTable.hh>
#include "CexmcSimpleRangeWithValue.hh"


std::ostream &  operator<<( std::ostream &  out,
                        const CexmcEnergyRangeWithDoubleValue &  range )
{
    out << "[ " << std::setw( 3 ) << G4BestUnit( range.bottom, "Energy" ) <<
           ", " << std::setw( 3 ) << G4BestUnit( range.top, "Energy" ) <<
           " )   " << range.value;

    return out;
}


std::ostream &  operator<<( std::ostream &  out,
                        const CexmcEnergyRangeWithDoubleValueList &  ranges )
{
    out << std::endl;
    for ( CexmcEnergyRangeWithDoubleValueList::const_iterator
                                  k( ranges.begin() ); k != ranges.end(); ++k )
    {
        out << "          " << *k << std::endl;
    }

    return out;
}

