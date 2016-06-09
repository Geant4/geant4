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
 *       Filename:  CexmcPrimitiveScorer.cc
 *
 *    Description:  primitive scorer with a messenger
 *
 *        Version:  1.0
 *        Created:  24.01.2011 18:34:30
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#include <G4ios.hh>
#include "CexmcPrimitiveScorer.hh"
#include "CexmcSensitiveDetectorMessenger.hh"


CexmcPrimitiveScorer::CexmcPrimitiveScorer( const G4String &  name ) :
    G4VPrimitiveScorer( name ), messenger( NULL )
{
}


CexmcPrimitiveScorer::~CexmcPrimitiveScorer()
{
    delete messenger;
}


void  CexmcPrimitiveScorer::InitializeMessenger( void )
{
    messenger = new CexmcSensitiveDetectorMessenger( this );
}


void  CexmcPrimitiveScorer::PrintHeader( G4int  nmbOfEntries ) const
{
    G4cout << " --- MultiFunctionalDet " << detector->GetName() << G4endl;
    G4cout << "     PrimitiveScorer " << primitiveName << G4endl;
    G4cout << "     Number of entries " << nmbOfEntries << G4endl;
}

