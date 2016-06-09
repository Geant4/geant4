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
// $Id: ExTGRCRegionData.cc,v 1.4 2010-11-05 08:52:34 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author:      P. Arce
// Changes:     creation   May 2007
// ---------------------------------------------------------------------------

#include "ExTGRCRegionData.hh"
#include "G4tgrUtils.hh"
#include "G4UIcommand.hh"

//-----------------------------------------------------------------------
ExTGRCRegionData::ExTGRCRegionData(const std::vector<G4String>& wl )
{
  if( wl.size() < 2 )
  { 
    G4Exception("ExTGRCRegionData::ExTGRCRegionData()",
                "InvalidArgument", FatalErrorInArgument,
                G4UIcommand::ConvertToString( G4int(wl.size()) ) );
  }
  theRegionName = wl[0];
  for( size_t ii = 1; ii < wl.size(); ii++ )
  {
    theLVNames.push_back( wl[ii] );
  }
  bCutsSet = false;
}

//-----------------------------------------------------------------------
ExTGRCRegionData::~ExTGRCRegionData()
{
}

//-----------------------------------------------------------------------
void ExTGRCRegionData::SetCutsData( const std::vector<G4String>& rc )
{
  if( (rc.size() != 3) && (rc.size() != 4) )
  { 
    G4cerr << "ERROR - ExTGRCRegionData::SetCutsData()" << G4endl
           << "        Must have 3 or 4 arguments : "
           << "REGION_NAME, gamma_CUT, e-_CUT (e+_CUT)." << G4endl
           << "        It has only " << rc.size() << " !" << G4endl; 
    G4Exception("ExTGRCRegionCutsMgr::AddRegionCuts()",
                "InvalidArgument", FatalErrorInArgument,
                G4UIcommand::ConvertToString( G4int(rc.size()) ) );
  }

  if( bCutsSet )
  {
    G4Exception("ExTGRCRegionData::SetCutsData()",
                "InvalidArgument", JustWarning,
    G4String("Cuts are already set for region " + theRegionName).c_str() );
  }

  theGammaCut = G4tgrUtils::GetDouble( rc[1] );
  theElectronCut = G4tgrUtils::GetDouble( rc[2] );
  if( rc.size() == 3 )
  {
    thePositronCut = theElectronCut;
  }
  else
  {
    thePositronCut = G4tgrUtils::GetDouble( rc[3] );
  }

  bCutsSet = true;
}
