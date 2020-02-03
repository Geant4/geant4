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
//
/// \file persistency/P03/src/ExTGRCRegionData.cc
/// \brief Implementation of the ExTGRCRegionData class

#include "ExTGRCRegionData.hh"
#include "G4tgrUtils.hh"
#include "G4UIcommand.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExTGRCRegionData::ExTGRCRegionData(const std::vector<G4String>& wl )
{
  if( wl.size() < 2 )
  { 
    G4Exception("ExTGRCRegionData::ExTGRCRegionData()",
                "InvalidArgument", FatalErrorInArgument,
                G4UIcommand::ConvertToString( G4int(wl.size()) ) );
  }
  fRegionName = wl[0];
  for( size_t ii = 1; ii < wl.size(); ii++ )
  {
    fLVNames.push_back( wl[ii] );
  }
  fbCutsSet = false;
  fGammaCut = 1.;
  fElectronCut = 1.;
  fPositronCut = 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExTGRCRegionData::~ExTGRCRegionData()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
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

  if( fbCutsSet )
  {
    G4Exception("ExTGRCRegionData::SetCutsData()",
                "InvalidArgument", JustWarning,
    G4String("Cuts are already set for region " + fRegionName).c_str() );
  }

  fGammaCut = G4tgrUtils::GetDouble( rc[1] );
  fElectronCut = G4tgrUtils::GetDouble( rc[2] );
  if( rc.size() == 3 )
  {
    fPositronCut = fElectronCut;
  }
  else
  {
    fPositronCut = G4tgrUtils::GetDouble( rc[3] );
  }

  fbCutsSet = true;
}
