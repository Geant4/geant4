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
/// \file persistency/P03/src/ExTGRCLineProcessor.cc
/// \brief Implementation of the ExTGRCLineProcessor class

#include "ExTGRCLineProcessor.hh"
#include "ExTGRCRegionCutsMgr.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExTGRCLineProcessor::ExTGRCLineProcessor() : G4tgrLineProcessor()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExTGRCLineProcessor::~ExTGRCLineProcessor()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool ExTGRCLineProcessor::ProcessLine( const std::vector<G4String>& wl )
{

  G4bool iret = G4tgrLineProcessor::ProcessLine( wl );

  G4String wl0 = wl[0];
  for( size_t ii = 0; ii < wl0.length(); ii++ )
  {
    wl0[ii] = toupper( wl0[ii] );
  }

  if( !iret )
  {
    //------------------------------- parameter number
    if( wl0 == ":REGION" )
    {
      std::vector<G4String>::const_iterator ite = wl.begin()+1;
      std::vector<G4String> wlc;
      for( ; ite != wl.end(); ite++ )   //loop skipping the first one
      {
        wlc.push_back( *ite );
      }
      //      wlc = wlc.erase( wlc.begin() );
      ExTGRCRegionCutsMgr::GetInstance()->AddRegionData( wlc );
      iret = 1; 

    }
    else if( wl0 == ":CUTS" )
    {
      std::vector<G4String>::const_iterator ite = wl.begin()+1;
      std::vector<G4String> wlc;
      for( ; ite != wl.end(); ite++ )   //loop skipping the first one
      {
        wlc.push_back( *ite );
      }
      ExTGRCRegionCutsMgr::GetInstance()->AddRegionCuts( wlc );
      iret = 1; 
    }
    else
    {
      iret = 0; 
    } 
  }

  return iret;
}
