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
/// \file persistency/P03/src/ExTGRCRegionCutsMgr.cc
/// \brief Implementation of the ExTGRCRegionCutsMgr class

#include "ExTGRCRegionCutsMgr.hh"
#include "ExTGRCRegionData.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include "G4LogicalVolume.hh"
#include "G4tgbVolumeMgr.hh"
#include "G4tgrUtils.hh"
#include "G4UIcommand.hh"

ExTGRCRegionCutsMgr* ExTGRCRegionCutsMgr::fInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExTGRCRegionCutsMgr* ExTGRCRegionCutsMgr::GetInstance()
{
  if( !fInstance )
  {
    fInstance = new ExTGRCRegionCutsMgr;
  }
  return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExTGRCRegionCutsMgr::ExTGRCRegionCutsMgr()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExTGRCRegionCutsMgr::~ExTGRCRegionCutsMgr()
{
  delete fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExTGRCRegionCutsMgr::AddRegionData( const std::vector<G4String>& rd )
{

  if( (rd.size() > 1) && (FindRegionData( rd[0] ).size() != 0) )
  {
    G4Exception("ExTGRCRegionCutsMgr::AddRegionData", "InvalidArgument",
                JustWarning,
                G4String("Region already exists: " + rd[0]).c_str() );
    return;
  }
  fRegionDatae.push_back( new ExTGRCRegionData( rd ) );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExTGRCRegionCutsMgr::AddRegionCuts( const std::vector<G4String>& rc )
{
  if( rc.size() == 0 )
  { 
    G4cerr << "ERROR - ExTGRCRegionCutsMgr::AddRegionCuts()" << G4endl
           << "        Must have 3 or 4 arguments : REGION_NAME, gamma_CUT,"
           << " e-_CUT (e+_CUT)." << G4endl
           << "        It has only " << rc.size() << " !" << G4endl; 
    G4Exception("ExTGRCRegionCutsMgr::AddRegionCuts()", "InvalidArgument",
                FatalErrorInArgument, G4UIcommand::ConvertToString(G4int(rc.size())) );
  }

  // Find region
  // std::vector<ExTGRCRegionData*>::const_iterator iter;
  std::vector<ExTGRCRegionData*> regs = FindRegionData(rc[0]);

  if( regs.size() == 0 )
  {
    G4Exception("ExTGRCRegionCutsMgr::AddRegionCuts()",
                "InvalidArgument", FatalErrorInArgument,
                G4String(" region does not exist: " + rc[0]).c_str());
  } 

  for( size_t ii = 0; ii < regs.size(); ii++)
  {
    regs[ii]->SetCutsData( rc );
  }
}

std::vector<ExTGRCRegionData*>
ExTGRCRegionCutsMgr::FindRegionData( const G4String& name)
{
  std::vector<ExTGRCRegionData*> regs;
  std::vector<ExTGRCRegionData*>::const_iterator iter;
  for( iter = fRegionDatae.begin(); iter != fRegionDatae.end(); iter++ )
  {
    if( G4tgrUtils::AreWordsEquivalent( name , (*iter)->GetRegionName()) )
    {
      regs.push_back(*iter);
    }
  }
  return regs; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExTGRCRegionCutsMgr::BuildRegions()
{
  std::vector<ExTGRCRegionData*>::const_iterator iter;
  std::vector<G4String>::const_iterator ites;
  // std::vector<G4LogicalVolume*>::const_iterator itelv;
  for( iter = fRegionDatae.begin(); iter != fRegionDatae.end(); iter++ )
  {
    G4Region* reg = new G4Region( (*iter)->GetRegionName() );
    std::vector<G4String> lvs = (*iter)->GetLVNames();
    for( ites = lvs.begin(); ites != lvs.end(); ites++ )
    {
      G4LogicalVolume* logVol =
        G4tgbVolumeMgr::GetInstance()->FindG4LogVol(*ites, true );
      reg->AddRootLogicalVolume( logVol );
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExTGRCRegionCutsMgr::BuildProductionCuts()
{
  std::vector<ExTGRCRegionData*>::const_iterator iter;
  G4RegionStore* regions = G4RegionStore::GetInstance();
  //----- loop to region datae
  for( iter = fRegionDatae.begin(); iter != fRegionDatae.end(); iter++ )
  {
    if( (*iter)->CutsAreSet() )
    {
      G4Region* reg = regions->GetRegion( (*iter)->GetRegionName() );
      if( !reg )
      { 
        G4Exception("ExTGRCRegionCutsMgr::BuildProductionCuts()",
                    "InvalidArgument", FatalErrorInArgument,
        G4String("Region not found: " + (*iter)->GetRegionName()).c_str() );
      }        
      G4ProductionCuts* cuts = new G4ProductionCuts ;

      cuts->SetProductionCut((*iter)->GetGammaCut(),"gamma");
      cuts->SetProductionCut((*iter)->GetElectronCut(),"e-");
      cuts->SetProductionCut((*iter)->GetPositronCut(),"e+");
      reg->SetProductionCuts(cuts);
    }
  }
}
