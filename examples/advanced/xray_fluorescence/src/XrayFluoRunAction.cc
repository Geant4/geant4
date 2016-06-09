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
// $Id: XrayFluoRunAction.cc
// GEANT4 tag $Name: xray_fluo-V03-02-00
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "XrayFluoRunAction.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "XrayFluoDataSet.hh"
#include "G4DataVector.hh"
#include "G4LogLogInterpolation.hh"
#include <fstream>
#include <sstream>
#include "XrayFluoNormalization.hh"
#ifdef G4ANALYSIS_USE
#include "XrayFluoAnalysisManager.hh"
#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE

XrayFluoRunAction::XrayFluoRunAction()
  :dataSet(0),dataGammaSet(0),dataAlphaSet(0)
{
  XrayFluoNormalization normalization;
  
  energies = new G4DataVector;
  data = new G4DataVector;
  
  
  ReadData(keV,"spec10");
  //ReadResponse("SILIresponse");
  
  G4double minGamma = 0.*keV;
  G4double maxGamma = 10. *keV;
  G4int nBinsGamma = 6;
  

  dataGammaSet = normalization.Normalize(minGamma, maxGamma, nBinsGamma,
				  "M_flare");
  

  //G4String fileName = "SILIefficiency";
  //G4VDataSetAlgorithm* interpolation4 = new G4LogLogInterpolation();
  //efficiencySet = new XrayFluoDataSet(1,fileName,interpolation4,keV,1);
  //delete interpolation4;  
  G4cout << "XrayFluoRunAction created" << G4endl;  
}
#else
XrayFluoRunAction::XrayFluoRunAction()
{
  G4cout << "XrayFluoRunAction created" << G4endl; 
}
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
#ifdef G4ANALYSIS_USE

XrayFluoRunAction::~XrayFluoRunAction()
{

  //std::map<G4int,G4DataVector*,std::less<G4int> >::iterator pos;
  
  // delete energies;
  // delete data;
  // G4cout << "energies and data deleted " << G4endl;

  //for (pos = energyMap.begin(); pos != energyMap.end(); pos++)
  //{
  //  G4DataVector* dataSet = (*pos).second;
  //  delete dataSet;
  //  dataSet = 0;
  //  }
  //for (pos = dataMap.begin(); pos != dataMap.end(); pos++)
  //  {
  //    G4DataVector* dataSet = (*pos).second;
  //    delete dataSet;
  //    dataSet = 0;
  //  }
  

  G4cout << "XrayFluoRunAction deleted" << G4endl; 

}

#else
XrayFluoRunAction::~XrayFluoRunAction()
{
  G4cout << "XrayFluoRunAction deleted" << G4endl;   
}
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoRunAction::BeginOfRunAction(const G4Run* aRun)
{
  
  G4cout << "### Run " << aRun << " start." << G4endl;
  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer(); 
      UI->ApplyCommand("/vis/scene/notifyHandlers");
    } 
#ifdef G4ANALYSIS_USE

  // Book histograms and ntuples
  XrayFluoAnalysisManager* analysis = XrayFluoAnalysisManager::getInstance();
  analysis->book();
  analysis->InitializePlotter();
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoRunAction::EndOfRunAction(const G4Run*)
{
#ifdef G4ANALYSIS_USE
  XrayFluoAnalysisManager* analysis = XrayFluoAnalysisManager::getInstance();
#endif
  // Run ended, update the visualization
  if (G4VVisManager::GetConcreteInstance()) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }
#ifdef G4ANALYSIS_USE
   analysis->finish();
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


const XrayFluoDataSet* XrayFluoRunAction::GetSet()
{
  return  dataSet;
}
const XrayFluoDataSet* XrayFluoRunAction::GetGammaSet()
{
  return  dataGammaSet;
}
const XrayFluoDataSet* XrayFluoRunAction::GetAlphaSet()
{
  return  dataAlphaSet;
}
//const XrayFluoDataSet* XrayFluoRunAction::GetEfficiencySet()
//{
//  return efficiencySet;
//}
G4DataVector* XrayFluoRunAction::GetEnergies()
{
  return energies;
}
G4DataVector* XrayFluoRunAction::GetData()
{
  return data;
}
G4double XrayFluoRunAction::GetDataSum()
{
  G4double sum = 0;
  size_t size = data->size();
  for (size_t i = 0; i <size; i++)
    {
      sum+=(*data)[i];
    }
  return sum;
}


void XrayFluoRunAction::ReadData(G4double unitE, G4String fileName)
{
  std::ostringstream ost;
  
  ost << fileName <<".dat";
  
  G4String name = ost.str();
  char* path;
  
  if (!(getenv("XRAYDATA"))) { 
    
    path = getenv("PWD");    
  }
  
  else {    
    path = getenv("XRAYDATA");
  }
  
  
  G4String pathString(path);
  name = pathString + "/" + name;
  
  
  std::ifstream file(name);
  std::filebuf* lsdp = file.rdbuf();
  
  if (! (lsdp->is_open()) )
    {
      G4String excep = "XrayFluoRunAction - data file: " + name + " not found";
      G4Exception(excep);
    }
  G4double a = 0;
  G4int k = 1;
  
  do
    {
      file >> a;
      G4int nColumns = 2;
      // The file is organized into two columns:
      // 1st column is the energy
      // 2nd column is the corresponding value
      // The file terminates with the pattern: -1   -1
      //                                       -2   -2
      if (a == -1 || a == -2)
	{
	  
	}
      else
	{
	  if (k%nColumns != 0)
	    {	
	      G4double e = a * unitE;
	      
	      energies->push_back(e);
	      
	      k++;
	      
	    }
	  else if (k%nColumns == 0)
	    {
	      G4double value = a;
	      data->push_back(value);
	      
	      k = 1;
	    }
	}
      
    } while (a != -2); // end of file
  
  file.close();
}














