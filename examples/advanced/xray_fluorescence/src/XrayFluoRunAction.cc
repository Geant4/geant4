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
#include "XrayFluoAnalysisManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoRunAction::XrayFluoRunAction()
  : isInitialized(false), dataSet(0), dataGammaSet(0), 
    dataAlphaSet(0)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoRunAction::~XrayFluoRunAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoRunAction::Initialise()
{
  //Only the master is initialized and keeps the data
  //(or in sequential mode)
  if (!IsMaster())
    return;

  XrayFluoNormalization normalization;
  
  energies = new G4DataVector;
  data = new G4DataVector;
  
  
  ReadData(keV,"M_flare");
  //ReadResponse("SILIresponse");
  
  G4double minGamma = 0.*keV;
  G4double maxGamma = 10. *keV;
  G4int nBinsGamma = 6;
  

  dataGammaSet = normalization.Normalize(minGamma, maxGamma, nBinsGamma,
				  "M_flare");
  isInitialized = true;
  G4cout << "XrayFluoRunAction initialized" << G4endl;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoRunAction::BeginOfRunAction(const G4Run* aRun)
{
 
  //Master mode or sequential
  if (IsMaster())
    {
      G4cout << "### Run " << aRun->GetRunID() << " starts (master)." << G4endl;
      if (!isInitialized)
	Initialise();
    }
  else    
    G4cout << "### Run " << aRun->GetRunID() << " starts (worker)." << G4endl;

  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer(); 
      UI->ApplyCommand("/vis/scene/notifyHandlers");
    } 

  // Book histograms and ntuples
  XrayFluoAnalysisManager* analysis = XrayFluoAnalysisManager::getInstance();
  analysis->book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoRunAction::EndOfRunAction(const G4Run*)
{
  XrayFluoAnalysisManager* analysis = XrayFluoAnalysisManager::getInstance();
  analysis->finish();
  // Run ended, update the visualization
  if (G4VVisManager::GetConcreteInstance()) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


const XrayFluoDataSet* XrayFluoRunAction::GetSet() const
{
  return  dataSet;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const XrayFluoDataSet* XrayFluoRunAction::GetGammaSet() const
{
  return  dataGammaSet;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const XrayFluoDataSet* XrayFluoRunAction::GetAlphaSet() const
{
  return  dataAlphaSet;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DataVector* XrayFluoRunAction::GetEnergies() const
{
  return energies;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DataVector* XrayFluoRunAction::GetData() const
{
  return data;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XrayFluoRunAction::GetDataSum() const
{
 
  G4double sum = 0;
  for (size_t i = 0; i < data->size(); i++)
    {
      sum+=(*data)[i];
    }
  return sum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoRunAction::ReadData(G4double unitE, G4String fileName)
{
  G4cout << "Reading data...";
  std::ostringstream ost;
  
  ost << fileName <<".dat";
  
  G4String name = ost.str();
  char* path;
  
  if (!(std::getenv("XRAYDATA"))) { 
    
    path = std::getenv("PWD");    
  }
  
  else {    
    path = std::getenv("XRAYDATA");
  }
  
  
  G4String pathString(path);
  name = pathString + "/" + name;
  
  
  std::ifstream file(name);
  std::filebuf* lsdp = file.rdbuf();
  
  if (! (lsdp->is_open()) )
    {
      G4ExceptionDescription execp;
      execp <<  "XrayFluoRunAction - data file: " + name + " not found";
      G4Exception("XrayFluoRunAction::ReadData()","example-xray_fluorescence04",
	  FatalException, execp);
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
  G4cout << " done" << G4endl;
}














