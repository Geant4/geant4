//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FluoTestRunAction.hh"
#include <stdlib.h>
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "FluoTestDataSet.hh"
#include "G4DataVector.hh"
#include "G4LogLogInterpolation.hh"
#include "g4std/fstream"
#include "g4std/strstream"
#include "FluoTestNormalization.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE
FluoTestRunAction::FluoTestRunAction(FluoTestAnalysisManager* aMgr)
  :analysisManager(aMgr)
{
 
}

FluoTestRunAction::FluoTestRunAction()
  :dataSet(0),dataGammaSet(0),dataAlphaSet(0),efficiencySet(0)
{
 
  /*
   G4double min = 0.*MeV;
   G4double max = 0.10000E+06 *MeV;
   G4int nBins = 835;
 */

  FluoTestNormalization* normalization = new FluoTestNormalization();
  
  energies = new G4DataVector;
  data = new G4DataVector;
  
  
    ReadData(MeV,  "mercury_flx_solmin");
    ReadResponse("response");
   
    /*  
    G4VDataSetAlgorithm* interpolation = new G4LogLogInterpolation();
    dataSet = 
    new FluoTestDataSet(1,"merc_flx_alp_min",interpolation,MeV,1);
  
    dataSet = normalization->Normalize(min, max, nBins,"mercury_flx_solmin");
    */
    G4double minGamma = 0.*keV;
    G4double maxGamma = 10. *keV;
    G4int nBinsGamma = 100;
  
  //G4VDataSetAlgorithm* interpolation2 = new G4LogLogInterpolation();
  // dataGammaSet = 
  //new FluoTestDataSet(1,"B_flare",interpolation2,MeV,1);
  dataGammaSet = normalization->Normalize(minGamma, maxGamma, nBinsGamma,"B_flare");
  /*
    G4double minAlpha = 0.*MeV;
    G4double maxAlpha = 100000. *MeV;
    G4int nBinsAlpha = 835;
    
    dataAlphaSet = normalization->Normalize(minAlpha, maxAlpha, nBinsAlpha,"merc_flx_alp_min");
    
    G4VDataSetAlgorithm* interpolation3 = new G4LogLogInterpolation();
 dataAlphaSet = 
 new FluoTestDataSet(1,"merc_flx_alp_min",interpolation3,MeV,1);
  */
  G4String fileName = "efficienza";
  G4VDataSetAlgorithm* interpolation4 = new G4LogLogInterpolation();
  efficiencySet = new FluoTestDataSet(1,fileName,interpolation4,keV,1);
  
  // delete normalization;
  
}
#else
FluoTestRunAction::FluoTestRunAction()
{
}
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestRunAction::~FluoTestRunAction()
{
    G4std::map<G4int,G4DataVector*,G4std::less<G4int> >::iterator pos;

  for (pos = energyMap.begin(); pos != energyMap.end(); pos++)
    {
      G4DataVector* dataSet = (*pos).second;
      delete dataSet;
    }
  for (pos = dataMap.begin(); pos != dataMap.end(); pos++)
    {
      G4DataVector* dataSet = (*pos).second;
      delete dataSet;
    }
  //delete responseFunction;
  //responseFunction = 0;
  // delete dataSet;
  //dataSet = 0;
  //delete energies;
  //delete data;
  //energies = 0;
  //data = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestRunAction::BeginOfRunAction(const G4Run* aRun)
{
  
G4cout << "### Run " << aRun << " start." << G4endl;

  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer(); 
      UI->ApplyCommand("/vis/scene/notifyHandlers");
    } 
#ifdef G4ANALYSIS_USE
  analysisManager->BeginOfRun();
#endif
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestRunAction::EndOfRunAction(const G4Run* aRun )
{
  // Run ended, update the visualization
if (G4VVisManager::GetConcreteInstance()) {
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }
 
 
// If analysis is used, print out the histograms
#ifdef G4ANALYSIS_USE
  analysisManager->EndOfRun(aRun->GetRunID());
#endif
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


const FluoTestDataSet* FluoTestRunAction::GetSet()
{
  return  dataSet;
}
const FluoTestDataSet* FluoTestRunAction::GetGammaSet()
{
  return  dataGammaSet;
}
const FluoTestDataSet* FluoTestRunAction::GetAlphaSet()
{
  return  dataAlphaSet;
}
const FluoTestDataSet* FluoTestRunAction::GetEfficiencySet()
{
  return efficiencySet;
}
G4DataVector* FluoTestRunAction::GetEnergies()
{
  return energies;
}
G4DataVector* FluoTestRunAction::GetData()
{
  return data;
}
G4double FluoTestRunAction::GetDataSum()
{
  G4double sum = 0;
  size_t size = data->size();
  for (G4int i = 0; i <size; i++)
    {
      sum+=(*data)[i];
    }
  return sum;
}
G4double FluoTestRunAction::GetInfData(G4double energy, G4double random)
{
  G4double value = 0.;
  G4int zMin = 1.;
  G4int zMax = 10.; 

  G4int Z = ((G4int)(energy/keV));
 
  if (Z<zMin) {Z=zMin;}
  if (Z>zMax) {Z=zMax;}
 
  if (Z >= zMin && Z <= zMax)
    {
      G4std::map<G4int,G4DataVector*,G4std::less<G4int> >::const_iterator pos;
      pos = energyMap.find(Z);
      G4std::map<G4int,G4DataVector*,G4std::less<G4int> >::const_iterator posData;
      posData = dataMap.find(Z);
      if (pos!= energyMap.end())
	{
	  G4DataVector energySet = *((*pos).second);
	  G4DataVector dataSet = *((*posData).second);
	  G4int nData = energySet.size();
	
	  G4double partSum = 0;
	  G4int index = 0;
	
	  while (random> partSum)
	    {
	      partSum += dataSet[index];
	      index++;
	    }
	  
	  
	  if (index >= 0 && index < nData)
	    {
	      value = energySet[index];
	  
	    }

	}
    }
  return value;
}

G4double FluoTestRunAction::GetSupData(G4double energy, G4double random)
{
  G4double value = 0.;
   G4int zMin = 1;
   G4int zMax = 10;
  G4int Z = ((G4int)(energy/keV)+1);

  if (Z<zMin) {Z=zMin;}
  if (Z>zMax) {Z=zMax;}
  if (Z >= zMin && Z <= zMax)
    {
      G4std::map<G4int,G4DataVector*,G4std::less<G4int> >::const_iterator pos;
      pos = energyMap.find(Z);
      G4std::map<G4int,G4DataVector*,G4std::less<G4int> >::const_iterator posData;
      posData = dataMap.find(Z);
      if (pos!= energyMap.end())
	{
	  G4DataVector energySet = *((*pos).second);
	  G4DataVector dataSet = *((*posData).second);
	  G4int nData = energySet.size();
	  G4double partSum = 0;
	  G4int index = 0;
	  
	  while (random> partSum)
	    {
	      partSum += dataSet[index];
	      index++;
	    }
	 
	  
	  if (index >= 0 && index < nData)
	    {
	      value = energySet[index];
	    
	      
	    }
	  
	}
    }
  return value;
}



void FluoTestRunAction::ReadData(G4double unitE, G4String fileName)
{
   char nameChar[100] = {""};
  G4std::ostrstream ost(nameChar, 100, G4std::ios::out);
  
  ost << fileName <<".dat";
  
  G4String name(nameChar);
  
  char* path = "/mnt/home/guardi/workdir/simpleFluo2/fluoTest";
 
  G4String pathString(path);
  G4String dirFile = pathString + "/" + name;
  G4std::ifstream file(dirFile);
  G4std::filebuf* lsdp = file.rdbuf();
  
  if (! (lsdp->is_open()) )
	{
	  G4String excep = "FluoTestRunAction - data file: " + dirFile + " not found";
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

void FluoTestRunAction::ReadResponse(const G4String& fileName)
{
   char nameChar[100] = {""};
  G4std::ostrstream ost(nameChar, 100, G4std::ios::out);
  
 
      ost << fileName<<".dat";
  
  G4String name(nameChar);
  
  char* path = "/mnt/home/guardi/workdir/simpleFluo2/fluoTest";
 
  G4String pathString(path);
  G4String dirFile = pathString + "/" + name;
  G4std::ifstream file(dirFile);
  G4std::filebuf* lsdp = file.rdbuf();
  
  if (! (lsdp->is_open()) )
	{
	  G4String excep = "FluoTestRunAction - data file: " + dirFile + " not found";
	  G4Exception(excep);
	}
  G4double a = 0;
  G4int k = 1;
  G4int s = 0;

  G4int Z = 1;
  G4DataVector* energies = new G4DataVector;
  G4DataVector* data = new G4DataVector;

  do
    {
      file >> a;
      G4int nColumns = 2;
      if (a == -1)
      {
	if (s == 0)
	  {
	    // End of a  data set
	    energyMap[Z] = energies;
            dataMap[Z] = data;
	    // Start of new shell data set
	    energies = new G4DataVector;
            data = new G4DataVector;
            Z++;	    
	  }      
	s++;
	if (s == nColumns)
	{
	  s = 0;
	}
      }
  else if (a == -2)
      {
	// End of file; delete the empty vectors created when encountering the last -1 -1 row
	delete energies;
	delete data;
	//nComponents = components.size();
      }
    else
      {
	// 1st column is energy
	if(k%nColumns != 0)
	  {	
	    G4double e = a * keV;
	    energies->push_back(e);
	    k++;
	  }
else if (k%nColumns == 0)
	  {
	    // 2nd column is data
	   
	    data->push_back(a);
	    k = 1;
	  }
      }
  } while (a != -2); // end of file
  file.close();    
}
