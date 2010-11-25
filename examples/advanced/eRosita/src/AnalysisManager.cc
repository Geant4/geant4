#include "AnalysisManager.hh"
#include <sstream>

AnalysisManager* AnalysisManager::instance = 0;


AnalysisManager* AnalysisManager::Instance() 
{
 
  // A new instance of AnalysisManager is created if it does not exist:
  if (instance == 0) 
    {
      instance = new AnalysisManager();
    }
  
  // The instance of AnalysisManager is returned:
  return instance;
}


void AnalysisManager::Destroy() 
{

  // The AnalysisManager instance is deleted if it exists:
  if (instance != 0) 
    {
      delete instance;
      instance = 0;
    }
}


AnalysisManager::AnalysisManager() 
{
  outFile.open("TrackerPhotonEnergy.out");
//   outFileT.open("eTot.out"); 
}


AnalysisManager::~AnalysisManager() 
{
  outFile.close();
  outFileT.close();
}


void AnalysisManager::Score(G4double eDep) 
{
//   outFile << id << " "
// 	  << eDep << " "
// 	  << x << " "
// 	  << y << " "
// 	  << z << " "
// 	  << std::endl;
  outFile << eDep << std::endl;
}

// void AnalysisManager::ScoreTot(G4double eTot) 
// {
//   outFileT << eTot << std::endl;
// }
