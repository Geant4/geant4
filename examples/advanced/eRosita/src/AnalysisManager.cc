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
  outFile.open("hits.out");
  outFileT.open("eTot.out"); 
}


AnalysisManager::~AnalysisManager() 
{
  outFile.close();
  outFileT.close();
}


void AnalysisManager::Score(G4int id, G4double eDep, G4double x, G4double y, G4double z) 
{
  outFile << id << " "
	  << eDep << " "
	  << x << " "
	  << y << " "
	  << z << " "
	  << std::endl;
}

void AnalysisManager::ScoreTot(G4double eTot) 
{
  outFileT << eTot << std::endl;
}
