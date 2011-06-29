
#include "DataManager.hh"
#include "StorageManager.hh"
#include <string>
#include <iostream>
#include <fstream>
#include <AIDA/AIDA.h>


DataManager::DataManager(double zLow, double zUp) : 
    zAxisLowerBound(zLow),
    zAxisUpperBound(zUp),
    eps(0.0000000001) {

  std::string fileNameBase = "Sandia";
  storageManager = StorageManager::Instance(fileNameBase);  

  histogramFactory = storageManager -> GetIHistogramFactory(); 
}


DataManager::~DataManager() {

  collector::iterator iter = dataCollectors.begin();
  collector::iterator iter_end = dataCollectors.end();
 
  for(;iter != iter_end; iter++) {
      delete iter -> second;
  }

  iter = garbage.begin();
  iter_end = garbage.end();
 
  for(;iter != iter_end; iter++) {
      delete iter -> second;
  }

  storageManager -> Destroy();
}


void DataManager::ScoreEnergyDeposit(double en, 
                                     double x, 
                                     double y, 
                                     double z,
                                     std::string type) {

  if(z < zAxisLowerBound || z > zAxisUpperBound) {
     return;
  }

  if(dataCollectors.begin() == dataCollectors.end()) return;

  collector::iterator iter = std::upper_bound(dataCollectors.begin(),
                                              dataCollectors.end(),
                                              z,zCompare());

  if(iter != dataCollectors.begin()) 
       (iter-1) -> second -> ScoreEnergyDeposit(en,x,y,z,type);
  
} 


void DataManager::ScoreParticleEnergy(double en, 
                                      double x, 
                                      double y, 
                                      double z,
                                      std::string type) {

  if(z < zAxisLowerBound || z >= zAxisUpperBound) return;

  collector::iterator iter = std::lower_bound(dataCollectors.begin(),
                                              dataCollectors.end(),
                                              z,zCompare());
  
  if(iter != dataCollectors.end()) {
     (iter-1) -> second -> ScoreParticleEnergy(en,x,y,z,type);
  }
} 


void DataManager::PrintResults() {

  collector::iterator iter = dataCollectors.begin();
  collector::iterator iter_end = dataCollectors.end();

  for(;iter != iter_end; iter++) {
     iter -> second -> PrintResults();
  }

}


void DataManager::AddDataCollector(DataManager* comp) {
 
  double zLow = comp -> GetLowerBound();
  double zUp  = comp -> GetUpperBound();

  if(!IsContained(zLow,zUp)) {
     std::cerr << "Error. Slab with lower bound z="     << zLow 
               << " exceeds boundary of mother volume." << std::endl;
     garbage.push_back(std::make_pair(zLow,comp));
     return;
  }

  if(!OverlapsWithOtherChild(zLow,zUp)) {
     dataCollectors.push_back(std::make_pair(zLow,comp));

     double thickn = zUp - zLow;
     double zPos = zLow + 0.5 * thickn;
     std::cerr << "INFORMATION. Slab with center at z="  << zPos 
               << " and thickness " << thickn <<" added." << std::endl;
    
  }
  else {
     std::cerr << "Error. Slab with lower bound z="  << zLow 
               << " overlaps with other slab." << std::endl;
     garbage.push_back(std::make_pair(zLow,comp));
     return;
  }

  std::sort(dataCollectors.begin(),dataCollectors.end(),zCompare());
}


bool DataManager::HasDataCollector(double z) {

  if(std::binary_search(dataCollectors.begin(),
                        dataCollectors.end(),z,zCompare())) {
     return true;
  }
 
  return false;
}


bool DataManager::OverlapsWithOtherChild(double zLow, double zUp) {

  collector::iterator iter = dataCollectors.begin();
  collector::iterator iter_end = dataCollectors.end();

  for(;iter != iter_end; iter++) {

     if(zLow < iter -> second -> GetLowerBound() && 
        zUp > (iter -> second -> GetLowerBound() + eps)) return true;

     if(zLow < (iter -> second -> GetUpperBound() - eps) && 
        zUp > iter -> second -> GetUpperBound()) return true;

     if(iter -> second -> IsContained(zLow,zUp)) return true;
  }

  return false;
}


bool DataManager::IsContained(double zLow, double zUp){

  if(zLow >= zAxisLowerBound && zUp <= zAxisUpperBound) return true;  

  return false;
}


DataManager* DataManager::MatchingChildDataCollector(double zLow, double zUp) {

  collector::iterator iter = dataCollectors.begin();
  collector::iterator iter_end = dataCollectors.end();

  for(;iter != iter_end; iter++) {
     if(iter -> second -> IsContained(zLow,zUp)) return iter -> second;
  }

  return 0;
}


std::ostream& DataManager::os() {

  return storageManager -> GetOutputFileStream();
}
