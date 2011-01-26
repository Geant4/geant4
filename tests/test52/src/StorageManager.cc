
#include "StorageManager.hh"
#include <string>
#include <iostream>
#include <AIDA/AIDA.h>


StorageManager* StorageManager::instance = 0;


StorageManager* StorageManager::Instance(std::string fileNameBase) {
 
  if(instance == 0) {
    instance = new StorageManager(fileNameBase);
  }
  
  return instance;
}


void StorageManager::Destroy() {

  if(!instance == 0) {

     delete instance;
     instance = 0;
  }
}


StorageManager::StorageManager(std::string fileNameBase) {

  analysisFactory = AIDA_createAnalysisFactory();

  AIDA::ITreeFactory* treeFactory = analysisFactory->createTreeFactory();

  std::string fileNameHisto = fileNameBase + "_histograms.xml";
  tree = treeFactory->create(fileNameHisto,"xml",false, true,"uncompressed"); 
  delete treeFactory;

  histogramFactory = analysisFactory->createHistogramFactory(*tree); 

  std::string fileNameTxt = fileNameBase + ".txt";

  os.open(fileNameTxt.c_str());
}


StorageManager::~StorageManager() {

  os.close();

  tree -> commit();
  tree -> close();

  delete histogramFactory;
  delete tree;
  delete analysisFactory;
}


AIDA::IHistogramFactory* StorageManager::GetIHistogramFactory() {

  return histogramFactory;
}


std::ostream& StorageManager::GetOutputFileStream() {

  return os;
}
