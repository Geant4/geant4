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
