//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
//
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//    ****************************************************
//    *      UltraAnalysisManager.cc
//    ****************************************************
//
//    Class used for analysis procedures if the environment variable 
//    G4ANALYSIS_USE is set. AIDA is supported.
//    Based on BrachyAnalysisManager class developed by S.Guatelli.
//
#ifdef G4ANALYSIS_USE 
#include <stdlib.h>
#include <fstream>
#include "UltraAnalysisManager.hh"
#include "G4ios.hh"

#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"

#include "AIDA/IManagedObject.h"
#include "AIDA/IAnalysisFactory.h"
#include "AIDA/IHistogramFactory.h"
#include "AIDA/ITupleFactory.h"
#include "AIDA/ITreeFactory.h"
#include "AIDA/ITree.h"
#include "AIDA/ITuple.h"

UltraAnalysisManager* UltraAnalysisManager::instance = 0;

UltraAnalysisManager::UltraAnalysisManager() : 
  aFact(0), theTree(0), histFact(0), h1(0),h2(0)
  
{
  //build up  the  factories
  aFact = AIDA_createAnalysisFactory();

  AIDA::ITreeFactory* treeFact = aFact->createTreeFactory();
  //parameters for the TreeFactory

  // "hbook" for HBOOK histograms
  // "xml" for AIDA histograms
  // "root" for ROOT histograms

  std::string fileName = "ultra.aida";
  theTree = treeFact->create(fileName,"xml",false, true);

  delete treeFact;
 
  histFact = aFact->createHistogramFactory( *theTree );

}

UltraAnalysisManager::~UltraAnalysisManager() 
{ 

  delete histFact;
  histFact = 0;

  delete theTree;
  histFact = 0;

  delete aFact;
  aFact = 0;
}

UltraAnalysisManager* UltraAnalysisManager::getInstance()
{
  if (instance == 0) instance = new UltraAnalysisManager;
  return instance;
}

void UltraAnalysisManager::book() 
{

  h1 = histFact->createHistogram1D("10","Optical photons energy (eV)", //histoID,histo name 
				  500,0.,5.); //bins' number, xmin, xmax
  
  h2 = histFact->createHistogram1D("20","Number of Detected Photons", 
				  10,0.,10.); //bins' number, xmin, xmax
}
 

void UltraAnalysisManager::FillHistogram(G4int i, G4double f){

  if(i == 1) h1->fill(f);
  if(i == 2) h2->fill(f);

}

void UltraAnalysisManager::finish() 
{  
  // write all histograms to file ...
  theTree->commit();

  // close (will again commit) ...
  theTree->close();
}


#endif







