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
#ifdef  G4ANALYSIS_USE
#include <stdlib.h>
#include <fstream>
#include "HadrontherapyAnalysisManager.hh"

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

HadrontherapyAnalysisManager* HadrontherapyAnalysisManager::instance = 0;

HadrontherapyAnalysisManager::HadrontherapyAnalysisManager() : 
  aFact(0), theTree(0), dataPointFactory(0), energyDepositDataPoint(0)
{
  //build up  the  factories
  aFact = AIDA_createAnalysisFactory();

  AIDA::ITreeFactory *treeFact = aFact -> createTreeFactory(); 
 
  //parameters for the TreeFactory
 
  G4String fileName = "hadrontherapy.xml";
  theTree = treeFact -> create(fileName,"xml",false, true, "uncompress");

  delete treeFact; 
}

HadrontherapyAnalysisManager::~HadrontherapyAnalysisManager() 
{ 
  delete energyDepositDataPoint;
  energyDepositDataPoint = 0;

  delete dataPointFactory;
  dataPointFactory = 0;

  delete theTree;
  theTree = 0;

  delete aFact;
  aFact = 0;
}

HadrontherapyAnalysisManager* HadrontherapyAnalysisManager::getInstance()
{
  if (instance == 0) instance = new HadrontherapyAnalysisManager;
  return instance;
}

void HadrontherapyAnalysisManager::book() 
{
  dataPointFactory = aFact -> createDataPointSetFactory(*theTree); 
  energyDepositDataPoint = dataPointFactory -> create("EnergyDep3D","x, y, z, E",4); 
}

void HadrontherapyAnalysisManager:: energyDeposit3D(G4int PointNumber,
                                                    G4int x,
                                                    G4int y,
                                                    G4int z, 
                                                    G4double energyDeposit)
{
  // Fill the DataSetPoint
  energyDepositDataPoint -> addPoint();  
  
  AIDA::IDataPoint* point = energyDepositDataPoint -> point(PointNumber);
  
  AIDA::IMeasurement* coordinateX = point -> coordinate( 0 );
  coordinateX->setValue(x);
  
  AIDA::IMeasurement* coordinateY = point -> coordinate( 1 );
  coordinateY -> setValue(y);

  AIDA::IMeasurement* coordinateZ = point -> coordinate( 2 );
  coordinateZ -> setValue(z);

  AIDA::IMeasurement* energy = point -> coordinate( 3 );
  energy -> setValue(energyDeposit);
}

void HadrontherapyAnalysisManager::finish() 
{  
  // write all histograms to file ...
  theTree->commit();

  // close (will again commit) ...
  theTree->close();
}
#endif











