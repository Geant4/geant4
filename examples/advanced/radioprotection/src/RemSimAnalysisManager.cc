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
//    *******************************
//    *                             *
//    *    RemSimAnalysisManager.cc *
//    *                             *
//    *******************************
//
// $Id: RemSimAnalysisManager.cc,v 1.2 2004-02-03 09:16:46 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------
#ifdef  G4ANALYSIS_USE 
#include <stdlib.h>
#include <fstream>
#include "RemSimAnalysisManager.hh"

#include "G4ios.hh"
#include <AIDA/AIDA.h>
#include "G4RunManager.hh"


RemSimAnalysisManager* RemSimAnalysisManager::instance = 0;

RemSimAnalysisManager::RemSimAnalysisManager() 
  :  aFact(0), treeFact(0),theTree(0),dataPointFactory(0), histogramFactory(0),
  stoppingPowerDataPoint(0),CSDARangeDataPoint(0)
 
{ 
  aFact = AIDA_createAnalysisFactory();
  treeFact = aFact -> createTreeFactory();
}

RemSimAnalysisManager::~RemSimAnalysisManager() 
{ 
  delete CSDARangeDataPoint;
  CSDARangeDataPoint = 0; 

  delete stoppingPowerDataPoint;
  stoppingPowerDataPoint = 0; 

  delete histogramFactory;
  histogramFactory = 0;

  delete treeFact;
  treeFact = 0;

  delete theTree;
  theTree = 0;

  delete aFact;
  aFact = 0;
}

RemSimAnalysisManager* RemSimAnalysisManager::getInstance()
{
  if (instance == 0) instance = new RemSimAnalysisManager;
  return instance;
}

void RemSimAnalysisManager::book() 
{
  theTree = treeFact -> create("test50.xml","xml",false, true,"uncompress");
  
  //Create the factories for dataPoint and histograms
  dataPointFactory = aFact -> createDataPointSetFactory(*theTree); 

  stoppingPowerDataPoint = dataPointFactory -> create("SP","Stopping Power test",2); 
  CSDARangeDataPoint = dataPointFactory -> create ("CSDARange","CSDA Range test",2);
  
 }
void RemSimAnalysisManager::StoppingPower(G4int PointNumber,
                                         G4double primaryParticleEnergy,
                                         G4double SP)
{
  stoppingPowerDataPoint->addPoint();
  AIDA::IDataPoint* point = stoppingPowerDataPoint->point(PointNumber);
  AIDA::IMeasurement* coordinateX = point->coordinate( 0 );
  coordinateX -> setValue(primaryParticleEnergy );
  coordinateX -> setErrorPlus( 0. );
  coordinateX -> setErrorMinus( 0. );
  AIDA::IMeasurement* coordinateY = point->coordinate( 1 );
  coordinateY -> setValue( SP);
  coordinateY -> setErrorPlus( 0. );
  coordinateY -> setErrorMinus( 0. );
}


void RemSimAnalysisManager::CSDARange(G4int PointNumber,
                                     G4double primaryParticleEnergy, 
                                     G4double range)
{
  CSDARangeDataPoint -> addPoint();
  AIDA::IDataPoint* point = CSDARangeDataPoint -> point(PointNumber);
  AIDA::IMeasurement* coordinateX = point -> coordinate( 0 );
  coordinateX -> setValue(primaryParticleEnergy );
  coordinateX -> setErrorPlus( 0. );
  coordinateX -> setErrorMinus( 0. );
  AIDA::IMeasurement* coordinateY = point -> coordinate( 1 );
  coordinateY -> setValue(range);
  coordinateY -> setErrorPlus( 0. );
  coordinateY -> setErrorMinus( 0. );
}


void RemSimAnalysisManager::finish() 
{  
  theTree -> commit();
  theTree -> close();
}
#endif











