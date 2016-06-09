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
// $Id: RemSimAnalysisManager.cc,v 1.7 2004/11/23 11:43:21 guatelli Exp $
//
// Author:Susanna Guatelli, guatelli@ge.infn.it 
//
#ifdef  G4ANALYSIS_USE 
#include <stdlib.h>
#include <fstream>
#include "RemSimAnalysisManager.hh"
#include "RemSimAnalysisMessenger.hh"
#include "G4ios.hh"
#include <AIDA/AIDA.h>
#include "G4RunManager.hh"


RemSimAnalysisManager* RemSimAnalysisManager::instance = 0;

RemSimAnalysisManager::RemSimAnalysisManager() 
  :  aFact(0), treeFact(0), theTree(0), dataPointFactory(0),
     histogramFactory(0), dataPoint(0), energyDeposit(0),
     primary(0), secondaryDeposit(0), primaryInitialE(0), 
     primaryInitialEout(0), initialE(0), 
     initialEout(0), shape(0), energyShape(0)
{ 
  
  aFact = AIDA_createAnalysisFactory();
  messenger = new  RemSimAnalysisMessenger(this); 
  fileFormat = "hbook";
}

RemSimAnalysisManager::~RemSimAnalysisManager() 
{ 
  delete messenger;

  delete energyShape;
  energyShape = 0;

  delete shape;
  shape = 0;
 
  delete initialEout;
  initialEout = 0;
  
  delete initialE;
  initialE = 0;
  
  delete primaryInitialEout;
  primaryInitialEout = 0;
 
  delete primaryInitialE;
  primaryInitialE = 0;
 
  delete secondaryDeposit;
  secondaryDeposit = 0;

  delete primary;
  primary = 0;

  delete energyDeposit;
  energyDeposit = 0;

  delete dataPoint;
  dataPoint = 0; 

  delete histogramFactory;
  histogramFactory = 0;

  delete dataPointFactory;
  dataPointFactory = 0;

  delete theTree;
  theTree = 0;

  delete treeFact;
  treeFact = 0;

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
 treeFact = aFact -> createTreeFactory();
  
  if (fileFormat == "hbook")
    { 
     theTree = treeFact -> create("remsim.hbk", "hbook", false, true);
     G4cout << "The format of the output file is hbook" << G4endl;
    }

  else if (fileFormat == "xml")
    {
     theTree = treeFact -> create("remsim.xml","xml",false, true,
                                  "uncompress");
     G4cout<< "The format of the output file is xml" << G4endl;
    }

  histogramFactory = aFact -> createHistogramFactory(*theTree);
  energyDeposit = histogramFactory -> createHistogram1D("10",
                                                        "Energy Deposit",
                                                        30,// number of bins
					                0.,//xmin
                                                        30.);//xmax 

  primary = histogramFactory -> createHistogram1D("20",
				                "Initial energy of primary particles", 
                                                200000,0.,100000.);
 
  secondaryDeposit = histogramFactory -> createHistogram1D("30",
					 "EnergyDeposit given by secondaries", 
					 30,0.,30.);

 primaryInitialE = histogramFactory -> createHistogram1D("40",
			   "Initial energy of primaries reaching the phantom", 
                           200000,0.,100000.);


 primaryInitialEout = histogramFactory -> createHistogram1D("50",
					       "Initial energy of primaries ougoing the phantom", 
                                                200000,0.,100000.);

 initialE = histogramFactory -> createHistogram1D("60",
					       "Energy of primaries reaching the phantom", 
                                                200000,0.,100000.);


 initialEout = histogramFactory -> createHistogram1D("70",
					       "Energy of primaries outgoing the phantom", 
                                                200000,0.,100000.);

 shape =  histogramFactory -> createHistogram2D("80",
					       "Shape", 
                                                300,-150.,150.,
                                                300, -150.,150.);

 energyShape = histogramFactory -> createHistogram2D("90", 
                                                        "energyDepShape",
                                                        300, -150.,150.,
                                                        300, -150.,150.);

}

void RemSimAnalysisManager::energyDepositStore(G4int layer, G4double eDep)
{
  energyDeposit -> fill(layer,eDep);
}

void RemSimAnalysisManager::primaryParticleEnergyDistribution(G4double energy)
{
  primary -> fill(energy);
}

void RemSimAnalysisManager::SecondaryEnergyDeposit(G4int i, G4double EnergyDep)
{
  secondaryDeposit -> fill(i,EnergyDep); 
}

void RemSimAnalysisManager::PrimaryInitialEnergyIn(G4double EnergyDep)
{
  primaryInitialE -> fill(EnergyDep); 
}

void RemSimAnalysisManager::PrimaryInitialEnergyOut(G4double EnergyDep)
{
  primaryInitialEout -> fill(EnergyDep); 
}

void RemSimAnalysisManager::PrimaryEnergyIn(G4double EnergyDep)
{
  initialE -> fill(EnergyDep); 
}

void RemSimAnalysisManager::PrimaryEnergyOut(G4double EnergyDep)
{
  initialEout -> fill(EnergyDep); 
}

void RemSimAnalysisManager::particleShape(G4double x, G4double y)
{
  shape -> fill(x,y); 
}

void RemSimAnalysisManager::energyDepShape(G4double x, G4double y, G4double energyDep)
{
  energyShape -> fill(x,y, energyDep); 
}

void RemSimAnalysisManager:: SetFormat(G4String format)
{ 
  fileFormat = format;
 
  if (fileFormat != "hbook" && fileFormat != "xml")
  G4cout << fileFormat << "is not available" << G4endl; 
}

void RemSimAnalysisManager::finish() 
{  
  theTree -> commit();
  theTree -> close();
}
#endif











