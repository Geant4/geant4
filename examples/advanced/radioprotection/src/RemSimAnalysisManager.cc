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
// $Id: RemSimAnalysisManager.cc,v 1.5 2004/05/22 12:57:06 guatelli Exp $
//
// Author:Susanna Guatelli, guatelli@ge.infn.it 
//
#ifdef  G4ANALYSIS_USE 
#include <stdlib.h>
#include <fstream>
#include "RemSimAnalysisManager.hh"

#include "G4ios.hh"
#include <AIDA/AIDA.h>
#include "G4RunManager.hh"


RemSimAnalysisManager* RemSimAnalysisManager::instance = 0;

RemSimAnalysisManager::RemSimAnalysisManager() 
  :  aFact(0), treeFact(0), theTree(0), dataPointFactory(0),
     histogramFactory(0), dataPoint(0), energyDeposit(0),
     primary(0), secondaryDeposit(0), primaryInitialE(0), 
     primaryInitialEout(0), initialE(0), 
     initialEout(0)
{ 
  aFact = AIDA_createAnalysisFactory();
  treeFact = aFact -> createTreeFactory();
}

RemSimAnalysisManager::~RemSimAnalysisManager() 
{ 
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
  // Define the hbook file
  theTree = treeFact -> create("remsim.hbk","hbook",false, true);
  
  // Create histogram factory
  histogramFactory = aFact -> createHistogramFactory(*theTree);

  // Histograms

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

void RemSimAnalysisManager::finish() 
{  
  theTree -> commit();
  theTree -> close();
}
#endif











