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
//    *********************************
//    *                               *
//    *      CellAnalysisManager.cc   *
//    *                               *
//    *********************************
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//	   Barbara Mascialino (Barbara.Mascialino@ge.infn.it)
//
// History:
// -----------
// 20 September 2006   S. Guatelli, B. Mascialino   1st implementation
//
// -------------------------------------------------------------------
 
#include <stdlib.h>
#include <fstream>
#include "CellAnalysisManager.hh"

#include "G4ios.hh"
#include <AIDA/AIDA.h>
#include "G4RunManager.hh"


CellAnalysisManager* CellAnalysisManager::instance = 0;

CellAnalysisManager::CellAnalysisManager() 
  :  aFact(0), treeFact(0),theTree(0), histogramFactory(0),
     primary_energy(0), histogramEnergyDeposit(0), primary_outgoing_energy(0),
     profile(0)
{ 
  aFact = AIDA_createAnalysisFactory();
  treeFact = aFact -> createTreeFactory();
}

CellAnalysisManager::~CellAnalysisManager() 
{ 
  delete profile;
  profile = 0;
  
  delete primary_outgoing_energy;
  primary_outgoing_energy = 0;

  delete  histogramEnergyDeposit;
  histogramEnergyDeposit = 0;

  delete primary_energy;
  primary_energy = 0;

  delete histogramFactory;
  histogramFactory = 0;

  delete treeFact;
  treeFact = 0;

  delete theTree;
  theTree = 0;

  delete aFact;
  aFact = 0;
}

CellAnalysisManager* CellAnalysisManager::getInstance()
{
  if (instance == 0) instance = new CellAnalysisManager;
  return instance;
}

void CellAnalysisManager::book(G4String name) 
{
  // The hbook file is created
  theTree = treeFact -> create(name + ".hbk", "hbook",false, true);
  
 histogramFactory = aFact -> createHistogramFactory( *theTree );

 // Definition of the histograms
 primary_energy =  histogramFactory -> createHistogram1D("10", 
							 "Initial energy of primary particles (MeV) ",
                                                         1000, 0., 100.);

 histogramEnergyDeposit = histogramFactory -> createHistogram2D
  ("20","Energy Deposit",
   100, -50. , 50., // x bin number, x min, x max
   100, -50. , 50.);  // y bin number, y min, y max


 primary_outgoing_energy =  histogramFactory -> createHistogram1D("30", 
							 "Energy (MeV) of the primary particles outgoing the target ",
                                                         1000, 0., 100.);

 profile = histogramFactory -> createHistogram1D("40", "Energy deposit (MeV) profile along z axis (mm)", 100, -50., 50.);

 }

void CellAnalysisManager::primaryparticle_energy(G4double energy)
{
  primary_energy -> fill(energy);
}


void CellAnalysisManager::FillEnergyDeposit(G4double xx, G4double yy, G4double energyDep)
{
  // the energy deposit is the weight of the 2D histo
  histogramEnergyDeposit -> fill(xx, yy, energyDep);
}

void CellAnalysisManager::primary_energy_outgoing(G4double energy)
{
  primary_outgoing_energy -> fill(energy);
}

void CellAnalysisManager::FillProfile(G4double zz, G4double energy)
{
  profile -> fill(zz, energy);
}

void CellAnalysisManager::finish() 
{  
  theTree -> commit();
  theTree -> close();
}












