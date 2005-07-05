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
//    *    Tst51AnalysisManager.cc *
//    *                             *
//    *******************************
//
// $Id: Tst51AnalysisManager.cc,v 1.1 2005-07-05 11:06:27 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------
 
#include <stdlib.h>
#include <fstream>
#include "Tst51AnalysisManager.hh"
#include "G4ios.hh"
#include <AIDA/AIDA.h>
#include "G4RunManager.hh"


Tst51AnalysisManager* Tst51AnalysisManager::instance = 0;

Tst51AnalysisManager::Tst51AnalysisManager() 
  :  aFact(0), treeFact(0),theTree(0), histogramFactory(0),
     h1(0), h2(0), h3(0), h4(0), tupFact(0), ntuple(0)  
{ 
  aFact = AIDA_createAnalysisFactory();
  treeFact = aFact -> createTreeFactory();
}

Tst51AnalysisManager::~Tst51AnalysisManager() 
{ 
   delete ntuple;
  ntuple = 0;

  delete tupFact;
  tupFact = 0;

  delete energyPostStep;
  energyPostStep = 0;

  delete angularPostStep;
  angularPostStep = 0;
 
  delete h4;
  h4 = 0; 

  delete h3;
  h3 = 0;

  delete h2;
  h2 = 0;

  delete h1;
  h1 = 0;

  delete histogramFactory;
  histogramFactory = 0;

  delete treeFact;
  treeFact = 0;

  delete theTree;
  theTree = 0;

  delete aFact;
  aFact = 0;
}

Tst51AnalysisManager* Tst51AnalysisManager::getInstance()
{
  if (instance == 0) instance = new Tst51AnalysisManager;
  return instance;
}

void Tst51AnalysisManager::book() 
{

  G4cout<<"Booking test51.hbk"<< G4endl;
 theTree = treeFact -> create("test51.hbk","hbook",false, true);
  
 histogramFactory = aFact -> createHistogramFactory( *theTree );
 
 tupFact = aFact -> createTupleFactory( *theTree );

 h1 = histogramFactory -> createHistogram1D("30", "Angular distribution of trans photons", 180, 0., 180.); 
 
 h2 = histogramFactory -> createHistogram1D("40", "Energy distribution of trans photons", 100, 0., 80.);
 
 h3 = histogramFactory -> createHistogram1D("50", "Angular distribution of back photons", 180, 0., 180.);
 
 h4 = histogramFactory -> createHistogram1D("60", "Energy distribution of back photons", 100, 0., 80.);
 
 angularPostStep = histogramFactory -> createHistogram1D("80", "angular distribution, post step", 180, 0., 180.);
 
 energyPostStep = histogramFactory -> createHistogram1D("90", "energy distribution post step", 100, 0., 80.);

 //ntuple
 G4String columnNames = "double energy; double angle";
 G4String options = "";

 if (tupFact) ntuple = tupFact -> create("1","1", columnNames,options);

}

void Tst51AnalysisManager::angularDistributionTransmittedGamma(G4double angle)
{
  h1 -> fill(angle);
}
void Tst51AnalysisManager::energyDistributionTransmittedGamma(G4double energy)
{
  h2 -> fill(energy);
}
void Tst51AnalysisManager::angularDistributionBackGamma(G4double angle)
{
  h3 -> fill(angle);
}

void Tst51AnalysisManager::energyDistributionBackGamma(G4double energy)
{
  h4 -> fill(energy);
}

void Tst51AnalysisManager::angularDistributionPostStep(G4double angle)
{
  angularPostStep -> fill(angle);
}

void Tst51AnalysisManager::energyDistributionPostStep(G4double energy)
{
  energyPostStep -> fill(energy);
}
void Tst51AnalysisManager::fillNtuple(G4double energy, G4double angle)
{
if (ntuple)
    {
      G4int iEnergy = ntuple -> findColumn("energy");
      G4int iAngle = ntuple -> findColumn("angle");
     
      ntuple -> fill(iEnergy, energy);
      ntuple -> fill(iAngle, angle); 
    }

  ntuple -> addRow(); 

}
void Tst51AnalysisManager::finish() 
{  
  theTree -> commit();
  theTree -> close();
  G4cout<<"Committing test51.hbk"<<G4endl;
}












