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
// $Id: RemSimAnalysisManager.cc,v 1.3 2004-03-12 10:55:55 guatelli Exp $
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
     dataPoint(0), energyDeposit(0),histo1(0),histo2(0),histo3(0),trasmission1(0), trasmission2(0), trasmission3(0),trasmission12(0), trasmission22(0), 
     trasmission32(0),trasmission13(0), trasmission23(0), trasmission33(0), 
     neutron(0), photon(0),electron(0),hadron(0), primary(0)
{ 
  aFact = AIDA_createAnalysisFactory();
  treeFact = aFact -> createTreeFactory();
}

RemSimAnalysisManager::~RemSimAnalysisManager() 
{
  delete primary;
  primary =0;

  delete hadron;
  hadron =0;

  delete electron;
  electron =0;

  delete photon;
  photon=0;
 
  delete neutron;
  neutron =0;

  delete trasmission33;
  trasmission33=0;
  
  delete trasmission23;
  trasmission23=0;

  delete trasmission13;
  trasmission13=0;
 
  delete trasmission32;
  trasmission32=0;
  
  delete trasmission22;
  trasmission22=0;

  delete trasmission12;
  trasmission12=0;

  delete trasmission3;
  trasmission3=0;
  
  delete trasmission2;
  trasmission2=0;

  delete trasmission1;
  trasmission1=0;

  delete histo3;
  histo3 =0;

  delete histo2;
  histo2 =0;
 
  delete histo1;
  histo1 =0;

  delete energyDeposit;
  energyDeposit =0;

  delete dataPoint;
  dataPoint = 0; 

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
  //theTree = treeFact -> create("remsim.xml","xml",false, true,"uncompress");
  theTree = treeFact -> create("remsim.hbk","hbook",false, true);

  histogramFactory = aFact -> createHistogramFactory( *theTree );
  energyDeposit = histogramFactory->createHistogram1D("10","Energy Deposit",50,0.,50.); 
   histo1 = histogramFactory->createHistogram2D("20","SpatialInfo at slab =0",
                                                500,-2.50,2.50,
                                                500,-2.50,2.50);
 
   histo2 = histogramFactory->createHistogram2D("30","SpatialInfo at slab=3",
                                                500,-2.50,2.50,
                                                500,-2.50,2.50);
 
   histo3 = histogramFactory->createHistogram2D("40","SpatialInfo3 at slab=6",
                                                500,-2.50,2.50,
                                                500,-2.50,2.50); 
  
   trasmission1 = histogramFactory->createHistogram1D("50",
                                                      "e+,e- impinging on the phantom ",
                                                      200000,0.,10000.);

   trasmission2 = histogramFactory->createHistogram1D("60",
                                                       "hadrons impinging on the phantom ",
                                                       200000,0.,10000.);
   trasmission3 = histogramFactory->createHistogram1D("70",
                                                      "photons impinging on the phantom ",
                                                      200000,0.,10000.);

   trasmission12 = histogramFactory->createHistogram1D("80",
                                                       "e+,e- impinging in the phantom at 2.cm",                                                             20000,0.,10000.);

   trasmission22 = histogramFactory->createHistogram1D("90",
                                                       "hadrons impinging in the phantom at 2.cm",
                                                       200000,0.,10000.);

   trasmission32 = histogramFactory->createHistogram1D("100",
                                                       "photons impinging in the phantom at 2.cm",
                                                      200000,0.,10000.);

   trasmission13 = histogramFactory->createHistogram1D("110",
                                                       "e+,e- outgoing the phantom",                                                                         200000,0.,10000.);

   trasmission23 = histogramFactory->createHistogram1D("120",
                                                       "hadrons outgoing the phantom", 
                                                       200000,0.,100000.);

   trasmission33 = histogramFactory->createHistogram1D("130",
                                                       "photons outgoing the phantom", 
                                                       200000,0.,100000.);

 neutron = histogramFactory->createHistogram1D("140",
					       "EnergyDistribution of secondary neutrons", 
                                                       200000,0.,100000.);
 photon = histogramFactory->createHistogram1D("150",
					       "EnergyDistribution of photons", 
                                                200000,0.,100000.);
 electron = histogramFactory->createHistogram1D("160",
					     "EnergyDistribution of electrons",
                                             200000,0.,100000.);

 hadron = histogramFactory->createHistogram1D("170",
					     "EnergyDistribution of secondary hadrons (pions, kaons)",
                                             200000,0.,100000.);

 primary = histogramFactory->createHistogram1D("180",
					       "Energy of primary particles", 
                                                200000,0.,100000.);


}

void RemSimAnalysisManager:: energyDeposit1(G4double xx,G4double yy)
{
  histo1 -> fill (xx,yy,1);
}
void RemSimAnalysisManager:: energyDeposit2(G4double xx,G4double yy)
{
  histo2 -> fill (xx,yy,1);
}
void RemSimAnalysisManager:: energyDeposit3(G4double xx,G4double yy)
{
  histo3 -> fill (xx,yy,1);
}

void RemSimAnalysisManager:: energyDepositStore(G4int layer, G4double eDep)
{
  energyDeposit -> fill(layer,eDep);
}

void RemSimAnalysisManager:: leptonsEnergySpectrum1(G4double energy)
{
  trasmission1 -> fill(energy);
}

void RemSimAnalysisManager:: hadronEnergySpectrum1(G4double energy)
{
  trasmission2 -> fill(energy);
}
void RemSimAnalysisManager:: gammaEnergySpectrum1(G4double energy)
{
  trasmission3 -> fill(energy);
}
void RemSimAnalysisManager:: leptonsEnergySpectrum2(G4double energy)
{
  trasmission12 -> fill(energy);
}

void RemSimAnalysisManager:: hadronEnergySpectrum2(G4double energy)
{
  trasmission22 -> fill(energy);
}
void RemSimAnalysisManager:: gammaEnergySpectrum2(G4double energy)
{
  trasmission32 -> fill(energy);
}
void RemSimAnalysisManager:: leptonsEnergySpectrum3(G4double energy)
{
  trasmission13 -> fill(energy);
}

void RemSimAnalysisManager:: hadronEnergySpectrum3(G4double energy)
{
  trasmission23 -> fill(energy);
}
void RemSimAnalysisManager:: gammaEnergySpectrum3(G4double energy)
{
  trasmission33 -> fill(energy);
}
void RemSimAnalysisManager::neutronEnergyDistribution(G4double energy)
{
  neutron-> fill(energy);
}
void RemSimAnalysisManager::photonEnergyDistribution(G4double energy)
{
  photon-> fill(energy);
}
void RemSimAnalysisManager::electronEnergyDistribution(G4double energy)
{
  electron-> fill(energy);
}
void RemSimAnalysisManager::hadronEnergyDistribution(G4double energy)
{
  hadron-> fill(energy);
}


void RemSimAnalysisManager::primaryParticleEnergyDistribution(G4double energy)
{
  primary-> fill(energy);
}

void RemSimAnalysisManager::finish() 
{  
  theTree -> commit();
  theTree -> close();
}
#endif











