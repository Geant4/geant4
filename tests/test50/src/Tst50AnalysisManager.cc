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
//    *******************************
//    *                             *
//    *    Tst50AnalysisManager.cc *
//    *                             *
//    *******************************
//
// $Id: Tst50AnalysisManager.cc,v 1.29 2010-06-07 10:08:39 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------
#ifdef G4ANALYSIS_USE

#include <stdlib.h>
#include <fstream>
#include "Tst50AnalysisManager.hh"

#include "G4ios.hh"
#include <AIDA/AIDA.h>
#include "G4RunManager.hh"


Tst50AnalysisManager* Tst50AnalysisManager::instance = 0;

Tst50AnalysisManager::Tst50AnalysisManager() 
  :  aFact(0), treeFact(0),theTree(0),dataPointFactory(0), histogramFactory(0),
  stoppingPowerDataPoint(0),CSDARangeDataPoint(0),
  particleTransmissionDataPoint(0),
  gammaAttenuationCoefficientDataPoint(0),histogramEnergyDeposit(0)
{ 
  aFact = AIDA_createAnalysisFactory();
  treeFact = aFact -> createTreeFactory();
}

Tst50AnalysisManager::~Tst50AnalysisManager() 
{ 
  delete gammaAttenuationCoefficientDataPoint;
  gammaAttenuationCoefficientDataPoint = 0;

  delete particleTransmissionDataPoint;
  particleTransmissionDataPoint = 0; 

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

Tst50AnalysisManager* Tst50AnalysisManager::getInstance()
{
  if (instance == 0) instance = new Tst50AnalysisManager;
  return instance;
}

void Tst50AnalysisManager::book(G4String name) 
{
  theTree = treeFact -> create(name + "output" + ".xml",
                               "xml",false, true,"uncompress");
  
  //Create the factories for dataPoint and histograms
  dataPointFactory = aFact -> createDataPointSetFactory(*theTree); 

  stoppingPowerDataPoint = dataPointFactory -> create("SP","Stopping Power test",2); 
  CSDARangeDataPoint = dataPointFactory -> create ("CSDARange","CSDA Range test",2);
  particleTransmissionDataPoint = dataPointFactory -> create ("Transmission test",3);

  // Energy of transmitted/backscattered particles
 particleTransmissionEnergyDataPoint = dataPointFactory -> create ("TransmissionEnergy test",2);

  gammaAttenuationCoefficientDataPoint = dataPointFactory -> create ("Gamma attenuation coefficient test",2);  

 histogramFactory = aFact -> createHistogramFactory( *theTree );
 histogramEnergyDeposit = histogramFactory->createHistogram1D
("10"," Energy Deposit per event",120,0.,1.2); 
 }
void Tst50AnalysisManager::bookHistograms() 
{
 std::string fileName = "test50.hbk";
 theTree = treeFact->create(fileName,"hbook",false, true);
 histogramFactory = aFact -> createHistogramFactory( *theTree );
 histogramEnergyDeposit = histogramFactory->createHistogram1D
("10"," Energy Deposit per event",120,0.,1.2); 
}

void Tst50AnalysisManager::AttenuationGammaCoeffiecient(G4int PointNumber,
                                                        G4double primaryParticleEnergy, 
                                                        G4double gammaAttenuationCoefficient,
                                                        G4double gammaAttenuationCoefficientError )
{
  gammaAttenuationCoefficientDataPoint -> addPoint();
  AIDA::IDataPoint* point = gammaAttenuationCoefficientDataPoint -> point(PointNumber);
  AIDA::IMeasurement* coordinateX = point->coordinate( 0 );
  coordinateX->setValue(primaryParticleEnergy );
  AIDA::IMeasurement* coordinateY = point -> coordinate( 1 );
  coordinateY -> setValue( gammaAttenuationCoefficient );
  coordinateY -> setErrorPlus(gammaAttenuationCoefficientError );
  coordinateY -> setErrorMinus(gammaAttenuationCoefficientError);
 }

void Tst50AnalysisManager::StoppingPower(G4int PointNumber,
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

void Tst50AnalysisManager::CSDARange(G4int PointNumber,
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

void Tst50AnalysisManager::ParticleTransmission(G4int PointNumber,
                                                G4double primaryParticleEnergy,
                                                G4double TransFraction,
                                                G4double BackFraction, 
                                                G4double TransError, 
                                                G4double BackError)
{
  particleTransmissionDataPoint -> addPoint();
  AIDA::IDataPoint* point = particleTransmissionDataPoint -> point(PointNumber);
  AIDA::IMeasurement* coordinateX = point -> coordinate( 0 );
  coordinateX -> setValue(primaryParticleEnergy );
  AIDA::IMeasurement* coordinateY = point -> coordinate( 1 );
  coordinateY -> setValue(TransFraction);
  coordinateY -> setErrorPlus(TransError);
  coordinateY -> setErrorMinus(TransError);
  AIDA::IMeasurement* coordinateZ = point->coordinate( 2 );
  coordinateZ -> setValue(BackFraction);
  coordinateZ -> setErrorPlus(BackError);
  coordinateZ -> setErrorMinus(BackError);
}
void Tst50AnalysisManager::TransmittedEnergy(G4int PointNumber, G4double primaryParticleEnergy,G4double energy)
{
  particleTransmissionEnergyDataPoint-> addPoint();
  AIDA::IDataPoint* point = particleTransmissionEnergyDataPoint -> point(PointNumber);
  AIDA::IMeasurement* coordinateX = point -> coordinate( 0 );
  coordinateX -> setValue(primaryParticleEnergy );
  AIDA::IMeasurement* coordinateY = point -> coordinate( 1 );
  coordinateY -> setValue(energy);
  coordinateY -> setErrorPlus(0.);
  coordinateY -> setErrorMinus(0.);
}

void Tst50AnalysisManager::FillEnergyDeposit(G4double energyDep)
{
  histogramEnergyDeposit ->fill(energyDep);
}

void Tst50AnalysisManager::finish() 
{  
  theTree -> commit();
  theTree -> close();
}

#endif







