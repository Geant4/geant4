
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
// Code developed by:
//  S.Guatelli
//
//    *******************************
//    *                             *
//    *    BrachyAnalysisManager.cc *
//    *                             *
//    *******************************
//

// $Id: Tst50AnalysisManager.cc,v 1.17 2003-05-17 11:59:43 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $


#ifdef  G4ANALYSIS_USE

#include <stdlib.h>
#include "g4std/fstream"
#include "Tst50AnalysisManager.hh"

#include "G4ios.hh"
#include <AIDA/AIDA.h>
#include "G4RunManager.hh"


Tst50AnalysisManager* Tst50AnalysisManager::instance = 0;

Tst50AnalysisManager::Tst50AnalysisManager() : 
  aFact(0), treeFact(0),theTree(0),dpsf(0),dpsa(0),dpsa1(0),dpsa2(0),dpsa3(0)
{ 
  aFact = AIDA_createAnalysisFactory();
  treeFact = aFact->createTreeFactory();
}

Tst50AnalysisManager::~Tst50AnalysisManager() 
{ 
  delete dpsa3;
  dpsa3=0;

  delete dpsa2;
  dpsa2=0; 

  delete dpsa1;
  dpsa1=0; 

  delete dpsa;
  dpsa=0; 

  delete treeFact;
  treeFact=0;

  delete theTree;
  theTree=0;

  delete aFact;
  aFact = 0;
}


Tst50AnalysisManager* Tst50AnalysisManager::getInstance()
{
  if (instance == 0) instance = new Tst50AnalysisManager;
  return instance;
}


void Tst50AnalysisManager::book() 
{ 
  theTree = treeFact->create("test50.xml","xml",false, true,"uncompress");
  dpsf = aFact->createDataPointSetFactory(*theTree); 
  dpsa = dpsf->create("Stopping Power test",2); 
  dpsa1 = dpsf->create ("CSDA Range test",2);
  dpsa2 = dpsf->create ("Transmission test",3);
  dpsa3 = dpsf->create ("Gamma attenuation coefficient test",2); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst50AnalysisManager::attenuation_coeffiecient(G4int PointNumber,G4double energy, G4double coeff, G4double coeff_error )
{
  dpsa3->addPoint();
  AIDA::IDataPoint* point = dpsa3->point(PointNumber);
  AIDA::IMeasurement* mX= point->coordinate( 0 );
  mX->setValue(energy );
  AIDA::IMeasurement* mY= point->coordinate( 1 );
  mY->setValue( coeff );
  mY->setErrorPlus(coeff_error );
  mY->setErrorMinus(coeff_error);
 }
void Tst50AnalysisManager::StoppingPower(G4int PointNumber,G4double energy, G4double SP)
{
  dpsa->addPoint();
  AIDA::IDataPoint* point = dpsa->point(PointNumber);
  AIDA::IMeasurement* mX= point->coordinate( 0 );
  mX->setValue(energy );
  mX->setErrorPlus( 0. );
  mX->setErrorMinus( 0. );
  AIDA::IMeasurement* mY= point->coordinate( 1 );
  mY->setValue( SP);
}
void Tst50AnalysisManager::CSDARange(G4int PointNumber,G4double energy, G4double range)
{
  dpsa1->addPoint();
  AIDA::IDataPoint* point = dpsa1->point(PointNumber);
  AIDA::IMeasurement* mX= point->coordinate( 0 );
  mX->setValue(energy );
  AIDA::IMeasurement* mY= point->coordinate( 1 );
  mY->setValue(range);
}
void Tst50AnalysisManager::trasmission(G4int PointNumber,G4double energy, G4double TransFraction, G4double BackFraction, G4double TransError, G4double BackError)
{
  dpsa2->addPoint();
  AIDA::IDataPoint* point = dpsa2->point(PointNumber);
  AIDA::IMeasurement* mX= point->coordinate( 0 );
  mX->setValue(energy );
  AIDA::IMeasurement* mY= point->coordinate( 1 );
  mY->setValue(TransFraction);
  mY->setErrorPlus(TransError);
  mY->setErrorMinus(TransError);
  AIDA::IMeasurement* mZ= point->coordinate( 2 );
  mZ->setValue(BackFraction);
  mZ->setErrorPlus(BackError);
  mZ->setErrorMinus(BackError);
}

void Tst50AnalysisManager::finish() 
{  
  theTree->commit();
  theTree->close();
}
#endif











