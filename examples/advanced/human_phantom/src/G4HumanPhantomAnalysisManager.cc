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
//
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
// 
#ifdef  G4ANALYSIS_USE
#include <stdlib.h>
#include <fstream>
#include "G4HumanPhantomAnalysisManager.hh"
#include "G4ios.hh"
#include <AIDA/AIDA.h>
#include "G4RunManager.hh"

G4HumanPhantomAnalysisManager* G4HumanPhantomAnalysisManager::instance = 0;

G4HumanPhantomAnalysisManager::G4HumanPhantomAnalysisManager() 
  :  aFact(0), treeFact(0),theTree(0), histogramFactory(0),
     histogramParticlePath(0), projectionXY(0), projectionYZ(0), 
     projectionZX(0), energy(0), innerBreast(0)
{ 
  aFact = AIDA_createAnalysisFactory();
  treeFact = aFact -> createTreeFactory();
}

G4HumanPhantomAnalysisManager::~G4HumanPhantomAnalysisManager() 
{
  delete innerBreast;
  innerBreast = 0;

  delete energy;
  energy = 0;

  delete projectionXY;
  projectionXY = 0;

  delete projectionYZ;
  projectionYZ = 0;

  delete projectionZX;
  projectionZX = 0;

  delete histogramParticlePath;
  histogramParticlePath = 0;
  
  delete histogramFactory;
  histogramFactory = 0;

  delete treeFact;
  treeFact = 0;

  delete theTree;
  theTree = 0;

  delete aFact;
  aFact = 0;
}

G4HumanPhantomAnalysisManager* G4HumanPhantomAnalysisManager::getInstance()
{
  if (instance == 0) instance = new G4HumanPhantomAnalysisManager;
  return instance;
}

void G4HumanPhantomAnalysisManager::book() 
{
  G4String fileName = "G4HumanPhantom.hbk"; 
  theTree = treeFact->create(fileName,"hbook",false, true);
      
  histogramFactory = aFact -> createHistogramFactory( *theTree );
   
  G4double histogramDimension = 20.;

  G4int bins = 3600; 

  G4int binsProj = 300; 

  histogramParticlePath = histogramFactory->createHistogram1D
   ("10","Particle Path Distribution",bins,0.,histogramDimension);

  projectionXY = histogramFactory->createHistogram2D
   ("20","Particle Projection XY",binsProj,-histogramDimension,histogramDimension,binsProj,-histogramDimension,histogramDimension);

  projectionYZ = histogramFactory->createHistogram2D
   ("30","Particle Projection YZ",binsProj,-histogramDimension,histogramDimension,binsProj,-histogramDimension,histogramDimension);

  projectionZX = histogramFactory->createHistogram2D
   ("40","Particle Projection ZX",binsProj,-histogramDimension,histogramDimension,binsProj,-histogramDimension,histogramDimension);

  energy = histogramFactory->createHistogram2D
   ("50","Energy Deposit in Body Part", 1000,0.,30.,500,0.,100.);

  innerBreast = histogramFactory->createHistogram2D("100", "Edep(MeV) in innerBreast, x= slice, y= sector",
						    11, -0.5, 10.5,
						    11, -0.5, 10.5); 

  G4cout<<"Booking the histograms"<<G4endl;

 }

void  G4HumanPhantomAnalysisManager::particlePath(G4double path)
{
  histogramParticlePath -> fill(path);
}

void G4HumanPhantomAnalysisManager::particleProjectionXY(G4double xx,G4double yy)
{  
 projectionXY -> fill(xx,yy,1.);
}

void G4HumanPhantomAnalysisManager::particleProjectionYZ(G4double yy,G4double zz)
{  
 projectionYZ -> fill(yy,zz,1.);
}

void G4HumanPhantomAnalysisManager::particleProjectionZX(G4double zz,G4double xx)
{  
 projectionZX -> fill(zz,xx,1.);
}

void G4HumanPhantomAnalysisManager::bodypartEnergyDep(G4double bpID, G4double eDep)
{
  energy -> fill(bpID,eDep);
}

void G4HumanPhantomAnalysisManager::innerBreastEnergyDep(G4int slice, G4int sector, G4double edep)
{
  // G4cout << "analisis " << slice << " "<< sector << " "<< edep << G4endl;
  innerBreast -> fill(slice,sector, edep);
}

void G4HumanPhantomAnalysisManager::finish() 
{  
  theTree -> commit();
  theTree -> close();
  G4cout<<"closing the hbk file"<<G4endl;
}
#endif
