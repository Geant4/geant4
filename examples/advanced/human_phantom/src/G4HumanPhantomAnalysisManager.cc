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
 
#include <stdlib.h>
#include <fstream>
#include "G4HumanPhantomAnalysisManager.hh"
#include "G4HumanPhantomConstruction.hh"
#include "G4ios.hh"
#include <AIDA/AIDA.h>
#include "G4RunManager.hh"

G4HumanPhantomAnalysisManager* G4HumanPhantomAnalysisManager::instance = 0;

G4HumanPhantomAnalysisManager::G4HumanPhantomAnalysisManager() 
  :  aFact(0), treeFact(0),theTree(0), histogramFactory(0),
     histogramParticlePath(0), projectionXY(0), projectionYZ(0), projectionZX(0)
{ 
  aFact = AIDA_createAnalysisFactory();
  treeFact = aFact -> createTreeFactory();
}

G4HumanPhantomAnalysisManager::~G4HumanPhantomAnalysisManager() 
{
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

  G4cout<<"booking the histograms"<<G4endl;

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


void G4HumanPhantomAnalysisManager::finish() 
{  
  theTree -> commit();
  theTree -> close();
  G4cout<<"closing the hbk file"<<G4endl;
}
