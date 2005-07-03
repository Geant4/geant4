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
// ********************************************************************
// $Id: MedLinacAnalysisManager.cc,v 1.7 2005-07-03 23:27:37 mpiergen Exp $
//
//
// Code developed by: M. Piergentili
//
//
#ifdef  G4ANALYSIS_USE
#include <stdlib.h>
#include <fstream>
#include "MedLinacAnalysisManager.hh"

#include "G4ios.hh"

#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"

#include "AIDA/IManagedObject.h"
#include "AIDA/IAnalysisFactory.h"
#include "AIDA/IHistogramFactory.h"
#include "AIDA/ITupleFactory.h"
#include "AIDA/ITreeFactory.h"
#include "AIDA/ITree.h"
#include "AIDA/ITuple.h"

MedLinacAnalysisManager* MedLinacAnalysisManager::instance = 0;

MedLinacAnalysisManager::MedLinacAnalysisManager(): 
  aFact(0), theTree(0), histFact(0),h1(0),h2(0),h3(0), h4(0),h5(0),h6(0),h7(0),h8(0)
  
{
  //build up  the  factories
  aFact = AIDA_createAnalysisFactory();

  AIDA::ITreeFactory *treeFact = aFact->createTreeFactory(); 
 
  //parameters for the TreeFactory
 
  std::string fileName = "MedLinacDiane.xml";
  theTree = treeFact->create(fileName,"xml",false, true,"uncompress");

  delete treeFact;
 
  histFact = aFact->createHistogramFactory( *theTree );

}
MedLinacAnalysisManager::~MedLinacAnalysisManager() 
{ 
  delete histFact;
  histFact = 0;

  delete theTree;
  histFact = 0;

  delete aFact;
  aFact = 0;
}

MedLinacAnalysisManager* MedLinacAnalysisManager::getInstance()
{
  if (instance == 0) instance = new MedLinacAnalysisManager;
  return instance;
}

void MedLinacAnalysisManager::book() 
{
   //creating an other 1D histogram ...
  h1 = histFact->createHistogram1D("14010","PDD5", //histoID,histo name
				    60,-150.0,150.0); //bins' number, zmin, zmax
  //creating an other 1D histogram ...
  h2 = histFact->createHistogram1D("14020","Flatness build-up5", //histoID,histo name
				    60,-150.0,150.0); //bins' number, zmin, zmax

  //creating an other 1D histogram ...
  h3 = histFact->createHistogram1D("14030","Flatness 50mm 5", //histoID,histo name
				    60,-150.0,150.0); //bins' number, zmin, zmax


  //creating an other 1D histogram ...
  h4 = histFact->createHistogram1D("14040","PDD", //histoID,histo name
				    150,-150.0,150.0); //bins' number, zmin, zmax
  //creating an other 1D histogram ...
  h5 = histFact->createHistogram1D("14050","Flatness build-up", //histoID,histo name
				    150,-150.0,150.0); //bins' number, zmin, zmax

 //creating an other 1D histogram ...
  h6 = histFact->createHistogram1D("14060","Flatness 50mm", //histoID,histo name
				    150,-150.0,150.0); //bins' number, zmin, zmax
 //creating an other 1D histogram ...
  h7 = histFact->createHistogram1D("14070","Flatness 100mm", //histoID,histo name
				    150,-150.0,150.0); //bins' number, zmin, zmax
 //creating an other 1D histogram ...
  h8 = histFact->createHistogram1D("14080","Flatness 200mm", //histoID,histo name
				    150,-150.0,150.0); //bins' number, zmin, zmax
}

void MedLinacAnalysisManager::FillHistogram1WithEnergy(G4double z, 
                                                    G4float energyDeposit)
{
  //G4cout << " fill HISTO1-------------"<<G4endl;
  //1DHistrogram: energy deposit in a voxel which center is fixed in position (0,0,z)  
  h1->fill(z,energyDeposit);
}

void MedLinacAnalysisManager::FillHistogram2WithEnergy(G4double x, 
                                                    G4float energyDeposit)
{
  //G4cout << " fill HISTO2-------------"<<G4endl;
  //1DHistrogram: energy deposit in a voxel which center is fixed in position (x,0,135mm)  
  h2->fill(x,energyDeposit);
} 
void MedLinacAnalysisManager::FillHistogram3WithEnergy(G4double x, 
                                                    G4float energyDeposit)
{
  //G4cout << " fill HISTO2-------------"<<G4endl;
  //1DHistrogram: energy deposit in a voxel which center is fixed in position (x,0,135mm)  
  h3->fill(x,energyDeposit);
} 
void MedLinacAnalysisManager::FillHistogram4WithEnergy(G4double z, 
                                                    G4float energyDeposit)
{
  //G4cout << " fill HISTO4-------------"<<G4endl;
  //1DHistrogram: energy deposit in a voxel which center is fixed in position (0,0,z)  
  h4->fill(z,energyDeposit);
}

void MedLinacAnalysisManager::FillHistogram5WithEnergy(G4double x, 
                                                    G4float energyDeposit)
{
  //G4cout << " fill HISTO5-------------"<<G4endl;
  //1DHistrogram: energy deposit in a voxel which center is fixed in position (x,0,135mm)  
  h5->fill(x,energyDeposit);
}

void MedLinacAnalysisManager::FillHistogram6WithEnergy(G4double x, 
                                                    G4float energyDeposit)
{
  //G4cout << " fill HISTO6-------------"<<G4endl;
  //1DHistrogram: energy deposit in a voxel which center is fixed in position (x,0,80mm)
  h6->fill(x,energyDeposit);
}

void MedLinacAnalysisManager::FillHistogram7WithEnergy(G4double x, 
                                                    G4float energyDeposit)
{
  //G4cout << " fill HISTO7-------------"<<G4endl; 
  //1DHistrogram: energy deposit in a voxel which center is fixed in position (x,0,30mm)
  h7->fill(x,energyDeposit);
}

void MedLinacAnalysisManager::FillHistogram8WithEnergy(G4double x, 
                                                    G4float energyDeposit)
{
  //G4cout << " fill HISTO8-------------"<<G4endl;
 //1DHistrogram: energy deposit in a voxel which center is fixed in position (x,0,-70mm)
  h8->fill(x,energyDeposit);
}


void MedLinacAnalysisManager::finish() 
{  
  // write all histograms to file ...
  theTree->commit();

  // close (will again commit) ...
  theTree->close();
}
#endif











