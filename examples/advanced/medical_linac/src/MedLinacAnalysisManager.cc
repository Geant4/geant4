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
  aFact(0), theTree(0), histFact(0),h1(0),h2(0),h3(0),h4(0),h5(0)
  
{
  //build up  the  factories
  aFact = AIDA_createAnalysisFactory();

  AIDA::ITreeFactory *treeFact = aFact->createTreeFactory(); 
 
  //parameters for the TreeFactory
 
  std::string fileName = "MedLinac.hbk";
  theTree = treeFact->create(fileName,"hbook",false, true);

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
  //creating a 2D histogram ...
  h1 = histFact->createHistogram2D("10","Energy, pos x z", //histoID,histo name
				    300 ,-170.,170.,   //bins'number,xmin,xmax 
                                    300,-170.,170.    );//bins'number,zmin,zmax 
  //creating a 1D histogram ...
  h2 = histFact->createHistogram1D("20","Initial Energy", //histoID,histo name 
				  300,4.0,8.0); //bins' number, xmin, xmax
  //creating an other 2D histogram ...
  h3 = histFact->createHistogram2D("30","Energy, pos x y", //histoID,histo name
				    300 ,-170.,170.,   //bins'number,xmin,xmax 
                                    300,-170.,170.    );//bins'number,ymin,ymax 

  //creating an other 1D histogram ...
  h4 = histFact->createHistogram1D("40","PDD", //histoID,histo name
				    100,-170.0,170.0); //bins' number, zmin, zmax
  //creating an other 1D histogram ...
  h5 = histFact->createHistogram1D("50","Flatness", //histoID,histo name
				    100,-170.0,170.0); //bins' number, zmin, zmax

}
 
void MedLinacAnalysisManager::FillHistogram1WithEnergy(G4double x,
                                                    G4double z, 
                                                    G4float energyDeposit)
{
  G4cout << " fill HISTO1-------------"<<G4endl;
  //2DHistrogram: energy deposit in a voxel which center is fixed in position (x,z)  
  h1->fill(x,z,energyDeposit);
}               
                                                                        
//Units:   the energy deposit is in MeV;
//       x, y, z in mm for both 2D-histogram and ntuple  

void MedLinacAnalysisManager::PrimaryParticleEnergySpectrum(G4double pEnergy)
{
  G4cout << " fill histo primaryparticle-------------"<<G4endl;
  //1DHistogram: energy spectrum of primary particles  
  h2->fill(pEnergy/MeV);
}
void MedLinacAnalysisManager::FillHistogram3WithEnergy(G4double x,
                                                    G4double y, 
                                                    G4float energyDeposit)
{
  G4cout << " fill HISTO3-------------"<<G4endl;
  //2DHistrogram: energy deposit in a voxel which center is fixed in position (x,y)  
  h3->fill(x,y,energyDeposit);
}


void MedLinacAnalysisManager::FillHistogram4WithEnergy(G4double z, 
                                                    G4float energyDeposit)
{
  G4cout << " fill HISTO4-------------"<<G4endl;
  //1DHistrogram: energy deposit in a voxel which center is fixed in position (0,0,z)  
  h4->fill(z,energyDeposit);
}

void MedLinacAnalysisManager::FillHistogram5WithEnergy(G4double x, 
                                                    G4float energyDeposit)
{
  G4cout << " fill HISTO5-------------"<<G4endl;
  //1DHistrogram: energy deposit in a voxel which center is fixed in position (x,0,295mm)  
  h5->fill(x,energyDeposit);
}


void MedLinacAnalysisManager::finish() 
{  
  // write all histograms to file ...
  theTree->commit();

  // close (will again commit) ...
  theTree->close();
}
#endif











