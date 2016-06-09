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
// Code developed by:
//  S.Guatelli
//
//    *******************************
//    *                             *
//    *    BrachyAnalysisManager.cc *
//    *                             *
//    *******************************
//
// $Id: BrachyAnalysisManager.cc,v 1.13 2004/03/11 15:38:42 guatelli Exp $
// GEANT4 tag $Name: geant4-06-01 $
//
#ifdef  G4ANALYSIS_USE
#include <stdlib.h>
#include <fstream>
#include "BrachyAnalysisManager.hh"

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

BrachyAnalysisManager* BrachyAnalysisManager::instance = 0;

BrachyAnalysisManager::BrachyAnalysisManager() : 
  aFact(0), theTree(0), histFact(0), tupFact(0),h1(0),h2(0),ntuple(0)
  
{
  //build up  the  factories
  aFact = AIDA_createAnalysisFactory();

  AIDA::ITreeFactory *treeFact = aFact->createTreeFactory(); 
 
  //parameters for the TreeFactory
 
  std::string fileName = "brachytherapy.hbk";
  theTree = treeFact->create(fileName,"hbook",false, true);

  delete treeFact;
 
  histFact = aFact->createHistogramFactory( *theTree );
  tupFact  = aFact->createTupleFactory    ( *theTree ); 
}

BrachyAnalysisManager::~BrachyAnalysisManager() 
{ 
  delete tupFact;
  tupFact = 0;

  delete histFact;
  histFact = 0;

  delete theTree;
  histFact = 0;

  delete aFact;
  aFact = 0;
}

BrachyAnalysisManager* BrachyAnalysisManager::getInstance()
{
  if (instance == 0) instance = new BrachyAnalysisManager;
  return instance;
}

void BrachyAnalysisManager::book() 
{
  //creating a 1D histogram ...
  h1 = histFact->createHistogram2D("10","Energy, pos", //histoID,histo name
				    300 ,-150.,150.,   //bins'number,xmin,xmax 
                                    300,-150.,150.    );//bins'number,ymin,ymax 
  //creating a 2D histogram ...
  h2 = histFact->createHistogram1D("20","Initial Energy", //histoID,histo name 
				  100,0.,1.); //bins' number, xmin, xmax
  
  h3 = histFact->createHistogram1D("30","Dose Distribution", 
				  300,-150.,150.); //bins' number, xmin, xmax

  //defining the ntuple columns' name ...
  std::string columnNames = "float energy; float x; float y; float z";
  std::string options = "";
  
  //creating a ntuple ...
  if (tupFact) ntuple = tupFact->create("1","1",columnNames, options);
  // check for non-zero ...
  if (ntuple) G4cout<<"The Ntuple is non-zero"<<G4endl;
}
 
void BrachyAnalysisManager::FillNtupleWithEnergy(G4double xx,
                                                 G4double yy, 
                                                 G4double zz,
                                                 G4float en)
{
  if (ntuple == 0) 
   {
     G4cout << "AAAAAAAGH..... The Ntuple is 0" << G4endl;
     return;
    }
  
  G4int indexX = ntuple->findColumn( "x" );
  G4int indexY = ntuple->findColumn( "y" );
  G4int indexZ = ntuple->findColumn( "z" );
  G4int indexEnergy = ntuple->findColumn( "energy" );
  ntuple->fill(indexEnergy, en);// fill ( int column, double value )
  ntuple->fill(indexX, xx);
  ntuple->fill(indexY, yy);
  ntuple->fill(indexZ, zz);

  ntuple->addRow();
}

void BrachyAnalysisManager::FillHistogramWithEnergy(G4double x,
                                                    G4double z, 
                                                    G4double energyDeposit)
{
  //2DHistrogram: energy deposit in a voxel which center is fixed in position (x,z)  
  h1->fill(x,z,energyDeposit);
}

void BrachyAnalysisManager::PrimaryParticleEnergySpectrum(G4double primaryParticleEnergy)
{
  //1DHisotgram: energy spectrum of primary particles  
  h2->fill(primaryParticleEnergy);
}
void BrachyAnalysisManager::DoseDistribution(G4double x,G4double energy)
{
  //1DHisotgram: energy spectrum of primary particles  
  h3->fill(x, energy);
}

void BrachyAnalysisManager::finish() 
{  
  // write all histograms to file ...
  theTree->commit();

  // close (will again commit) ...
  theTree->close();
}
#endif











