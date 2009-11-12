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
// Code developed by:
//  S.Guatelli
//
//    *******************************
//    *                             *
//    *    BrachyAnalysisManager.cc *
//    *                             *
//    *******************************
//
// $Id: BrachyAnalysisManager.cc,v 1.22 2009-11-12 10:32:59 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifdef G4ANALYSIS_USE

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
aFact(0), theTree(0), histFact(0), tupFact(0), h1(0),
  h2(0), h3(0), ntuple(0)

{
  
  // Instantiate the factories
  // The factories manage the analysis objects
  aFact = AIDA_createAnalysisFactory();

  AIDA::ITreeFactory *treeFact = aFact -> createTreeFactory(); 
  
  // Definition of the output file
  G4String fileName = "brachytherapy.root";
  theTree = treeFact -> create(fileName,"ROOT",false, true);

  delete treeFact;
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
  // Instantiate the histogram and ntuple factories
  histFact = aFact -> createHistogramFactory( *theTree );
  tupFact  = aFact -> createTupleFactory    ( *theTree ); 
 
  // Creating a 2D histogram
  // Energy deposit in the plane containing the source
  h1 = histFact -> createHistogram2D("10","Energy, pos", //histoID,histo name
				     300 ,-150.,150.,    //bins'number,xmin,xmax 
                                     300,-150.,150.);    //bins'number,ymin,ymax 
  //creating a 1D histograms
  // Histogram containing the initial energy (MeV) of the photons delivered by the radioactive core
  h2 = histFact -> createHistogram1D("20","Initial Energy", //histoID, histo name 
				     1000,0.,1.);            //bins' number, xmin, xmax
   
  // Histogram containing the energy deposit in the plane containing the source, along the axis 
  // perpendicular to the source main axis
  h3 = histFact -> createHistogram1D("30","Energy deposit  Distribution", 
				     300,-150.,150.); //bins' number, xmin, xmax

  //defining the ntuple columns' name 
  std::string columnNames = "double energy, x , y , z ";
  std::string options = "";
  
  //creating a ntuple
  if (tupFact) ntuple = tupFact -> create("1","1",columnNames, options);
  // check for non-zero ...
  if (ntuple) G4cout<<"The Ntuple is non-zero"<<G4endl;
}
 
void BrachyAnalysisManager::FillNtupleWithEnergy(G4double xx,
                                                 G4double yy, 
                                                 G4double zz,
                                                 G4double en)
{
  if (ntuple == 0) 
   {
     G4cout << "AAAAAAAGH..... The Ntuple is 0" << G4endl;
     return;
    }
 
  // Fill the ntuple
  
  G4int indexX = ntuple -> findColumn( "x" );
  G4int indexY = ntuple -> findColumn( "y" );
  G4int indexZ = ntuple -> findColumn( "z" );
  G4int indexEnergy = ntuple -> findColumn( "energy" );

  ntuple -> fill(indexEnergy, en);// method: fill ( int column, double value )
  ntuple -> fill(indexX, xx);
  ntuple -> fill(indexY, yy);
  ntuple -> fill(indexZ, zz);

  ntuple->addRow();
}

void BrachyAnalysisManager::FillHistogramWithEnergy(G4double x,
                                                    G4double z, 
                                                    G4double energyDeposit)
{
  // 2DHistogram: energy deposit in a voxel which center is fixed in position (x,z)  
  h1 -> fill(x,z,energyDeposit);
}

void BrachyAnalysisManager::PrimaryParticleEnergySpectrum(G4double primaryParticleEnergy)
{
 // 1DHistogram: energy spectrum of primary particles  
  h2 -> fill(primaryParticleEnergy);
  return;
}

void BrachyAnalysisManager::DoseDistribution(G4double x,G4double energy)
{
  // 1DHistogram: energy spectrum of primary particles  
  h3 -> fill(x, energy);
}

void BrachyAnalysisManager::finish() 
{  

 // write all histograms to file ...
  theTree -> commit();

  // close (will again commit) ...
  theTree -> close();
}
#endif











