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
  :  aFact(0), treeFact(0),theTree(0), histogramFactory(0),tupFact(0),
     ntuple(0), voxelLeftBreast(0), voxelRightBreast(0)
{ 
  aFact = AIDA_createAnalysisFactory();
  treeFact = aFact -> createTreeFactory();
}

G4HumanPhantomAnalysisManager::~G4HumanPhantomAnalysisManager() 
{
  delete voxelRightBreast;
  voxelRightBreast = 0;

  delete voxelLeftBreast;
  voxelLeftBreast = 0;

  delete ntuple;
  ntuple = 0;
  
  delete tupFact;
  tupFact =0;

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
  tupFact  = aFact -> createTupleFactory( *theTree ); 


  voxelLeftBreast = histogramFactory->createHistogram2D("100", 
			      "Edep(MeV) in LeftBreast, x= slice, y= sector",
						    11, -0.5, 10.5,
						    11, -0.5, 10.5); 

  voxelRightBreast = histogramFactory->createHistogram2D("110", 
			      "Edep(MeV) in RightBreast, x= slice, y= sector",
						    11, -0.5, 10.5,
						    11, -0.5, 10.5); 

  
  // Defining the ntuple columns' name 
  std::string columnNames = "int id; float energy";
  std::string options = "";
  
  // Creating a ntuple
  if (tupFact) ntuple = tupFact -> create("1","1",columnNames, options);       
  
  G4cout<<"Booking !!!!"<<G4endl;
 }

void G4HumanPhantomAnalysisManager::bodyPartEnergyDeposit(G4int bodyPartID, 
						      G4double eDep)
{
 if (ntuple == 0) 
   {
     G4cout << "AAAAAAAGH..... The Ntuple is 0" << G4endl;
     return;
    }
 
  // Fill the ntuple
  
  // Each organ is identified with an integer  
  G4int indexX = ntuple -> findColumn( "id" );
  G4int indexEnergy = ntuple -> findColumn( "energy" );

  ntuple -> fill(indexX, bodyPartID);
  ntuple -> fill(indexEnergy, eDep);

  ntuple->addRow();
}

void G4HumanPhantomAnalysisManager::voxelLeftBreastEnergyDeposit(G4int slice, G4int sector, G4double edep)
{
  // G4cout << "analisis " << slice << " "<< sector << " "<< edep << G4endl;
  voxelLeftBreast -> fill(slice,sector, edep);
}

void G4HumanPhantomAnalysisManager::voxelRightBreastEnergyDeposit(G4int slice, G4int sector, G4double edep)
{
  // G4cout << "analisis " << slice << " "<< sector << " "<< edep << G4endl;
  voxelRightBreast -> fill(slice,sector, edep);
}

void G4HumanPhantomAnalysisManager::finish() 
{  
  theTree -> commit();
  theTree -> close();
  G4cout<<"closing the hbk file"<<G4endl;
}
#endif
