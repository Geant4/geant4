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
// $Id: HadrontherapyAnalysisManager.cc,v 1.0
//
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
// --------------------------------------------------------------

#ifdef  G4ANALYSIS_USE
#include <stdlib.h>
#include <fstream>
#include "HadrontherapyAnalysisManager.hh"
#include "G4ios.hh"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IManagedObject.h"
#include "AIDA/IAnalysisFactory.h"
#include "AIDA/IHistogramFactory.h"
#include "AIDA/ITupleFactory.h"
#include "AIDA/ITreeFactory.h"
#include "AIDA/ITree.h"
#include "AIDA/ITuple.h"

HadrontherapyAnalysisManager* HadrontherapyAnalysisManager::instance = 0;

HadrontherapyAnalysisManager::HadrontherapyAnalysisManager() : 
  aFact(0), theTree(0), histFact(0), tupFact(0), h1(0), ntuple(0)
{
  //build up  the  factories
   aFact = AIDA_createAnalysisFactory();
   AIDA::ITreeFactory *treeFact = aFact->createTreeFactory();

 //parameters for the TreeFactory
 
  G4String fileName="protontherapy.hbk";
  theTree = treeFact->create(fileName,"hbook",false, true);
  delete treeFact;
  histFact = aFact->createHistogramFactory( *theTree );
  tupFact  = aFact->createTupleFactory    ( *theTree );
}

HadrontherapyAnalysisManager::~HadrontherapyAnalysisManager() 
{ 
  delete ntuple;
  ntuple=0;

  delete h1;
  h1=0;
  
  delete tupFact;
  tupFact=0;

  delete histFact;
  histFact=0;

  delete theTree;
  theTree=0;

  delete aFact;
  aFact = 0;
}

HadrontherapyAnalysisManager* HadrontherapyAnalysisManager::getInstance()
{
  if (instance == 0) instance = new HadrontherapyAnalysisManager;
  return instance;
}

void HadrontherapyAnalysisManager::book() 
{
  
  h1 = histFact -> createHistogram1D("10","slice, energy",
                                     40*4 ,-2.,40. );
 G4String columnNames = "int slice; double En;";
 G4String options = "";
 if (tupFact) ntuple = tupFact->create("1","1",columnNames, options);
}

void HadrontherapyAnalysisManager::Energy_Dep(G4double slice, G4double En)
{
  if (ntuple)
    {
      G4int islice = ntuple -> findColumn("slice" );
      G4int iEn = ntuple -> findColumn("En" );
      
      ntuple -> fill(islice,slice);
      ntuple -> fill(iEn, En);
    }
  ntuple -> addRow(); 
}

void HadrontherapyAnalysisManager::Energy_Event(G4int slice, G4double En)
{
  h1 -> fill(slice,En);
}

void HadrontherapyAnalysisManager::finish() 
{  
 // write all histograms to file
 theTree -> commit();
 G4cout<<" commit done"<<G4endl;
 
 // close (will again commit)
 theTree ->close();
 G4cout<<" close() done"<<G4endl;
}
#endif











