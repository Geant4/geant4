
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
// $Id: BrachyAnalysisManager.cc,v 1.7 2002-11-18 15:18:37 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include <stdlib.h>
#include "g4std/fstream"
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
  aFact(0), theTree(0), histFact(0), tupFact(0)
  

{
  //build up  the  factories
   aFact = AIDA_createAnalysisFactory();

   AIDA::ITreeFactory *treeFact = aFact->createTreeFactory();
 
 
  
 //parameters for the TreeFactory
 
  std::string fileName="brachytherapy.hbk";
  theTree = treeFact->create(fileName,"hbook",false, true);

  delete treeFact;
 
  histFact = aFact->createHistogramFactory( *theTree );
  tupFact  = aFact->createTupleFactory    ( *theTree );
 
}



BrachyAnalysisManager::~BrachyAnalysisManager() 
{ 
  delete tupFact;
  tupFact=0;

   delete histFact;
  histFact=0;

  delete theTree;
  histFact=0;

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
  
  h1 = histFact->createHistogram2D("10","Energy, pos",300 ,-150.,150.,300,-150.,150.);


 h2= histFact->createHistogram1D("20","Initial Energy", 100,0.,4.);

 std::string columnNames = "float energy, float x,float y, float z";
 std::string options = "";
 if (tupFact) ntuple = tupFact->create("1","1",columnNames, options);
 // check for non-zero ...
 if (ntuple) G4cout<<"The Ntuple is non-zero"<<G4endl;
}
 

 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void BrachyAnalysisManager::fill_Tuple(G4double xx,G4double yy, G4double zz,G4float en)
{

  //ntuple = dynamic_cast<ITuple *> ( theTree->find("1") );

  if (ntuple == 0) {
    cout << "AAAAAAAGH" << endl;
    return;
  }

  ntuple->fill(1, en);// fill ( int column, double value )
  ntuple->fill(2, xx);
  ntuple->fill(3, yy);
  ntuple->fill(4, zz);

  ntuple->addRow();

}

void BrachyAnalysisManager::hist(G4double x,G4double z, G4float enn)
{ 
  h1->fill(x,z,enn);
 
}

void BrachyAnalysisManager::Spectrum(G4double Init_En)
{
  h2->fill(Init_En);
}

void BrachyAnalysisManager::finish() 
{  
  // write all histograms to file
  theTree->commit();

  // close (will again commit)
  theTree->close();

}












