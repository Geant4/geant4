
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
// $Id: Tst50AnalysisManager.cc,v 1.3 2002-12-16 13:50:08 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifdef  G4ANALYSIS_USE

#include <stdlib.h>
#include "g4std/fstream"
#include "Tst50AnalysisManager.hh"

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

Tst50AnalysisManager* Tst50AnalysisManager::instance = 0;

Tst50AnalysisManager::Tst50AnalysisManager() : 
  aFact(0), theTree(0), histFact(0)
  // tupFact(0)
  

{
  //build up  the  factories
   aFact = AIDA_createAnalysisFactory();

   AIDA::ITreeFactory *treeFact = aFact->createTreeFactory();
 
 
  
 //parameters for the TreeFactory
 
  std::string fileName="Test_500kev.hbk";
  theTree = treeFact->create(fileName,"hbook",false, true);

  delete treeFact;
 
  histFact = aFact->createHistogramFactory( *theTree );
  //  tupFact  = aFact->createTupleFactory    ( *theTree );
 
}



Tst50AnalysisManager::~Tst50AnalysisManager() 
{ 
 // delete tupFact;
  //tupFact=0;

   delete histFact;
  histFact=0;

  delete theTree;
  histFact=0;

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
  
 


 h1= histFact->createHistogram1D("10","Energy Deposit", 3000000,0.,0.03);

 // in questo istogramma  metto il deposito di energia di ogni evento nel target
 h2= histFact->createHistogram1D("20","primary_processes", 6000,0.,6.);

 h3= histFact->createHistogram1D("30","transmitted gamma",200,0.,2.);
}
 
void Tst50AnalysisManager::primary_processes(G4int process)
{
 h2->fill(process) ;
}


 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void Tst50AnalysisManager::energy_deposit(G4double En)
{
  h1->fill(En);
}


void Tst50AnalysisManager::trans_particles()
{
 h3->fill(1);

}
void Tst50AnalysisManager::finish() 
{  
  // write all histograms to file
  theTree->commit();

  // close (will again commit)
  theTree->close();

}
#endif











