
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
// $Id: Tst50AnalysisManager.cc,v 1.6 2003-01-28 08:57:49 guatelli Exp $
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
#include "Tst50PrimaryGeneratorAction.hh"

Tst50AnalysisManager* Tst50AnalysisManager::instance = 0;

Tst50AnalysisManager::Tst50AnalysisManager() : 
  aFact(0), theTree(0), histFact(0), tupFact(0)
  // tupFact(0)
  

{
  //build up  the  factories
   aFact = AIDA_createAnalysisFactory();

   AIDA::ITreeFactory *treeFact = aFact->createTreeFactory();
 
 
  
 //parameters for the TreeFactory
 
  std::string fileName="Test50.hbk";
  theTree = treeFact->create(fileName,"hbook",false, true);

  delete treeFact; 

 
  histFact = aFact->createHistogramFactory( *theTree );
  tupFact  = aFact->createTupleFactory    ( *theTree );


}



Tst50AnalysisManager::~Tst50AnalysisManager() 
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

Tst50AnalysisManager* Tst50AnalysisManager::getInstance()
{
  if (instance == 0) instance = new Tst50AnalysisManager;
  return instance;
}


void Tst50AnalysisManager::book() 
{
  
 


  // h1= histFact->createHistogram1D("10","Energy Deposit",200. ,0.,2.);

 h2=histFact->createHistogram1D("20","Primary transmitted particle energy/initial_energy",1000. ,0.,1.);
 h3=histFact->createHistogram1D("30","Primary backscattered  particle energy/initial_energy",1000. ,0.,1.);
 
h4=histFact->createHistogram1D("40","angle of backscattered particles",80.*2, 80.,190.);

 h5=histFact->createHistogram1D("50","angle of transmitted  particles",100.*2,0.,100.);

 // in questo istogramma  metto il deposito di energia di ogni evento nel target

 std::string columnNames = "double initial_energy; double energy; double angle";
 std::string options = "";
 if (tupFact) ntuple = tupFact->create("1","1",columnNames, options);
 // check for non-zero ...
 if (ntuple) G4cout<<"The Ntuple is non-zero"<<G4endl;

 
}

 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void Tst50AnalysisManager::energy_deposit(G4double En)
{
  //  h1->fill(En);
}
void Tst50AnalysisManager::energy_transmitted(G4double En2)
{
  h2->fill(En2);
}

void Tst50AnalysisManager::energy_backscatter(G4double En_back)
{

  h3->fill(En_back);
}
void Tst50AnalysisManager::angleB(G4double angle)
{

  h4->fill(angle);
}
void Tst50AnalysisManager::angleT(G4double angle)
{

  h5->fill(angle);
}



void Tst50AnalysisManager::fill_data(G4double initial_en, G4double en, G4double angle)
{

  //ntuple = dynamic_cast<ITuple *> ( theTree->find("1") );

  if (ntuple == 0) {
    cout << "AAAAAAAGH" << endl;
    return;
  }
  if (ntuple)
    {

  G4int ien = ntuple->findColumn("initial_energy" );
  G4int ien2 = ntuple->findColumn( "energy" );
  G4int iangle = ntuple->findColumn( "angle" );
  ntuple->fill(ien, initial_en);// fill ( int column, double value )
  ntuple->fill(ien2,en);
  ntuple->fill(iangle,angle);

  ntuple->addRow();
    }
}

void Tst50AnalysisManager::finish() 
{  
  // write all histograms to file
  theTree->commit();

  // close (will again commit)
  theTree->close();

}
#endif











