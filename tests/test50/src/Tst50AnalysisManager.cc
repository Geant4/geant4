
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
// $Id: Tst50AnalysisManager.cc,v 1.11 2003-02-07 13:27:49 guatelli Exp $
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
#include "G4RunManager.hh"
Tst50AnalysisManager* Tst50AnalysisManager::instance = 0;

Tst50AnalysisManager::Tst50AnalysisManager() : 
  aFact(0), theTree(0), histFact(0), tupFact(0),p_Primary(0)
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
  delete p_Primary;
}

Tst50AnalysisManager* Tst50AnalysisManager::getInstance()
{
  if (instance == 0) instance = new Tst50AnalysisManager;
  return instance;
}


void Tst50AnalysisManager::book() 
{
  
  G4RunManager* runManager = G4RunManager::GetRunManager();
  p_Primary =
(Tst50PrimaryGeneratorAction*)(runManager->GetUserPrimaryGeneratorAction());

 G4double initial_energy= p_Primary->GetInitialEnergy();

 h1= histFact->createHistogram1D("10","Energy Deposit X event",100.*initial_energy ,0.,initial_energy*2.);

 h2=histFact->createHistogram1D("20","Primary transmitted particle energy/initial_energy",1000. ,0.,1.);
 h3=histFact->createHistogram1D("30","Primary backscattered  particle energy/initial_energy",1000. ,0.,1.);
 
h4=histFact->createHistogram1D("40","angle of backscattered particles",80.*2, 80.,190.);

 h5=histFact->createHistogram1D("50","angle of transmitted  particles",100.*2,0.,100.);
 h6=histFact->createHistogram2D("60","angle, energy of bremmstrahlung gamma",180.*4,0., 180.,100.*initial_energy,0.,initial_energy);
h7= histFact->createHistogram1D("70","Primary  Energy Deposit X event",100.*initial_energy ,0.,initial_energy*2.);

h8= histFact->createHistogram1D("80","Secondary Energy Deposit",100.*initial_energy ,0.,initial_energy*2.);
 // in questo istogramma  metto il deposito di energia di ogni evento nel target

 std::string columnNames = "double initial_energy; double energy; double angle";
 std::string options = "";
 if (tupFact) ntuple = tupFact->create("1","1",columnNames, options);
 // check for non-zero ...

 std::string columnNames2 = "double energy; double angle"; 
 if (tupFact) ntuple2 = tupFact->create("2","2",columnNames2, options);
 // check for non-zero ...
 if (ntuple && ntuple2) G4cout<<"The Ntuple is non-zero"<<G4endl;

 
}

 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void Tst50AnalysisManager::energy_deposit(G4double En)
{ 
   h1->fill(En);
}
void Tst50AnalysisManager::energy_depositPrimary(G4double En3)
{
  h7->fill(En3);
}
void Tst50AnalysisManager::energy_depositSecondary(G4double En4)
{
  h8->fill(En4);
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
void Tst50AnalysisManager::angle_energy_gamma(G4double angle, G4double energy)
{

  h6->fill(angle, energy,1);
}
void Tst50AnalysisManager::fill_dataBrem(G4double energy,G4double angle)
{
  if (ntuple2 == 0) {
    cout << "AAAAAAAGH" << endl;
    return;
  }
  if (ntuple2)
    {
  G4int ien2 = ntuple2->findColumn("energy" );
  G4int iangle2 = ntuple2->findColumn("angle" );
 
  ntuple2->fill(ien2,energy);
  ntuple2->fill(iangle2,angle);
    }
 ntuple2->addRow();
}
void Tst50AnalysisManager::fill_data(G4double initial_en, G4double en, G4double angle)
{
 
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











