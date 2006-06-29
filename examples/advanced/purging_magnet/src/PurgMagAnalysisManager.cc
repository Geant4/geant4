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
//  S.Larsson
//
//    ************************************
//    *                                  *
//    *    PurgMagAnalysisManager.cc     *
//    *                                  *
//    ************************************
//
// $Id: PurgMagAnalysisManager.cc,v 1.3 2006-06-29 16:06:09 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifdef  G4ANALYSIS_USE
#include <stdlib.h>
#include <fstream>
#include "PurgMagAnalysisManager.hh"

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

PurgMagAnalysisManager* PurgMagAnalysisManager::instance = 0;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagAnalysisManager::PurgMagAnalysisManager() : 
  aFact(0), theTree(0), tupFact(0)
  

{
  // Build the factories
  aFact = AIDA_createAnalysisFactory();
  
  AIDA::ITreeFactory *treeFact = aFact->createTreeFactory();
  
  
  // Parameters for the TreeFactory
  
  std::string fileName="purgmag.hbk";
  theTree = treeFact->create(fileName,"hbook",false, true);
  
  delete treeFact;
  
  tupFact  = aFact->createTupleFactory    ( *theTree );
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagAnalysisManager::~PurgMagAnalysisManager() 
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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagAnalysisManager* PurgMagAnalysisManager::getInstance()
{
  if (instance == 0) instance = new PurgMagAnalysisManager;
  return instance;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PurgMagAnalysisManager::book() 
{
  
  // N-tuple
  std::string options = "";

  // Electrons
  std::string columnNames1 = "double ex; double ey; double ez; double ee; double epx; double epy; double epz";
  if (tupFact) ntuple1 = tupFact->create("1","1",columnNames1, options);
  // check for non-zero ...
  if (ntuple1) G4cout<<"N-tuple 1 is non-zero"<<G4endl;

  // Photons (gamma)
  std::string columnNames2 = "double gx; double gy; double gz; double ge; double gpx; double gpy; double gpz";
  if (tupFact) ntuple2 = tupFact->create("2","2",columnNames2, options);
  // check for non-zero ...
  if (ntuple2) G4cout<<"N-tuple 2 is non-zero"<<G4endl;

  // Positrons
  std::string columnNames3 = "double px; double py; double pz; double pe; double ppx; double ppy; double ppz";
  if (tupFact) ntuple3 = tupFact->create("3","3",columnNames3, options);
  // check for non-zero ...
  if (ntuple3) G4cout<<"N-tuple 3 is non-zero"<<G4endl;

   }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function fills a N-tuple with position, energy and momentum of 
// electrons entering the measurement volume. 
void PurgMagAnalysisManager::fill_Tuple_Electrons(G4double ex, G4double ey, G4double ez,     // Position
					     G4double ee,                               // Energy
					     G4double epx, G4double epy, G4double epz)  // Momentum
{

  if (ntuple1 == 0) {
    G4cout << "N-tuple 1 is zero " << "\n";
    return;
  }
  
  int iex = ntuple1->findColumn( "ex" );
  int iey = ntuple1->findColumn( "ey" );
  int iez = ntuple1->findColumn( "ez" );
  int iee = ntuple1->findColumn( "ee" );
  int iepx = ntuple1->findColumn( "epx" );
  int iepy = ntuple1->findColumn( "epy" );
  int iepz = ntuple1->findColumn( "epz" );
  ntuple1->fill(iex, ex);                  // fill ( int column, double value )
  ntuple1->fill(iey, ey);
  ntuple1->fill(iez, ez);
  ntuple1->fill(iee, ee);
  ntuple1->fill(iepx, epx);
  ntuple1->fill(iepy, epy);
  ntuple1->fill(iepz, epz);

  ntuple1->addRow();
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function fills a N-tuple with position, energy and momentum of 
// photons entering the measurement volume. 
void PurgMagAnalysisManager::fill_Tuple_Gamma(G4double gx, G4double gy, G4double gz,     // Position 
					 G4double ge,                               // Energy
					 G4double gpx, G4double gpy, G4double gpz)  // Momentum
{

  if (ntuple2 == 0) {
    G4cout << "N-tuple 2 is zero" << "\n";
    return;
  }

  int igx = ntuple2->findColumn( "gx" );
  int igy = ntuple2->findColumn( "gy" );
  int igz = ntuple2->findColumn( "gz" );
  int ige = ntuple2->findColumn( "ge" );
  int igpx = ntuple2->findColumn( "gpx" );
  int igpy = ntuple2->findColumn( "gpy" );
  int igpz = ntuple2->findColumn( "gpz" );
  ntuple2->fill(igx, gx);                   // fill ( int column, double value )
  ntuple2->fill(igy, gy);
  ntuple2->fill(igz, gz);
  ntuple2->fill(ige, ge);
  ntuple2->fill(igpx, gpx);
  ntuple2->fill(igpy, gpy);
  ntuple2->fill(igpz, gpz);

  ntuple2->addRow();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function fills a N-tuple with position, energy and momentum of 
// positrons entering the measurement volume. 
void PurgMagAnalysisManager::fill_Tuple_Positrons(G4double px, G4double py, G4double pz,     // Position 
					     G4double pe,                               // Energy
					     G4double ppx, G4double ppy, G4double ppz)  // Momentum
{

  if (ntuple3 == 0) {
    G4cout << "N-tuple 3 is zero" << "\n";
    return;
  }

  int ipx = ntuple3->findColumn( "px" );
  int ipy = ntuple3->findColumn( "py" );
  int ipz = ntuple3->findColumn( "pz" );
  int ipe = ntuple3->findColumn( "pe" );
  int ippx = ntuple3->findColumn( "ppx" );
  int ippy = ntuple3->findColumn( "ppy" );
  int ippz = ntuple3->findColumn( "ppz" );
  ntuple3->fill(ipx, px);                  // fill ( int column, double value )
  ntuple3->fill(ipy, py);
  ntuple3->fill(ipz, pz);
  ntuple3->fill(ipe, pe);
  ntuple3->fill(ippx, ppx);
  ntuple3->fill(ippy, ppy);
  ntuple3->fill(ippz, ppz);

  ntuple3->addRow();

}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void PurgMagAnalysisManager::finish() 
{  
  // Writes all histograms to file
  theTree->commit();

  // Close (will again commit)
  theTree->close();
}
#endif

