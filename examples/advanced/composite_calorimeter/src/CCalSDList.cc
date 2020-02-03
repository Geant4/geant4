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
///////////////////////////////////////////////////////////////////////////////
// File: CCalSDList.cc
// Description: Records names of all SD objects
///////////////////////////////////////////////////////////////////////////////
#include "CCalSDList.hh"

CCalSDList* CCalSDList::theList = 0;

CCalSDList::CCalSDList(){}
CCalSDList::~CCalSDList(){ delete theList;}

CCalSDList* CCalSDList::getInstance(){

  if (theList == 0) 
    theList = new CCalSDList;
  return theList;
}
     
void CCalSDList::addCalo(nameType name){
    
  theList->caloSD.push_back(name);
}

void CCalSDList::addTracker(nameType name){

  theList->trackerSD.push_back(name);
} 

nameType CCalSDList::getCaloSDName(G4int i){
  
  if (i>=theList->getNumberOfCaloSD() || i<0) {
    G4cout << "CCalSDList invalid calo SD no: " << i << " max is "
           << theList->getNumberOfCaloSD() << G4endl;
    return " ";
  } else 
    return theList->caloSD[i];
}

nameType CCalSDList::getTrackerSDName(G4int i){

  if (i>=theList->getNumberOfTrackerSD() || i<0) {
    G4cout << "CCalSDList invalid tracker SD no: " << i << " max is "
           << theList->getNumberOfTrackerSD() << G4endl;
    return " ";
  }   
  else 
    return theList->trackerSD[i];
}

      
G4int CCalSDList::getNumberOfCaloSD(){
  
  return theList->caloSD.size();
}

G4int CCalSDList::getNumberOfTrackerSD(){

  return theList->trackerSD.size();
}
     

CCalSDList& CCalSDList::operator=(CCalSDList&){

  return *this;
}  
