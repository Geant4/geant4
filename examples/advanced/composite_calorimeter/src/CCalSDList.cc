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

nameType CCalSDList::getCaloSDName(int i){
  
  if (i>=theList->getNumberOfCaloSD() || i<0) {
    G4cout << "CCalSDList invalid calo SD no: " << i << " max is "
	   << theList->getNumberOfCaloSD() << G4endl;
    return " ";
  } else 
    return theList->caloSD[i];
}

nameType CCalSDList::getTrackerSDName(int i){

  if (i>=theList->getNumberOfTrackerSD() || i<0) {
    G4cout << "CCalSDList invalid tracker SD no: " << i << " max is "
	   << theList->getNumberOfTrackerSD() << G4endl;
    return " ";
  }   
  else 
    return theList->trackerSD[i];
}

      
int CCalSDList::getNumberOfCaloSD(){
  
  return theList->caloSD.size();
}

int CCalSDList::getNumberOfTrackerSD(){

  return theList->trackerSD.size();
}
     

CCalSDList& CCalSDList::operator=(CCalSDList&){

  return *this;
}  
