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
    cout << "CCalSDList invalid calo SD no: " << i << " max is "
	 << theList->getNumberOfCaloSD() << endl;
    return " ";
  } else 
    return theList->caloSD[i];
}

nameType CCalSDList::getTrackerSDName(int i){

  if (i>=theList->getNumberOfTrackerSD() || i<0) {
    cout << "CCalSDList invalid tracker SD no: " << i << " max is "
	 << theList->getNumberOfTrackerSD() << endl;
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
