///////////////////////////////////////////////////////////////////////////////
// File: SDList.cc
// Date: 10.1999 V.Lefebure
// Modifications:
///////////////////////////////////////////////////////////////////////////////

#include "SDList.hh"

SDList* SDList::theList = 0;

SDList::SDList(){}
SDList::~SDList(){ delete theList;}

SDList* SDList::getInstance(){

  if (theList == 0) 
    theList = new SDList;
  return theList;
}
     
void SDList::addCalo(nameType name){
    
  theList->caloSD.push_back(name);
}

void SDList::addTracker(nameType name){

  theList->trackerSD.push_back(name);
} 

nameType SDList::getCaloSDName(int i){
  
  if (i>=theList->getNumberOfCaloSD() || i<0) {
    cout << "SDList invalid calo SD no: " << i << " max is "
	 << theList->getNumberOfCaloSD() << endl;
    return " ";
  } else 
    return theList->caloSD[i];
}

nameType SDList::getTrackerSDName(int i){

  if (i>=theList->getNumberOfTrackerSD() || i<0) {
    cout << "SDList invalid tracker SD no: " << i << " max is "
	 << theList->getNumberOfTrackerSD() << endl;
    return " ";
  }   
  else 
    return theList->trackerSD[i];
}

      
int SDList::getNumberOfCaloSD(){
  
  return theList->caloSD.size();
}

int SDList::getNumberOfTrackerSD(){

  return theList->trackerSD.size();
}
     

SDList& SDList::operator=(SDList&){

  return *this;
}  
