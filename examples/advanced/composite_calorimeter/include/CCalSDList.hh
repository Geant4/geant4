///////////////////////////////////////////////////////////////////////////////
// File: CCalSDList.hh
// Description: Records name of all SD objects and classify them into
//              CALO and Tracker SD
//              Singleton class
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalSDList_h
#define CCalSDList_h 1

#include <vector>
#include "globals.hh"

typedef G4String nameType;

class CCalSDList{
private: 
  CCalSDList();
  ~CCalSDList();

public:
  static CCalSDList* getInstance();
      
public:
     
  void addCalo(nameType name);
  void addTracker(nameType name);
  
  nameType getCaloSDName(int i);
  nameType getTrackerSDName(int i);
      
  int getNumberOfCaloSD();
  int getNumberOfTrackerSD();
  
private:
  static CCalSDList* theList;
  vector<nameType> caloSD;
  vector<nameType> trackerSD;
  
private:
  CCalSDList& operator=(CCalSDList&);   

};

#endif
