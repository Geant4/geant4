///////////////////////////////////////////////////////////////////////////////
// File: SDList.hh
// Date: 10.1999 
// Description: Records name of all SD objects and classify them into
//              CALO and Tracker SD
//              Singleton class
// Modifications:
///////////////////////////////////////////////////////////////////////////////
#ifndef SDList_h
#define SDList_h 1

#include <vector>
#include "globals.hh"

typedef G4String nameType;

class SDList{
private: 
  SDList();
  ~SDList();

public:
  static SDList* getInstance();
      
public:
     
  void addCalo(nameType name);
  void addTracker(nameType name);
  
  nameType getCaloSDName(int i);
  nameType getTrackerSDName(int i);
      
  int getNumberOfCaloSD();
  int getNumberOfTrackerSD();
  
private:
  static SDList* theList;
  vector<nameType> caloSD;
  vector<nameType> trackerSD;
  
private:
  SDList& operator=(SDList&);   

};

#endif
