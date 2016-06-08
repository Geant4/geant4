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
// File: CCalSDList.hh
// Description: Records name of all SD objects and classify them into
//              CALO and Tracker SD
//              Singleton class
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalSDList_h
#define CCalSDList_h 1

#include "g4std/vector"
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
  G4std::vector<nameType> caloSD;
  G4std::vector<nameType> trackerSD;
  
private:
  CCalSDList& operator=(CCalSDList&);   

};

#endif
