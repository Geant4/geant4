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

class CCalSDList
{
private: 
  CCalSDList();
  ~CCalSDList();

public:
  static CCalSDList* getInstance();
      
public:
     
  void addCalo(nameType name);
  void addTracker(nameType name);
  
  nameType getCaloSDName(G4int i);
  nameType getTrackerSDName(G4int i);
      
  G4int getNumberOfCaloSD();
  G4int getNumberOfTrackerSD();
  
private:
  static CCalSDList* theList;
  std::vector<nameType> caloSD;
  std::vector<nameType> trackerSD;
  
private:
  CCalSDList& operator=(CCalSDList&);   

};

#endif
