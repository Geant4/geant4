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
//
#ifndef G4NeutronHPList_h
#define G4NeutronHPList_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>

class G4NeutronHPList
{
  public:
  
  G4NeutronHPList()
  {
    theData = new G4double[2]; 
    nPoints=2;
    nEntries=0;
  }
  
  ~G4NeutronHPList()
  {
    delete [] theData;
  }
  
  inline void SetValue(G4int i, G4double y) 
  { 
    Check(i);
    theData[i]=y;
  }
  G4double GetValue(G4int i);
  
  inline G4int GetListLength() {return nEntries;}

  void Dump();
  
  void Init(std::ifstream & aDataFile, G4int nPar, G4double unit=1.);
  
  void Init(std::ifstream & aDataFile, G4double unit=1.);

  inline void SetLabel(G4double aLabel) { theLabel = aLabel; }
  
  inline G4double GetLabel() { return theLabel; }

  private:
  
  void Check(G4int i);
 
  G4double theLabel;

  G4double * theData;
  G4int nEntries;
  G4int nPoints;
};

#endif
