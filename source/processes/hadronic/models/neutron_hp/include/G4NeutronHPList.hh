// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPList.hh,v 1.4 1999-12-15 14:53:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPList_h
#define G4NeutronHPList_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/fstream"

class G4NeutronHPList
{
  public:
  
  G4NeutronHPList()
  {
    theData = new G4double[100]; 
    nPoints=100;
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
  
  void Init(G4std::ifstream & aDataFile, G4int nPar, G4double unit=1.);
  
  void Init(G4std::ifstream & aDataFile, G4double unit=1.);

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
