// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPField.hh,v 1.2 1999-06-29 18:43:58 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPField_h
#define G4NeutronHPField_h 1

#include "G4NeutronHPFieldPoint.hh"
#include "G4PhysicsVector.hh"

class G4NeutronHPField
{
  public:
  
  G4NeutronHPField();
  ~G4NeutronHPField();
  
  inline void InitY(G4int i, G4int n)
  {
    Check(i);
    theData[i].InitY(n);
  }
  inline void SetData(G4int i, G4double x, G4int j, G4double y) 
  { 
    Check(i);
    theData[i].SetData(x, j, y);
  }
  inline void SetEnergy(G4int i, G4double e)
  {
    Check(i);
    theData[i].SetX(e);
  }
  inline void SetX(G4int i, G4double e)
  {
    Check(i);
    theData[i].SetX(e);
  }
  inline void SetY(G4int i, G4int j, G4double x)
  {
    Check(i);
    theData[i].SetY(j, x);
  }
  inline G4double GetEnergy(G4int i) { return theData[i].GetX(); }
  inline G4double GetX(G4int i) { return theData[i].GetX(); }
  inline G4double GetY(G4int i, G4int j) { return theData[i].GetY(j); }
  inline G4NeutronHPFieldPoint & GetPoint(G4int i) { return theData[i]; }
  
  G4double GetY(G4double e, G4int j);

  inline G4int GetFieldLength() {return nEntries;}

  void Dump();

  private:
  
  void Check(G4int i);
  
  G4NeutronHPFieldPoint * theData;
  G4int nEntries;
  G4int nPoints;
};

#endif
