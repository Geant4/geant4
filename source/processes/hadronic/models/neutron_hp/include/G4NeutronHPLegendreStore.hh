// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPLegendreStore.hh,v 1.2 1999-06-29 18:44:04 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPLegendreStore_h
#define G4NeutronHPLegendreStore_h 1

#include "G4NeutronHPLegendreTable.hh"
#include "G4InterpolationManager.hh"
#include "G4ios.hh"
#include <fstream.h>

class G4NeutronHPLegendreStore
{
  public:
  
  G4NeutronHPLegendreStore(G4int n);
  
  ~G4NeutronHPLegendreStore();
  inline void Init(G4int i, G4double e, G4int n)
  {
    theCoeff[i].Init(e, n);
  }
  inline void SetNPoints(G4int n) { nEnergy = n; }
  inline void SetEnergy(G4int i, G4double energy) { theCoeff[i].SetEnergy(energy); }
  inline void SetTemperature(G4int i, G4double temp) { theCoeff[i].SetTemperature(temp); }
  inline void SetCoeff(G4int i, G4int l, G4double coeff) {theCoeff[i].SetCoeff(l, coeff); }
  inline void SetCoeff(G4int i, G4NeutronHPLegendreTable * theTable)
  {
    if(i>nEnergy) G4Exception("LegendreTableIndex out of range");
    theCoeff[i] = *theTable;
// not here -- see G4NeutronHPPhotonDist.cc line 275
//    delete theTable;
  }
  
  inline G4double GetCoeff(G4int i, G4int l) {return theCoeff[i].GetCoeff(l);}
  inline G4double GetEnergy(G4int i){return theCoeff[i].GetEnergy();}
  inline G4double GetTemperature(G4int i){return theCoeff[i].GetTemperature();}
  inline G4int GetNumberOfPoly(G4int i) {return theCoeff[i].GetNumberOfPoly();}
  
  G4double SampleElastic (G4double anEnergy);
  G4double Sample (G4double energy);
  G4double SampleMax (G4double energy);
  G4double Integrate(G4int k, G4double costh);
  
  void InitInterpolation(ifstream & aDataFile)
  {
    theManager.Init(aDataFile);
  }
  
  void SetManager(G4InterpolationManager & aManager)
  {
    theManager = aManager;
  }

  private:
  
  G4int nEnergy;
  G4NeutronHPLegendreTable * theCoeff;
  G4InterpolationManager theManager; // interpolate between different Tables
};
#endif
