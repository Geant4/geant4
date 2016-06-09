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
//
// $Id$
//
#ifndef G4NeutronHPLegendreStore_h
#define G4NeutronHPLegendreStore_h 1

#include "G4NeutronHPLegendreTable.hh"
#include "G4InterpolationManager.hh"
#include "G4ios.hh"
#include <fstream>

class G4NeutronHPLegendreStore
{
  public:
  
  G4NeutronHPLegendreStore(G4int n)
  {
    theCoeff = new G4NeutronHPLegendreTable[n];
    nEnergy = n;
  }
  
  ~G4NeutronHPLegendreStore()
  {
    delete [] theCoeff;
  }
  
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
    if(i>nEnergy) throw G4HadronicException(__FILE__, __LINE__, "LegendreTableIndex out of range");
    theCoeff[i] = *theTable;
// not here -- see G4NeutronHPPhotonDist.cc line 275
//    delete theTable;
  }
  
  inline G4double GetCoeff(G4int i, G4int l) {return theCoeff[i].GetCoeff(l);}
  inline G4double GetEnergy(G4int i){return theCoeff[i].GetEnergy();}
  inline G4double GetTemperature(G4int i){return theCoeff[i].GetTemperature();}
  inline G4int GetNumberOfPoly(G4int i) {return theCoeff[i].GetNumberOfPoly();}
  
  G4double SampleDiscreteTwoBody (G4double anEnergy);
  G4double SampleElastic (G4double anEnergy);
  G4double Sample (G4double energy);
  G4double SampleMax (G4double energy);
  G4double Integrate(G4int k, G4double costh);
  
  void InitInterpolation(std::ifstream & aDataFile)
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
