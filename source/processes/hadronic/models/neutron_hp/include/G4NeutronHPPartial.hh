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
// 070618 Comment out unused private member leaking by T. Koi
//
#ifndef G4NeutronHPPartial_h
#define G4NeutronHPPartial_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "G4InterpolationManager.hh"
#include "G4NeutronHPInterpolator.hh"

class G4NeutronHPPartial
{
  public:
  
  G4NeutronHPPartial(G4int n)
  { 
    X = new G4double[n];
    data = new G4NeutronHPVector[n];
    nData = n;
    T=0;
  }
  
  G4NeutronHPPartial(G4int n1, G4int n2)
  {
    T = new G4double[n2];
    X = new G4double[n1];
    data = new G4NeutronHPVector[n1];
    nData = std::max(n1,n2);
  }  
  
  void InitInterpolation(G4int i, std::ifstream & aDataFile)
  {
    data[i].InitInterpolation(aDataFile);
  }
  
  void InitInterpolation(std::ifstream & aDataFile)
  {
    theManager.Init(aDataFile);
  }
  
  void Init(std::ifstream & aDataFile)
  {
    G4int i;
    G4double e;
    for( i=0; i<nData; i++)
    {
      aDataFile >> e;
      e *= CLHEP::eV;
      SetX(i,e);
      InitData(i, aDataFile, CLHEP::eV); // energy and probability for gammas
    }  
  }
  
  void InitData(G4int i, std::ifstream & aDataFile, G4double unit=1.)
  {
    G4int ii;
    G4double eg, pg;
    G4int neg;
    aDataFile >> neg;
    data[i].InitInterpolation(aDataFile);
    for (ii=0; ii<neg; ii++)
    { 
      aDataFile >> eg >> pg;
      eg *= unit;
      SetX(i,ii,eg);
      SetY(i,ii,pg);
    }          
  }
  
  ~G4NeutronHPPartial()
  {
    delete [] X;
    if(T!=0) delete [] T;
    delete [] data;
  }
  inline G4int GetNumberOfEnergies() {return nData;}
  
  inline void SetX(G4int i, G4double x) {X[i]=x;}
  inline void SetT(G4int i, G4double x) {T[i]=x;}
  inline void SetX(G4int i, G4int j, G4double x) {data[i].SetX(j,x);}
  inline void SetY(G4int i, G4int j, G4double y) {data[i].SetY(j,y);}
  
  inline G4double GetX(G4int i) {return X[i];}
  inline G4double GetT(G4int i) {return T[i];}
  inline G4double GetX(G4int i, G4int j) {return data[i].GetX(j);}
  inline G4double GetY(G4int i, G4int j) {return data[i].GetY(j);}
  inline G4double GetY(G4int i, G4double e) {return data[i].GetY(e);}
  inline G4int GetNEntries(G4int i) {return data[i].GetVectorLength();}
  G4NeutronHPVector * GetY(G4double e1);
  G4double Sample(G4double x);
  
  private:
  
  G4double * X;
  G4double * T;
  G4NeutronHPVector * data;
  // TKDB
  //G4NeutronHPVector * theBuffer;
  G4int nData;
  G4InterpolationManager theManager; // interpolate between different data[i]
  G4NeutronHPInterpolator theInt;
};

#endif
