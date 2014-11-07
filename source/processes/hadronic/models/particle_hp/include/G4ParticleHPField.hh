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
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPField_h
#define G4ParticleHPField_h 1

#include "G4ParticleHPFieldPoint.hh"
#include "G4PhysicsVector.hh"

class G4ParticleHPField
{
  public:
  
  G4ParticleHPField();
  
  ~G4ParticleHPField();
  
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
  inline G4ParticleHPFieldPoint & GetPoint(G4int i) { return theData[i]; }
  
  G4double GetY(G4double e, G4int j);

  inline G4int GetFieldLength() {return nEntries;}

  void Dump();

  private:
  
  void Check(G4int i);
  
  G4ParticleHPFieldPoint * theData;
  G4int nEntries;
  G4int nPoints;
};

#endif
