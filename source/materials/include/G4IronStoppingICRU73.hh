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
// $Id: G4IronStoppingICRU73.hh,v 1.3 2008/11/02 12:22:19 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $

#ifndef G4IronStoppingICRU73_h
#define G4IronStoppingICRU73_h 1

//---------------------------------------------------------------------------
//
// ClassName:   G4IronStoppingICRU73
//
// Description: Data on stopping powers for light ions in compounds
//
// Author:      A.Ivantchenko 8.08.2008
//
// Modifications:
//
//----------------------------------------------------------------------------
//
// Class Description:
//
// Data on Stopping Powers from the ICRU73 report
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "globals.hh"
#include "G4LPhysicsFreeVector.hh"
#include <vector>

class G4IronStoppingICRU73
{
public:

  G4IronStoppingICRU73(G4bool splineFlag = true);

  ~G4IronStoppingICRU73();

  G4double GetDEDX(G4int idxMaterial, G4double kinEnergy);

  inline G4double GetDEDX(const G4String& NameMaterial, G4double kinEnergy);

  inline G4int GetMaterialIndex(const G4String& NameMaterial);

  // Function returns an unique index (>=0) for each ion-material couple (the 
  // return value is -1 if the couple is not found):
  inline G4int GetIonMaterialCoupleIndex(
                        G4int atomicNumber,            // Atomic number of ion 
                        const G4String& materialName); // Material name

  inline G4double GetDensity(G4int idx);

  inline G4String GetMaterialName(G4int idx);

  inline G4PhysicsVector* GetPhysicsVector(G4int idx);

  inline G4PhysicsVector* GetPhysicsVector(const G4String& NameMaterial);

  inline G4double GetLowerEnergyBoundary();

  inline G4double GetUpperEnergyBoundary();

private:

  void AddData(G4double* energy, G4double* stoppower, G4double factor);

  void Initialise();

  // hide assignment operator
  G4IronStoppingICRU73 & operator=(const  G4IronStoppingICRU73 &right);
  G4IronStoppingICRU73(const  G4IronStoppingICRU73&);

  G4bool spline;
  G4String MatName[16];
  G4double Density[16];

  // Lower and upper energy boundaries for dE/dx vectors:
  G4double lowerEnergyBoundary;
  G4double upperEnergyBoundary;

  std::vector<G4LPhysicsFreeVector*>  dedx;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4IronStoppingICRU73::GetDEDX(const G4String& NameMaterial, 
					      G4double kinEnergy)
{
  return GetDEDX(GetMaterialIndex(NameMaterial), kinEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4int
G4IronStoppingICRU73::GetMaterialIndex(const G4String& NameMaterial)
{
  G4int idx = -1;
  for (G4int idxMaterial=0; idxMaterial<16; idxMaterial++){
    if(MatName[idxMaterial] == NameMaterial) {
      idx = idxMaterial;
      break;
    }
  }
  return idx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4int
G4IronStoppingICRU73::GetIonMaterialCoupleIndex(G4int atomicNumber,
                                                const G4String& materialName) 
{
  G4int idx = -1;
  if(atomicNumber == 26) idx = GetMaterialIndex(materialName);
  return idx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4IronStoppingICRU73::GetDensity(G4int idxMaterial)
{
  G4double d = 0.0;
  if( idxMaterial >= 0 && idxMaterial <= 15) d = Density[idxMaterial];
  return d;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4String G4IronStoppingICRU73::GetMaterialName(G4int idxMaterial)
{
  G4String s = "";
  if( idxMaterial >= 0 && idxMaterial <= 15) s = MatName[idxMaterial];
  return s;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline 
G4PhysicsVector* G4IronStoppingICRU73::GetPhysicsVector(G4int idxMaterial)
{
  G4PhysicsVector* v = 0;
  if(idxMaterial >= 0 && idxMaterial <= 15) v = dedx[idxMaterial];
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4PhysicsVector* 
G4IronStoppingICRU73::GetPhysicsVector(const G4String& NameMaterial)
{
  return GetPhysicsVector(GetMaterialIndex(NameMaterial));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double
G4IronStoppingICRU73::GetLowerEnergyBoundary() {

  return lowerEnergyBoundary;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double
G4IronStoppingICRU73::GetUpperEnergyBoundary() {

  return upperEnergyBoundary;
}


#endif
 
