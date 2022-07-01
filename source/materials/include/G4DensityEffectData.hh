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
//---------------------------------------------------------------------------
//
// GEANT4 Class file
//
// Description: Data on density effect
//
// Authors:   A.Bagulya, A.Ivanchenko 28.10.2009
//
//----------------------------------------------------------------------------
//
//  Data are taken from:  
//  R.M. Sternheimer et al. Density Effect for the Ionization Loss of Charged 
//  Particles in Various Substances. Atom. Data Nucl. Data Tabl. 30 (1984) 261-271. 

#ifndef DensityEffectData_h
#define DensityEffectData_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "G4Material.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4int NDENSDATA  = 278;
const G4int NDENSARRAY = 10;
const G4int NDENSELEM  = 98;

class G4DensityEffectData 
{
public:

  explicit G4DensityEffectData();

  ~G4DensityEffectData() = default;

  // return index by Z, -1 if material is not in the table 
  G4int GetElementIndex(G4int Z, G4State st = kStateUndefined) const;

  // return index by material name, -1 if material is not in the table 
  G4int GetIndex(const G4String& matName) const;

  // printout data for material
  void PrintData(const G4String& matName) const;

  // printout all data
  void DumpData() const;

  // Access to the data via index
  inline G4double GetPlasmaEnergy(G4int idx) const; 
  inline G4double GetAdjustmentFactor(G4int idx) const; 
  inline G4double GetCdensity(G4int idx) const; 
  inline G4double GetX0density(G4int idx) const; 
  inline G4double GetX1density(G4int idx) const; 
  inline G4double GetAdensity(G4int idx) const; 
  inline G4double GetMdensity(G4int idx) const; 
  inline G4double GetDelta0density(G4int idx) const; 
  inline G4double GetErrorDensity(G4int idx) const; 
  inline G4double GetMeanIonisationPotential(G4int idx) const; 

  // Assignment operator and copy constructor
  G4DensityEffectData & operator=(const G4DensityEffectData &right) = delete;
  G4DensityEffectData(const G4DensityEffectData&) = delete;

private:

  void Initialize();

  void AddMaterial(G4double* val, const G4String& matName);

  G4double data[NDENSDATA][NDENSARRAY];
  std::vector<G4String> names;

  // indexes defined only for pure materials 
  G4int indexZ[NDENSELEM];
  G4State state[NDENSELEM];

  G4int index;
};

inline G4double G4DensityEffectData::GetPlasmaEnergy(G4int idx) const
{
  return (idx >= 0 && idx < NDENSDATA) ? data[idx][0] : DBL_MAX; 
} 

inline G4double G4DensityEffectData::GetAdjustmentFactor(G4int idx) const
{
  return (idx >= 0 && idx < NDENSDATA) ? data[idx][1] : DBL_MAX; 
} 

inline G4double G4DensityEffectData::GetCdensity(G4int idx) const
{
  return (idx >= 0 && idx < NDENSDATA) ? data[idx][2] : DBL_MAX;
} 

inline G4double G4DensityEffectData::GetX0density(G4int idx) const
{
  return (idx >= 0 && idx < NDENSDATA) ?  data[idx][3] : DBL_MAX;
} 

inline G4double G4DensityEffectData::GetX1density(G4int idx) const
{
  return (idx >= 0 && idx < NDENSDATA) ? data[idx][4] : DBL_MAX; 
} 

inline G4double G4DensityEffectData::GetAdensity(G4int idx) const
{
  return (idx >= 0 && idx < NDENSDATA) ?  data[idx][5] : DBL_MAX;
}
 
inline G4double G4DensityEffectData::GetMdensity(G4int idx) const
{
  return (idx >= 0 && idx < NDENSDATA) ?  data[idx][6] : DBL_MAX;
} 

inline G4double G4DensityEffectData::GetDelta0density(G4int idx) const
{
  return (idx >= 0 && idx < NDENSDATA) ?  data[idx][7] : DBL_MAX;
} 

inline G4double G4DensityEffectData::GetErrorDensity(G4int idx) const
{
  return (idx >= 0 && idx < NDENSDATA) ?  data[idx][8] : DBL_MAX;
} 

inline G4double G4DensityEffectData::GetMeanIonisationPotential(G4int idx) const
{
  return (idx >= 0 && idx < NDENSDATA) ?  data[idx][9] : DBL_MAX;
} 


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
 
