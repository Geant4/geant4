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

#ifndef G4ICRU90StoppingData_h
#define G4ICRU90StoppingData_h 1

//----------------------------------------------------------------------
//
// File name:    G4ICRU90StoppingData
//
// Description:  Data on electroninc stopping power from ICRU 90
//
// Author:       Lucas Norberto Burigo
//
// Creation date: 03.09.2018
//
// Modifications: 25.09.2018 V.Ivanchenko adopted for material sub-library
// 
//----------------------------------------------------------------------------
//
// Class Description:
//
// Data on electonic stopping powers from ICRU 90 report
//   
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "G4Material.hh"
#include "G4LPhysicsFreeVector.hh"

class G4ICRU90StoppingData 
{ 
public: 

  explicit G4ICRU90StoppingData();

  ~G4ICRU90StoppingData();

  void Initialise();

  G4double 
  GetElectronicDEDXforProton(const G4Material*, G4double kinEnergy) const;

  G4double 
  GetElectronicDEDXforAlpha(const G4Material*, G4double scaledKinEnergy) const;

  inline G4int GetIndex(const G4Material*) const;

  inline G4int GetIndex(const G4String&) const;

  inline G4double 
  GetElectronicDEDXforProton(G4int idx, G4double kinEnergy) const;

  inline G4double 
  GetElectronicDEDXforAlpha(G4int idx, G4double scaledKinEnergy) const;
 
  inline G4bool IsApplicable(const G4Material*) const;

private:

  inline G4double GetDEDX(G4LPhysicsFreeVector*, G4double e) const;

  void FillData();

  G4LPhysicsFreeVector* AddData(G4int n, const G4double* e, const G4float* dedx);

  // hide assignment operator
  G4ICRU90StoppingData & operator=
  (const  G4ICRU90StoppingData &right) = delete;
  G4ICRU90StoppingData(const G4ICRU90StoppingData&) = delete;

  static constexpr G4int nvectors = 3;
  const G4Material* materials[nvectors];
  G4LPhysicsFreeVector* sdata_proton[nvectors];
  G4LPhysicsFreeVector* sdata_alpha[nvectors];
  G4bool isInitialized;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4bool G4ICRU90StoppingData::IsApplicable(const G4Material* mat) const
{
  return (mat == materials[0] || mat == materials[1] || mat == materials[2]);
}

inline G4int G4ICRU90StoppingData::GetIndex(const G4Material* mat) const
{  
  G4int idx = -1;
  if (mat == materials[0]) { idx = 0; }
  else if(mat == materials[1]) { idx = 1; }
  else if(mat == materials[2]) { idx = 2; }
  return idx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4int G4ICRU90StoppingData::GetIndex(const G4String& nam) const
{
  G4int idx = -1;
  if (nam == materials[0]->GetName()) { idx = 0; }
  else if(nam == materials[1]->GetName()) { idx = 1; }
  else if(nam == materials[2]->GetName()) { idx = 2; }
  return idx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double 
G4ICRU90StoppingData::GetDEDX(G4LPhysicsFreeVector* data, G4double e) const
{
  G4double emin = data->Energy(0);
  return (e <= emin) ? (*data)[0]*std::sqrt(e/emin) : data->Value(e);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4ICRU90StoppingData::GetElectronicDEDXforProton(
       G4int idx, G4double kinEnergy) const
{
  return (idx < 0 || idx >= nvectors) ? 0.0 
    : GetDEDX(sdata_proton[idx], kinEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4ICRU90StoppingData::GetElectronicDEDXforAlpha(
       G4int idx, G4double scaledKinEnergy) const
{
  return (idx < 0 || idx >= nvectors) ? 0.0 
    : GetDEDX(sdata_alpha[idx], scaledKinEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
