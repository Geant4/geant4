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

#ifndef G4ASTARStopping_h
#define G4ASTARStopping_h 1

//---------------------------------------------------------------------------
//
// ClassName:   G4ASTARStopping
//
// Description: Data on stopping power
//
// Author:      Anton Ivantchenko 21.04.2006
//
// Organisation:        QinetiQ Ltd, UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            CSMAN-5288
//
// Modifications:
// 
//----------------------------------------------------------------------------
//
// Class Description:
//
// Data on Stopping Powers of the He atoms from the NIST Data Base  
// http://physics.nist.gov/PhysRefData/STAR/Text/ASTAR.html
//
// Current ASTAR database is available via url:
// http://physics.nist.gov/PhysRefData/Star/Text/ASTAR.html

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "G4Material.hh"
#include "G4LPhysicsFreeVector.hh"
#include <vector>

class G4ASTARStopping 
{ 
public: 

  explicit G4ASTARStopping();

  ~G4ASTARStopping();

  void Initialise();

  inline G4int GetIndex(const G4Material*) const;

  inline G4int GetIndex(const G4String&) const;

  inline G4double GetElectronicDEDX(G4int idx, G4double energy) const;

  inline G4double GetElectronicDEDX(const G4Material*, G4double energy) const;

private:

  void AddData(const G4float* s, const G4Material*);

  void FindData(G4int idx, const G4Material*);

  void PrintWarning(G4int idx) const;

  // hide assignment operator
  G4ASTARStopping & operator=(const  G4ASTARStopping &right) = delete;
  G4ASTARStopping(const  G4ASTARStopping&) = delete;

  size_t nvectors;
  G4double emin;
  std::vector<const G4Material*> materials;
  std::vector<G4LPhysicsFreeVector*> sdata;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4int G4ASTARStopping:: GetIndex (const G4Material* mat) const
{  
  G4int idx = -1;
  for (size_t i=0; i<nvectors; ++i){
    if (mat == materials[i]){
      idx = i;
      break;
    }
  }
  return idx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4int G4ASTARStopping::GetIndex(const G4String& nam) const
{
  G4int idx = -1;
  for (size_t i=0; i<nvectors; ++i){
    if (nam == materials[i]->GetName()){
      idx = i;
      break;
    }
  }
  return idx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double 
G4ASTARStopping::GetElectronicDEDX(G4int idx, G4double energy) const
{
  G4double res = 0.0;
  if (idx<0 || idx>= G4int(nvectors)) { PrintWarning(idx); }
  if(energy < emin) { res = (*(sdata[idx]))[0]*std::sqrt(energy/emin); } 
  else              { res = sdata[idx]->Value(energy); }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double 
G4ASTARStopping::GetElectronicDEDX(const G4Material* mat, G4double energy) const
{
  return GetElectronicDEDX(GetIndex(mat), energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
