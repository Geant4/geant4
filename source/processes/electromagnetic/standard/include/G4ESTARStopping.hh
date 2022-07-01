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

#ifndef G4ESTARStopping_h
#define G4ESTARStopping_h 1

//---------------------------------------------------------------------------
//
// ClassName:   G4ESTARStopping
//
// Description: Data on stopping power
//
// Author:      Anton Ivantchenko 18.04.2006
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
// Data on Stopping Powers from the NIST Data Base  
// http://physics.nist.gov/PhysRefData/STAR/Text/ESTAR.html
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "globals.hh"
#include "G4PhysicsFreeVector.hh"
#include <vector>

class G4Material;

class G4ESTARStopping 
{ 
public: 

  explicit G4ESTARStopping(const G4String& datatype = "");

  ~G4ESTARStopping();

  G4int GetIndex(const G4Material*);

  G4double GetElectronicDEDX(G4int idx, G4double energy);

  inline G4double GetElectronicDEDX(const G4Material*, G4double energy);

  // hide assignment operator
  G4ESTARStopping & operator=(const  G4ESTARStopping &right) = delete;
  G4ESTARStopping(const  G4ESTARStopping&) = delete;

private:

  void Initialise();

  void AddData(const G4double* energy, const G4double* xs, G4int idx);

  const char* dirPath = nullptr;
  const G4Material* currentMaterial = nullptr;

  G4int type = 0;
  G4int matIndex = -1;
  std::vector<G4PhysicsFreeVector*> sdata;
  std::vector<G4String> name;
};

inline G4double G4ESTARStopping::GetElectronicDEDX(const G4Material* mat, 
						   G4double energy)
{
  return GetElectronicDEDX(GetIndex(mat), energy);
}

#endif
