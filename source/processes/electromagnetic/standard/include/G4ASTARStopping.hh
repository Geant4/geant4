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
// $Id$

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "G4LPhysicsFreeVector.hh"
#include <vector>

class G4Material;

class G4ASTARStopping 
{ 
public: 

  G4ASTARStopping();

  ~G4ASTARStopping();

  G4int GetIndex(const G4Material*);
  G4double GetEffectiveZ(G4int idx);

  G4double GetElectronicDEDX(G4int idx, G4double energy);
  G4double GetElectronicDEDX(const G4Material*, G4double energy);

private:

  void Initialise();

  void AddData(G4double* e, G4double* s, G4int idx);

  // hide assignment operator
  G4ASTARStopping & operator=(const  G4ASTARStopping &right);
  G4ASTARStopping(const  G4ASTARStopping&);

  G4int matIndex;
  const G4Material* currentMaterial;
  G4double emin;
  std::vector<G4String> name;
  std::vector<G4double> effZ; 
  std::vector<G4LPhysicsFreeVector*> sdata;
};

inline G4double G4ASTARStopping::GetElectronicDEDX(const G4Material* mat, 
						   G4double energy)
{
  return GetElectronicDEDX(GetIndex(mat), energy);
}

inline G4double G4ASTARStopping::GetEffectiveZ(G4int idx)
{
  return effZ[idx];
} 

#endif
