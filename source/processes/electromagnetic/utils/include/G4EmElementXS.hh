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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4EmElementXS
//
// Author:        V. Ivanchenko 
// 
// Creation date: 22 August 2024
//
// Modifications: 
//
// Class Description: 
//
// Class provide access to the cross section data stored in the G4LEDATA
// data structure in form of G4PhysicsVector. For that data are uploaded
// to the G4ElementData via G4ElementDataRegistry allowing sharing of
// data between threads.
//
// When retrieve data, this class set mutex lock.
//
// Class Description: End 

// -------------------------------------------------------------------
//

#ifndef G4EmElementXS_h
#define G4EmElementXS_h 1

#include "globals.hh"

class G4ElementData;
class G4EmParameters;
class G4PhysicsVector;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4EmElementXS
{
public:

  explicit G4EmElementXS(G4int zmin, G4int zmax, const G4String& name,
			 const G4String& subname);

  ~G4EmElementXS() = default;

  G4PhysicsVector* Retrieve(G4int Z) const;

  G4double GetXS(G4int Z, G4double ekin) const;

  // hide assignment operator 
  G4EmElementXS & operator=(const G4EmElementXS &right) = delete;
  G4EmElementXS(const G4EmElementXS&) = delete;

private:

  G4int Zmin;
  G4int Zmax;
  G4ElementData* fData;
  G4EmParameters* fParameters;
  const G4String fSubName;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

