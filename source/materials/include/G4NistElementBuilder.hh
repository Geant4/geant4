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

#ifndef G4NistElementBuilder_h
#define G4NistElementBuilder_h 1

//---------------------------------------------------------------------------
//
// ClassName:   G4NistElementBuilder
//
// Description: Utility class to hold and manipulate G4Elements defined from
//              Nist data base
//
// Author:      V.Ivanchenko 21.11.2004
//
// Modifications:
// 27.02.06 V.Ivanchenko Return m=0 if Z&N combination is out of NIST  
// 27.02.06 V.Ivanchenko add GetAtomicMassAmu 
// 17.10.06 V.Ivanchenko add GetAtomicMass and GetNistElementNames methods
// 02.05.07 V.Ivanchenko add GetNistFirstIsotopeN and GetNumberOfNistIsotopes 
// 06.08.08 V.Ivanchenko add binding energy parameterisation and use isotope 
//                       mass in G4 units 
//
//----------------------------------------------------------------------------
//
// Class Description:
//
// Element data from the NIST DB on Atomic Weights and Isotope Compositions
// http://physics.nist.gov/PhysRefData/Compositions/index.html
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <vector>
#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"
#include "G4Element.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4int maxNumElements = 108;
const G4int maxAbundance   = 3500;

class G4NistElementBuilder
{
public:

  explicit G4NistElementBuilder(G4int vb);
  ~G4NistElementBuilder() = default;

  // Find or build a G4Element by atomic number
  inline G4Element* FindElement (G4int Z) const;
  G4Element* FindOrBuildElement (G4int Z, G4bool buildIsotopes = true);

  // Find  or build a G4Element by symbol
  G4Element* FindOrBuildElement (const G4String& symb,	
			         G4bool buildIsotopes = true);
  // print element information
  void PrintElement (G4int Z) const;

  // Access to the vector of Geant4 predefined element names 
  const std::vector<G4String>& GetElementNames() const;

  // Get atomic number by element symbol
  G4int GetZ(const G4String& symb) const;

  // Get atomic weight in atomic units by element symbol
  G4double GetAtomicMassAmu(const G4String& symb) const;

  // Get atomic weight in atomic units - mean mass in units of amu of an atom 
  // with electron shell for the natural isotope composition 
  inline G4double GetAtomicMassAmu(G4int Z) const;

  // Get mass of isotope without electron shell in Geant4 energy units 
  inline G4double GetIsotopeMass(G4int Z, G4int N) const;

  // Get mass in Geant4 energy units of an atom of a particular isotope 
  // with the electron shell  
  inline G4double GetAtomicMass(G4int Z, G4int N) const;

  // Get total ionisation energy of an atom
  inline G4double GetTotalElectronBindingEnergy(G4int Z) const;

  // Get natural isotope abundance
  inline G4double GetIsotopeAbundance (G4int Z, G4int N) const;

  // Get N for the first natural isotope
  inline G4int GetNistFirstIsotopeN(G4int Z) const;

  // Get number of natural isotopes
  inline G4int GetNumberOfNistIsotopes(G4int Z) const;

  // Get max Z in the Geant4 element database
  inline G4int GetMaxNumElements() const; 

  inline void SetVerbose(G4int);

private:

  void Initialise();

  // Add element parameters to internal G4 database: 
  // Z - atomic number, N - number of nucleons, A - atomic mass (amu),
  // sigmaA - accuracy of mass in last digits, W - natural abundances (percent) 
  void AddElement(const G4String& symbol, G4int Z, G4int NumberOfIsotopes,
                  const G4int& N, const G4double& A, const G4double& sigmaA, 
		  const G4double& W);

  // Build a G4Element from the G4 dataBase
  G4Element* BuildElement(G4int Z);

private:

  G4String   elmSymbol     [maxNumElements];
  G4double   atomicMass    [maxNumElements];  // amu
  G4double   bindingEnergy [maxNumElements];
  G4int      nIsotopes     [maxNumElements];
  G4int      nFirstIsotope [maxNumElements];
  G4int      idxIsotopes   [maxNumElements];

  G4int      elmIndex      [maxNumElements];

  G4double   massIsotopes  [maxAbundance];    // G4 units
  G4double   sigMass       [maxAbundance];    // G4 units
  G4double   relAbundance  [maxAbundance];

  G4int      index;
  G4int      verbose;

  std::vector<G4String>    elmNames;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4NistElementBuilder::GetAtomicMassAmu(G4int Z) const
{
  return (Z>0 && Z<maxNumElements) ? atomicMass[Z] : 0.0; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4NistElementBuilder::GetIsotopeMass(G4int Z, G4int N) const
{
  G4double mass = 0.0;
  if(Z > 0 && Z < maxNumElements) {
    G4int i = N - nFirstIsotope[Z];
    if(i >= 0 && i <nIsotopes[Z]) {mass = massIsotopes[i + idxIsotopes[Z]];}
  }
  return mass;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4NistElementBuilder::GetAtomicMass(G4int Z, G4int N) const
{
  G4double mass = 0.0;
  if(Z > 0 && Z < maxNumElements) {
    G4int i = N - nFirstIsotope[Z];
    if(i >= 0 && i <nIsotopes[Z]) {
      mass = massIsotopes[i + idxIsotopes[Z]] + 
	Z*CLHEP::electron_mass_c2 - bindingEnergy[Z]; 
    }
  }
  return mass;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline 
G4double G4NistElementBuilder::GetTotalElectronBindingEnergy(G4int Z) const
{
  return (Z > 0 && Z < maxNumElements) ? bindingEnergy[Z] : 0.0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline 
G4double G4NistElementBuilder::GetIsotopeAbundance(G4int Z, G4int N) const
{
  G4double x = 0.0;
  if(Z > 0 && Z < maxNumElements) {
    G4int i = N - nFirstIsotope[Z];
    if(i >= 0 && i <nIsotopes[Z]) { x = relAbundance[i + idxIsotopes[Z]]; }
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4int G4NistElementBuilder::GetNistFirstIsotopeN(G4int Z) const
{
  return (Z > 0 && Z < maxNumElements) ? nFirstIsotope[Z] : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4int G4NistElementBuilder::GetNumberOfNistIsotopes(G4int Z) const
{
  return (Z > 0 && Z < maxNumElements) ? nIsotopes[Z] : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline 
const std::vector<G4String>& G4NistElementBuilder::GetElementNames() const
{
  return elmNames;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4int G4NistElementBuilder::GetMaxNumElements() const
{
  return maxNumElements-1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4NistElementBuilder::SetVerbose(G4int val) 
{
  verbose = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4Element* G4NistElementBuilder::FindElement(G4int Z) const
{
  const G4ElementTable* theElementTable = G4Element::GetElementTable();
  return (Z > 0 && Z < maxNumElements && elmIndex[Z] >= 0) ?
    (*theElementTable)[elmIndex[Z]] : nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif
