//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4NistElementBuilder.hh,v 1.5 2005/05/12 17:29:08 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $

#ifndef G4NistElementBuilder_h
#define G4NistElementBuilder_h 1

//---------------------------------------------------------------------------
//
// ClassName:   G4NistElementBuilder
//
// Description: Utility class to hold and manipulate G4Elements defined from
//              Nist data base
//
// Author:      V.Ivanchenko 21-11-2004
//
//----------------------------------------------------------------------------
//
// Class Description:
//
// Element data from the NIST DB on Atomic Weights and Isotope Compositions
// http://physics.nist.gov/PhysRefData/Compositions/index.html
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4int maxNumElements = 108;
const G4int maxAbundance   = 3500;

class G4Element;

class G4NistElementBuilder
{
public:
  G4NistElementBuilder(G4int vb);
 ~G4NistElementBuilder();
 
  G4int    GetZ           (const G4String& symb);
  G4double GetA           (G4int Z);
  G4double GetIsotopeMass (G4int Z, G4int N);
  
  G4int    GetMaxNumElements() {return maxNumElements-1;};
    
  void SetVerbose   (G4int vb) {verbose = vb;};
  void PrintElement (G4int Z);
    
  // Find or build a G4Element by atomic number  
  G4Element* FindOrBuildElement (G4int Z, G4bool buildIsotopes = true);

  // Find  or build a G4Element by symbol
  G4Element* FindOrBuildElement (const G4String& symb,
                                          G4bool buildIsotopes = true);

private:

  void Initialise();

  void AddElement(const G4String& symbol, G4int Z, G4int nc, const G4int& N,
                  const G4double& A, const G4double& sA, const G4double& W);

  // Build a G4Element from dataBase
  G4Element* BuildElement(G4int Z, G4bool buildIsotopes);

private:

  G4String   elmSymbol     [maxNumElements];
  G4double   atomicMass    [maxNumElements];
  G4int      nIsotopes     [maxNumElements];
  G4int      nFirstIsotope [maxNumElements];
  G4int      idxIsotopes   [maxNumElements];

  G4double   massIsotopes   [maxAbundance];
  G4double   sigMass        [maxAbundance];
  G4double   relAbundance   [maxAbundance];

  G4int      index;
  G4int      verbose;
  G4bool     first;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4int G4NistElementBuilder::GetZ(const G4String& name)
{
  G4int Z = maxNumElements;
  do {Z--;} while( Z>0 && elmSymbol[Z] != name);
  return Z;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4NistElementBuilder::GetA(G4int Z)
{
  G4double a = 0.0;
  if(Z>0 && Z<maxNumElements) a = atomicMass[Z];
  return a;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4NistElementBuilder::GetIsotopeMass(G4int Z, G4int N)
{
  G4int i = N - nFirstIsotope[Z];
  if(i < 0) i = 0;
  else if(i >= nIsotopes[Z]) i = nIsotopes[Z] -1;
  i += idxIsotopes[Z];
  return massIsotopes[i];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif
