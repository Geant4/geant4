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
// $Id: G4ElementData.hh,v 1.10 2010-05-15 15:37:33 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// GEANT4 Class file
//
// Description: Data structure for cross sections, shell cross sections
//              isotope cross sections
//
// Author:      V.Ivanchenko 10.03.2011
//
// Modifications:
//
//----------------------------------------------------------------------------
//

#ifndef ElementData_h
#define ElementData_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "G4NistElementBuilder.hh"
#include "G4PhysicsVector.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ElementData 
{
public:

  G4ElementData(const G4String& nam = "");

  ~G4ElementData();

  // add cross section for the element
  void InitialiseForElement(G4int Z, G4PhysicsVector* v);

  // prepare vector of components
  void InitialiseForComponent(G4int Z, G4int nComponents);

  // prepare vector of components
  void AddComponent(G4int Z, G4int id, size_t idx, G4PhysicsVector* v);

  // check if cross section for the element is available
  inline G4bool IsInitializedForElement(G4int Z);

  // check if the cross section per shell or per isotope
  // is available
  inline G4bool IsInitializedForComponent(G4int Z, size_t idx);

  // return cross section per element 
  // if not available return zero
  inline G4double GetValueForElement(G4int Z, G4double kinEnergy);

  // return cross section per element 
  // if not available return zero
  inline G4double GetValueForComponent(G4int Z, size_t idx, G4double kinEnergy);

  inline G4int GetComponentID(G4int Z, size_t idx);

private:

  // Assignment operator and copy constructor
  G4ElementData & operator=(const G4ElementData &right);
  G4ElementData(const G4ElementData&);

  G4PhysicsVector* elmData[maxNumElements];
  std::vector<G4PhysicsVector*> compData[maxNumElements];
  std::vector<G4int> compID[maxNumElements];
  size_t compLength[maxNumElements];

  G4String name;
};

inline 
G4bool G4ElementData::IsInitializedForElement(G4int Z)
{
  return G4bool(elmData[Z]);
}

inline 
G4bool G4ElementData::IsInitializedForComponent(G4int Z, size_t idx)
{
  return (idx < compLength[Z]);
}

inline 
G4double G4ElementData::GetValueForElement(G4int Z, G4double kinEnergy)
{
  return elmData[Z]->Value(kinEnergy);
}

inline 
G4double G4ElementData::GetValueForComponent(G4int Z, size_t idx, G4double kinEnergy)
{
  return ((compData[Z])[idx])->Value(kinEnergy);
}

inline G4int G4ElementData::GetComponentID(G4int Z, size_t idx)
{
  G4int i = -1;
  if((compID[Z]).size() < idx) { i = (compID[Z])[idx]; }
  return i;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
 
