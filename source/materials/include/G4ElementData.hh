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
// $Id: G4ElementData.hh 96794 2016-05-09 10:09:30Z gcosmo $
//
//---------------------------------------------------------------------------
//
// GEANT4 Class file
//
// Description: Data structure for cross sections, shell cross sections,
//              isotope cross sections. Control of vector size should be
//              performed in user code, no protection in this class
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
#include "G4Physics2DVector.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ElementData 
{
public:

  explicit G4ElementData();

  ~G4ElementData();

  // add cross section for the element
  void InitialiseForElement(G4int Z, G4PhysicsVector* v);

  // add 2D cross section for the element
  void InitialiseForElement(G4int Z, G4Physics2DVector* v);

  // reserve vector of components
  void InitialiseForComponent(G4int Z, G4int nComponents=0);

  // prepare vector of components
  void AddComponent(G4int Z, G4int id, G4PhysicsVector* v);

  // set name of the dataset
  void SetName(const G4String& nam);

  // get vector for the element 
  inline G4PhysicsVector* GetElementData(G4int Z);

  // get 2-D vector for the element 
  inline G4Physics2DVector* GetElement2DData(G4int Z);

  // get number of components for the element 
  inline size_t GetNumberOfComponents(G4int Z);

  // get component ID which may be number of nucleons, 
  // or shell number, or any other integer
  inline G4int GetComponentID(G4int Z, size_t idx);

  // get vector per shell or per isotope
  inline G4PhysicsVector* GetComponentDataByIndex(G4int Z, size_t idx);

  // get vector per shell or per isotope
  inline G4PhysicsVector* GetComponentDataByID(G4int Z, G4int id);

  // return cross section per element 
  // if not available return zero
  inline G4double GetValueForElement(G4int Z, G4double kinEnergy);

  // return cross section per element 
  // if not available return zero
  inline G4double GetValueForComponent(G4int Z, size_t idx, G4double kinEnergy);

private:

  // Assignment operator and copy constructor
  G4ElementData & operator=(const G4ElementData &right) = delete;
  G4ElementData(const G4ElementData&) = delete;

  G4PhysicsVector* elmData[maxNumElements]; 
  G4Physics2DVector* elm2Data[maxNumElements];
  std::vector<G4PhysicsVector*> compData[maxNumElements];
  std::vector<G4int> compID[maxNumElements];
  size_t compLength[maxNumElements];
  G4String name;
};

inline void G4ElementData::SetName(const G4String& nam)
{
  name = nam;
}

inline 
G4PhysicsVector* G4ElementData::GetElementData(G4int Z)
{
  return elmData[Z];
}

inline 
G4Physics2DVector* G4ElementData::GetElement2DData(G4int Z)
{
  return elm2Data[Z];
}

inline 
size_t G4ElementData::GetNumberOfComponents(G4int Z)
{
  return compLength[Z];
}

inline G4int G4ElementData::GetComponentID(G4int Z, size_t idx)
{
  return (compID[Z])[idx];
}

inline 
G4PhysicsVector* G4ElementData::GetComponentDataByIndex(G4int Z, size_t idx)
{
  return (compData[Z])[idx];
}

inline 
G4PhysicsVector* G4ElementData::GetComponentDataByID(G4int Z, G4int id)
{
  G4PhysicsVector* v = 0;
  for(size_t i=0; i<compLength[Z]; ++i) {
    if(id == (compID[Z])[i]) {
      v = (compData[Z])[i];
      break;
    }
  }
  return v;
}

inline 
G4double G4ElementData::GetValueForElement(G4int Z, G4double kinEnergy)
{
  return elmData[Z]->Value(kinEnergy);
}

inline G4double 
G4ElementData::GetValueForComponent(G4int Z, size_t idx, G4double kinEnergy)
{
  return ((compData[Z])[idx])->Value(kinEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
 
