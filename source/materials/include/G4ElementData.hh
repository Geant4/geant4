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

//---------------------------------------------------------------------------
//
// GEANT4 Class file
//
// Description: Data structure for cross sections, shell cross sections,
//              isotope cross sections. Data access via integer variable
//              Z (atomic number in majority of applications), which may
//              be in the interval 0 <= Z < length. For isotope like
//              data a second parameter idx or data ID code are used.
//              In most cases ID = A - atomic weight number.
//              There are run time const methods, in which input is not checked
//              assuming responsibility of consumer code. Another run time
//              check input and may throwgh a fatal exception
//
// Author:      V.Ivanchenko 10.03.2011
//
// Modifications: 30.09.2023 Extended functionality, data size defined in constructor
//
//----------------------------------------------------------------------------
//

#ifndef G4ElementData_h
#define G4ElementData_h 1

#include "G4Physics2DVector.hh"
#include "G4PhysicsVector.hh"
#include "globals.hh"

#include <vector>

class G4ElementData
{
 public:
  explicit G4ElementData(G4int length = 99);

  ~G4ElementData();

  // Assignment operator and copy constructor
  G4ElementData& operator=(const G4ElementData& right) = delete;
  G4ElementData(const G4ElementData&) = delete;

  // add cross section for the element
  void InitialiseForElement(G4int Z, G4PhysicsVector* v);

  // add 2D cross section for the element
  void InitialiseForElement(G4int Z, G4Physics2DVector* v);

  // reserve vector of components
  void InitialiseForComponent(G4int Z, G4int nComponents = 0);

  // reserve vector of 2D components
  void InitialiseFor2DComponent(G4int Z, G4int nComponents = 0);

  // prepare vector of components
  void AddComponent(G4int Z, G4int id, G4PhysicsVector* v);

  // prepare vector of 2D components
  void Add2DComponent(G4int Z, G4int id, G4Physics2DVector* v);

  // set name of the dataset (optional)
  inline void SetName(const G4String& nam);

  //--------------------------------------------------------------
  // run time const methods - no check on validity of input
  // it is a responsibility of the consume code to check the input
  //--------------------------------------------------------------

  // get name of the dataset
  inline const G4String& GetName() const;

  // get vector for the element
  inline G4PhysicsVector* GetElementData(G4int Z) const;

  // get 2-D vector for the element
  inline G4Physics2DVector* GetElement2DData(G4int Z) const;

  // get vector per shell or per isotope
  inline G4PhysicsVector* GetComponentDataByID(G4int Z, G4int id) const;

  // get vector per shell or per isotope
  inline G4Physics2DVector* Get2DComponentDataByID(G4int Z, G4int id) const;

  // return cross section per element
  inline G4double GetValueForElement(G4int Z, G4double kinEnergy) const;

  //--------------------------------------------------------------
  // run time const methods with input parameters control
  //--------------------------------------------------------------

  // get number of components for the element
  inline std::size_t GetNumberOfComponents(G4int Z) const;

  // get number of 2D components for the element
  inline std::size_t GetNumberOf2DComponents(G4int Z) const;

  // get component ID which may be number of nucleons,
  // or shell number, or any other integer
  inline G4int GetComponentID(G4int Z, std::size_t idx) const;

  // get vector per shell or per isotope
  inline G4PhysicsVector*
  GetComponentDataByIndex(G4int Z, std::size_t idx) const;

  // get vector per shell or per isotope
  inline G4Physics2DVector*
  Get2DComponentDataByIndex(G4int Z, std::size_t idx) const;

  // return cross section per element
  // if not available return zero
  inline G4double
  GetValueForComponent(G4int Z, std::size_t idx, G4double kinEnergy) const;

 private:

  void DataError(G4int Z, const G4String&);

  const G4int maxNumElm;

  std::vector<G4PhysicsVector*> elmData;
  std::vector<G4Physics2DVector*> elm2Data;
  std::vector<std::vector<std::pair<G4int, G4PhysicsVector*> >* > compData;
  std::vector<std::vector<std::pair<G4int, G4Physics2DVector*> >* > comp2D;

  G4String name{""};
};

//--------------------------------------------------------------
// run time const methods without check on validity of input
//--------------------------------------------------------------

inline void G4ElementData::SetName(const G4String& nam)
{
  name = nam;
}

inline const G4String& G4ElementData::GetName() const
{
  return name;
}

inline G4PhysicsVector* G4ElementData::GetElementData(G4int Z) const
{
  return elmData[Z];
}

inline G4Physics2DVector* G4ElementData::GetElement2DData(G4int Z) const
{
  return elm2Data[Z];
}

inline G4PhysicsVector*
G4ElementData::GetComponentDataByID(G4int Z, G4int id) const
{
  G4PhysicsVector* v = nullptr;
  for (auto const & p : *(compData[Z])) {
    if (id == p.first) {
      v = p.second;
      break;
    }
  }
  return v;
}

inline G4Physics2DVector*
G4ElementData::Get2DComponentDataByID(G4int Z, G4int id) const
{
  G4Physics2DVector* v = nullptr;
  for (auto const & p : *(comp2D[Z])) {
    if (id == p.first) {
      v = p.second;
      break;
    }
  }
  return v;
}

inline G4double
G4ElementData::GetValueForElement(G4int Z, G4double kinEnergy) const
{
  return elmData[Z]->Value(kinEnergy);
}

//--------------------------------------------------------------
// run time const methods with check on validity of input
//--------------------------------------------------------------

inline std::size_t G4ElementData::GetNumberOfComponents(G4int Z) const
{
  return (nullptr != compData[Z]) ? compData[Z]->size() : 0;
}

inline std::size_t G4ElementData::GetNumberOf2DComponents(G4int Z) const
{
  return (nullptr != comp2D[Z]) ? comp2D[Z]->size() : 0;
}

inline G4int G4ElementData::GetComponentID(G4int Z, std::size_t idx) const
{
  return (idx < GetNumberOfComponents(Z)) ? (*(compData[Z]))[idx].first : 0;
}

inline G4PhysicsVector*
G4ElementData::GetComponentDataByIndex(G4int Z, std::size_t idx) const
{
  return
    (idx < GetNumberOfComponents(Z)) ? (*(compData[Z]))[idx].second : nullptr;
}

inline G4Physics2DVector*
G4ElementData::Get2DComponentDataByIndex(G4int Z, std::size_t idx) const
{
  return
    (idx < GetNumberOf2DComponents(Z)) ? (*(comp2D[Z]))[idx].second : nullptr;
}

inline G4double
G4ElementData::GetValueForComponent(G4int Z, std::size_t idx, G4double e) const
{
  return (idx < GetNumberOfComponents(Z)) ?
	  (*(compData[Z]))[idx].second->Value(e) : 0.0;
}

#endif
