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
// Generic dimensioned type.
//
// Jane Tinslay, September 2006
//
#ifndef G4DIMENSIONEDTYPE_HH
#define G4DIMENSIONEDTYPE_HH

#include "globals.hh"
#include "G4ConversionFatalError.hh"
#include "G4String.hh"
#include "G4UnitsTable.hh"
#include <ostream>

namespace G4DimensionedTypeUtils
{
  G4bool GetUnitValue(const G4String& unit, G4double& value);
}

// Default error handling done through G4ConversionFatalError
template <typename T, typename ConversionErrorPolicy = G4ConversionFatalError>
class G4DimensionedType : public ConversionErrorPolicy {

public:

  // Constructors
  G4DimensionedType();
  G4DimensionedType(const T& value, const G4String& unit);

  // Destructor
  virtual ~G4DimensionedType();

  // Accessors

  // Raw, undimensioned value
  T RawValue() const;
  
  // Unit string
  G4String Unit() const;
  
  // Dimensioned value - rawValue*converted unit
  T DimensionedValue() const;

  // Operators
  T operator()() const;
  bool operator == (const G4DimensionedType<T>& rhs) const;
  bool operator != (const G4DimensionedType<T>& rhs) const;
  bool operator < (const G4DimensionedType<T>& rhs) const;
  bool operator > (const G4DimensionedType<T>& rhs) const;

private:

  // Data members
  T fValue;
  G4String fUnit;
  T fDimensionedValue;

};

template <typename T, typename ConversionErrorPolicy>
G4DimensionedType<T, ConversionErrorPolicy>::G4DimensionedType()
  :fValue(0)
  ,fUnit("Undefined")
  ,fDimensionedValue(0) 
{}

template <typename T, typename ConversionErrorPolicy>
G4DimensionedType<T, ConversionErrorPolicy>::G4DimensionedType(const T& value, const G4String& unit)
  :fValue(value)
  ,fUnit(unit)
{
  G4double unitValue(0);

  // Convert unit string to unit value
  if (!G4DimensionedTypeUtils::GetUnitValue(unit, unitValue)) ConversionErrorPolicy::ReportError(unit, "Invalid unit");

  fDimensionedValue = value*unitValue;
}

template <typename T, typename ConversionErrorPolicy>
G4DimensionedType<T, ConversionErrorPolicy>::~G4DimensionedType() = default;

template <typename T, typename ConversionErrorPolicy>
T
G4DimensionedType<T, ConversionErrorPolicy>::RawValue() const
{
  return fValue;
}

template <typename T, typename ConversionErrorPolicy>
G4String
G4DimensionedType<T, ConversionErrorPolicy>::Unit() const
{
  return fUnit;
}

template <typename T, typename ConversionErrorPolicy>
T
G4DimensionedType<T, ConversionErrorPolicy>::DimensionedValue() const
{
  return fDimensionedValue;
}

template <typename T, typename ConversionErrorPolicy>
T 
G4DimensionedType<T, ConversionErrorPolicy>::operator()() const 
{
  return fDimensionedValue;
}

template <typename T, typename ConversionErrorPolicy>
bool
G4DimensionedType<T, ConversionErrorPolicy>::operator == (const G4DimensionedType<T>& rhs) const 
{
  return fDimensionedValue == rhs.fDimensionedValue;
}

template <typename T, typename ConversionErrorPolicy>
bool
G4DimensionedType<T, ConversionErrorPolicy>::operator != (const G4DimensionedType<T>& rhs) const 
{
  return fDimensionedValue != rhs.fDimensionedValue;
}

template <typename T, typename ConversionErrorPolicy>
bool
G4DimensionedType<T, ConversionErrorPolicy>::operator < (const G4DimensionedType<T>& rhs) const 
{
  return fDimensionedValue < rhs.fDimensionedValue;
}

template <typename T, typename ConversionErrorPolicy>
bool
G4DimensionedType<T, ConversionErrorPolicy>::operator > (const G4DimensionedType<T>& rhs) const 
{
  return fDimensionedValue > rhs.fDimensionedValue;
}

template <typename M>
std::ostream& operator << (std::ostream& os, const G4DimensionedType<M>& obj) {
  os << obj.RawValue()<<" "<<obj.Unit();
  return os;
}

#endif
