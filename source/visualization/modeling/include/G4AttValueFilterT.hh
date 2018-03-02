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
// $Id: G4AttValueFilterT.hh 108542 2018-02-16 09:35:32Z gcosmo $
//
// Templated class for G4AttValue filters.
//
// Jane Tinslay, September 2006
//
#ifndef G4ATTVALUEFILTERT_HH
#define G4ATTVALUEFILTERT_HH

#include "G4AttValue.hh"
#include "G4VAttValueFilter.hh"
#include "G4ConversionFatalError.hh"
#include "G4ConversionUtils.hh"

namespace {
  
  // Helper classes
  template <typename T>
  class IsEqual{
  public:
    IsEqual(const T& value): fValue(value) {};
    bool operator()(const std::pair<const G4String, T>& myPair) const
    {
      return myPair.second == fValue;
    }
  private:
    T fValue;
  };
  
  template <typename T>
  class InInterval{
  public:
    InInterval(const T& value): fValue(value) {};
    bool operator()(const std::pair<const G4String, std::pair<T, T> >& myPair) const
    {
      T min = myPair.second.first;
      T max = myPair.second.second;
      return ((fValue > min || fValue == min) && (fValue < max));
    }
  private:
    T fValue;
  };

}

template <typename T, typename ConversionErrorPolicy = G4ConversionFatalError>
class G4AttValueFilterT : public ConversionErrorPolicy, public G4VAttValueFilter {
public:

  // Constructor
  G4AttValueFilterT();

  // Destructor
  virtual ~G4AttValueFilterT();

  // Filter methods
  G4bool Accept(const G4AttValue& attVal) const;
  G4bool GetValidElement(const G4AttValue& input, G4String& interval) const;

  // Print configuration
  virtual void PrintAll(std::ostream& ostr) const;
  
  // Reset 
  virtual void Reset();

  void LoadIntervalElement(const G4String& input);
  void LoadSingleValueElement(const G4String& input);

private:

  typedef std::pair<T, T> Pair;
  typedef typename std::map<G4String, Pair> IntervalMap;
  typedef std::map<G4String, T> SingleValueMap;


  // Data members  
  IntervalMap fIntervalMap;
  SingleValueMap fSingleValueMap;
  
};

template <typename T, typename ConversionErrorPolicy> 
G4AttValueFilterT<T, ConversionErrorPolicy>::G4AttValueFilterT() {}

template <typename T, typename ConversionErrorPolicy> 
G4AttValueFilterT<T, ConversionErrorPolicy>::~G4AttValueFilterT() {}

template <typename T, typename ConversionErrorPolicy>
G4bool 
G4AttValueFilterT<T, ConversionErrorPolicy>::GetValidElement(const G4AttValue& attValue, G4String& element) const 
{
  T value;
  
  G4String input = attValue.GetValue();
  if (!G4ConversionUtils::Convert(input, value)) ConversionErrorPolicy::ReportError(input, "Invalid format. Was the input data formatted correctly ?");	
  
  typename SingleValueMap::const_iterator iterValues = 
    std::find_if(fSingleValueMap.begin(), fSingleValueMap.end(), IsEqual<T>(value));

  if (iterValues != fSingleValueMap.end()) {
    element = iterValues->first;
    return true;
  }
  
  typename IntervalMap::const_iterator iterIntervals = 
    std::find_if(fIntervalMap.begin(), fIntervalMap.end(), InInterval<T>(value));

  if (iterIntervals != fIntervalMap.end()) {
    element = iterIntervals->first;
    return true;
  }

  return false;
}

template <typename T, typename ConversionErrorPolicy>
G4bool 
G4AttValueFilterT<T, ConversionErrorPolicy>::Accept(const G4AttValue& attValue) const 
{
  T value;

  G4String input = attValue.GetValue();
  if (!G4ConversionUtils::Convert(input, value)) ConversionErrorPolicy::ReportError(input, "Invalid format. Was the input data formatted correctly ?");	
  
  typename SingleValueMap::const_iterator iterValues = 
    std::find_if(fSingleValueMap.begin(), fSingleValueMap.end(), IsEqual<T>(value));

  if (iterValues != fSingleValueMap.end()) return true;
  
  typename IntervalMap::const_iterator iterIntervals = 
    std::find_if(fIntervalMap.begin(), fIntervalMap.end(), InInterval<T>(value));

  if (iterIntervals != fIntervalMap.end()) return true;

  return false;
}

template <typename T, typename ConversionErrorPolicy>
void 
G4AttValueFilterT<T, ConversionErrorPolicy>::LoadIntervalElement(const G4String& input) 
{
  T min;
  T max;

  if (!G4ConversionUtils::Convert(input, min, max)) ConversionErrorPolicy::ReportError(input, "Invalid format. Was the input data formatted correctly ?");

  std::pair<T, T> myPair(min, max);
  fIntervalMap[input] = myPair;
}

template <typename T, typename ConversionErrorPolicy>
void 
G4AttValueFilterT<T, ConversionErrorPolicy>::LoadSingleValueElement(const G4String& input) 
{
  T output;
  
  if (!G4ConversionUtils::Convert(input, output)) ConversionErrorPolicy::ReportError(input, "Invalid format. Was the input data formatted correctly ?");

  fSingleValueMap[input] = output;
}

template <typename T, typename ConversionErrorPolicy>
void 
G4AttValueFilterT<T, ConversionErrorPolicy>::PrintAll(std::ostream& ostr) const
{
  ostr<<"Printing data for filter: "<<Name()<<std::endl;

  ostr<<"Interval data:"<<std::endl;

  typename IntervalMap::const_iterator iterIntervals = fIntervalMap.begin();
 
  while (iterIntervals != fIntervalMap.end()) {
    ostr<<iterIntervals->second.first<<" : "<<iterIntervals->second.second<<std::endl;
    iterIntervals++;
  }

  ostr<<"Single value data:"<<std::endl;

  typename SingleValueMap::const_iterator iterValues = fSingleValueMap.begin();
 
  while (iterValues != fSingleValueMap.end()) {
    ostr<<iterValues->second<<std::endl;
    iterValues++;
  }
}

template <typename T, typename ConversionErrorPolicy>
void 
G4AttValueFilterT<T, ConversionErrorPolicy>::Reset() 
{
  fIntervalMap.clear();
  fSingleValueMap.clear();
}

#endif
