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
// Generic attribute filter.
//
// Jane Tinslay, May 2006
//
#ifndef G4ATTRIBUTEFILTERT_HH
#define G4ATTRIBUTEFILTERT_HH

#include "G4AttDef.hh"
#include "G4AttFilterUtils.hh"
#include "G4AttUtils.hh"
#include "G4AttValue.hh"
#include "G4SmartFilter.hh"
#include "G4VAttValueFilter.hh"
#include <vector>

template <typename T>
class G4AttributeFilterT : public G4SmartFilter<T> {

public:
 
  // Construct with filter name
  G4AttributeFilterT(const G4String& name = "Unspecified");
  
  // Destructor
  virtual ~G4AttributeFilterT();

  // Evaluate
  virtual bool Evaluate(const T&) const;

  // Print configuration
  virtual void Print(std::ostream& ostr) const;

  // Clear filter
  virtual void Clear();

  // Configuration functions
  void Set(const G4String& name);
  void AddInterval(const G4String&);
  void AddValue(const G4String&);

private:

  enum Config {Interval, SingleValue};

  typedef std::pair<G4String, Config> Pair;
  typedef std::vector<Pair> ConfigVect;

  // Data members
  G4String fAttName;
  ConfigVect fConfigVect;

  // Caching
  mutable G4bool fFirst;
  mutable G4VAttValueFilter* filter;

};

template <typename T>
G4AttributeFilterT<T>::G4AttributeFilterT(const G4String& name)
  :G4SmartFilter<T>(name)
  ,fAttName("")
  ,fFirst(true)
  ,filter(0)
{}

template <typename T>
G4AttributeFilterT<T>::~G4AttributeFilterT() 
{
  delete filter;
}

template <typename T>
G4bool
G4AttributeFilterT<T>::Evaluate(const T& object) const
{
  // Return true (i.e., do not filter out) if attribute name has not yet been set.
  if (fAttName.empty()) return true;
  
  // ...or required attribute value has not yet been set
  if (fConfigVect.size() == 0) return true;

  if (fFirst) {

    fFirst = false;

    // Get attribute definition
    G4AttDef attDef;
    
    // Expect definition to exist    
    if (!G4AttUtils::ExtractAttDef(object, fAttName, attDef)) {
      static G4bool warnedUnableToExtract = false;
      if (!warnedUnableToExtract) {
	G4ExceptionDescription ed;
	ed <<"Unable to extract attribute definition named "<<fAttName<<'\n'
        << "Available attributes:\n"
        << *object.GetAttDefs();
	G4Exception
	  ("G4AttributeFilterT::Evaluate", "modeling0102", JustWarning, ed, "Invalid attribute definition");
	warnedUnableToExtract = true;
      }
      return false;
    }
    
    // Get new G4AttValue filter
    filter = G4AttFilterUtils::GetNewFilter(attDef);

    // Load both interval and single valued data.
    typename ConfigVect::const_iterator iter = fConfigVect.begin();
    
    while (iter != fConfigVect.end()) {
      if (iter->second == G4AttributeFilterT<T>::Interval) {filter->LoadIntervalElement(iter->first);}
      else if (iter->second == G4AttributeFilterT<T>::SingleValue) {filter->LoadSingleValueElement(iter->first);}
      iter++;
    }
  }

  // Get attribute value
  G4AttValue attVal;

  // Expect value to exist
  if (!G4AttUtils::ExtractAttValue(object, fAttName, attVal)) {
    static G4bool warnedUnableToExtract = false;
    if (!warnedUnableToExtract) {
      G4ExceptionDescription ed;
      ed <<"Unable to extract attribute definition named "<<fAttName<<'\n'
      << "Available attributes:\n"
      << *object.GetAttDefs();
      G4Exception
	("G4AttributeFilterT::Evaluate", "modeling0103", JustWarning, ed, "InvalidAttributeValue");
      warnedUnableToExtract = true;
    }
    return false;
  }

  if (G4SmartFilter<T>::GetVerbose()) {
    G4cout<<"G4AttributeFilterT processing attribute named "<<fAttName;
    G4cout<<" with value "<<attVal.GetValue()<<G4endl;
  }

  // Pass subfilter
  return (filter->Accept(attVal));
}

template <typename T>
void
G4AttributeFilterT<T>::Clear()
{
  fConfigVect.clear();
  if (0 != filter) filter->Reset();
}

template <typename T>
void
G4AttributeFilterT<T>::Print(std::ostream& ostr) const
{
  ostr<<"Printing data for G4Attribute filter named: "<<G4VFilter<T>::Name()<<std::endl;
  ostr<<"Filtered attribute name: "<<fAttName<<std::endl;
  ostr<<"Printing sub filter data:"<<std::endl;
  if (0 != filter) filter->PrintAll(ostr);
}

template <typename T>
void
G4AttributeFilterT<T>::Set(const G4String& name)
{
  fAttName = name;
}

template <typename T>
void
G4AttributeFilterT<T>::AddInterval(const G4String& interval)
{
  std::pair<G4String, Config> myPair(interval, G4AttributeFilterT<T>::Interval);

  typename ConfigVect::iterator iter = std::find(fConfigVect.begin(), fConfigVect.end(), myPair);
  
  if (iter != fConfigVect.end()) {
    G4ExceptionDescription ed;
    ed <<"Interval "<< interval <<" already exists";
    G4Exception
      ("G4AttributeFilterT::AddInterval", "modeling0104", JustWarning, ed);
    return;
  }

  fConfigVect.push_back(myPair);
}

template <typename T>
void
G4AttributeFilterT<T>::AddValue(const G4String& value)
{
  std::pair<G4String, Config> myPair(value, G4AttributeFilterT<T>::SingleValue);

  typename ConfigVect::iterator iter = std::find(fConfigVect.begin(), fConfigVect.end(), myPair);
  
  if (iter != fConfigVect.end()) {
    G4ExceptionDescription ed;
    ed <<"Single value "<< value <<" already exists";
    G4Exception
      ("G4AttributeFilterT::AddValue", "modeling0105", JustWarning, ed);
    return;
  }
  fConfigVect.push_back(myPair);
}

#endif
