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
//
// $Id: G4AttCheck.cc,v 1.2 2005-03-26 22:37:14 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4AttCheck.hh"

#include "globals.hh"

#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"
#include "G4UIcommand.hh"

#include <algorithm>
#include <cctype>

G4AttCheck::G4AttCheck
(const std::vector<G4AttValue>* values,
 const std::map<G4String,G4AttDef>* definitions):
  fpValues(values),
  fpDefinitions(definitions)
{
  // Legal Unit Category Types...
  fUnitCategories.insert("Length");
  fUnitCategories.insert("Energy");

  // Corresponding Standard Units...
  fStandardUnits["Length"] = "m";
  fStandardUnits["Energy"] = "MeV";

  // Legal Value Types...
  fValueTypes.insert("G4String");
  fValueTypes.insert("G4BestUnit");
  // ...plus any legal unit symbol ("MeV", "km", etc.)...
  G4UnitsTable& units = G4UnitDefinition::GetUnitsTable();
  for (size_t i = 0; i < units.size(); ++i) {
    if (fUnitCategories.find(units[i]->GetName()) !=
	fUnitCategories.end()) {
      //G4cout << units[i]->GetName() << G4endl;
      G4UnitsContainer& container = units[i]->GetUnitsList();
      for (size_t j = 0; j < container.size(); ++j) {
	//G4cout << container[j]->GetName() << ' '
	//       << container[j]->GetSymbol() << G4endl;
	fValueTypes.insert(container[j]->GetSymbol());
      }
   }
  }
}

G4AttCheck::~G4AttCheck() {}

std::set<G4String> G4AttCheck::fUnitCategories;

std::map<G4String,G4String> G4AttCheck::fStandardUnits;

std::set<G4String> G4AttCheck::fValueTypes;

void G4AttCheck::Check() const {
  // Check only.  Silent unless error - then G4cerr.
  using namespace std;
  vector<G4AttValue>::const_iterator iValue;
  for (iValue = fpValues->begin(); iValue != fpValues->end(); ++iValue) {
    map<G4String,G4AttDef>::const_iterator iDef =
      fpDefinitions->find(iValue->GetName());
    if (iDef == fpDefinitions->end()) {
      G4cerr <<
	"\n*******************************************************"
	"\nERROR: No G4AttDef for G4AttValue \""
	     <<  iValue->GetName()
	     << "\": "
	     << iValue->GetValue() <<
	"\n*******************************************************"
	     << G4endl;
    } else {
      if (fValueTypes.find(iDef->second.GetValueType()) == fValueTypes.end()) {
	G4cerr <<
	  "\n*******************************************************"
	  "\nERROR: Illegal value type \""
	       << iDef->second.GetValueType()
	       << "\" for G4AttValue \""
	       <<  iValue->GetName()
	       << "\": "
	       << iValue->GetValue() <<
	  "\n  Possible value types:";
	std::set<G4String>::iterator i;
	for (i = fValueTypes.begin(); i != fValueTypes.end(); ++i) {
	  G4cerr << ' ' << *i;
	}
	G4cerr <<
	  "\n*******************************************************"
	       << G4endl;
      }
    }
  }
}

void G4AttCheck::AddValuesAndDefs
(std::vector<G4AttValue>* pValues,
 std::map<G4String,G4AttDef>* pDefinitions,
 const G4String& oldName,
 const G4String& name,
 const G4String& value,
 const G4String& valueType,
 const G4String& description) const {
  // Add new G4AttDeff...
  pValues->push_back(G4AttValue(name,value,""));
  // Copy original G4AttDef...
  (*pDefinitions)[name] = fpDefinitions->find(oldName)->second;
  // ...and make appropriate changes...
  (*pDefinitions)[name].SetName(name);
  (*pDefinitions)[name].SetValueType(valueType);
  if (description != "") (*pDefinitions)[name].SetDesc(description);
}

G4AttCheck G4AttCheck::Standard() const {
  // Returns standard versions on the heap.
  using namespace std;

  vector<G4AttValue>* pValues = new vector<G4AttValue>;
  map<G4String,G4AttDef>* pDefinitions = new map<G4String,G4AttDef>;

  vector<G4AttValue>::const_iterator iValue;
  for (iValue = fpValues->begin(); iValue != fpValues->end(); ++iValue) {
    const G4String& valueName = iValue->GetName();
    //G4cout << "G4AttValue \"" << valueName << "\" found..." << G4endl;
    map<G4String,G4AttDef>::const_iterator iDef =
      fpDefinitions->find(valueName);
    if (iDef != fpDefinitions->end()) {
      const G4String& valueType = iDef->second.GetValueType();
      const G4String& value = iValue->GetValue();
      if (fValueTypes.find(valueType) != fValueTypes.end()) {
	if (valueType == "G4String") {  // Simply copy...
	  pValues->push_back(*iValue);
	  (*pDefinitions)[valueName] =
	    fpDefinitions->find(valueName)->second;
	}
	else {
	  G4String valueAndUnit;
	  G4String unit;
	  if (valueType == "G4BestUnit") {
	    valueAndUnit = value;
	    valueAndUnit = valueAndUnit.strip();
	    unit = valueAndUnit.substr(valueAndUnit.rfind(' ')+1);
	  } else {
	    valueAndUnit = value + ' ' + valueType;
	    unit = valueType;
	  }
	  G4String category = G4UnitDefinition::GetCategory(unit);
	  if (fUnitCategories.find(category) != fUnitCategories.end()) {
	    G4String standardUnit = fStandardUnits[category];
	    G4String standardValueType = standardUnit;
	    G4double valueOfUnit = G4UnitDefinition::GetValueOf(standardUnit);
	    G4String newValue;
	    G4String extra = iDef->second.GetExtra();
	    /*
	    G4cout << "valueName = \"" << valueName
		   << "\", valueAndUnit = \"" << valueAndUnit
		   << "\", unit = \"" << unit
		   << "\", extra = \"" << extra
		   << "\", category \"" << G4UnitDefinition::GetCategory(unit)
		   << G4endl;
	    */
	    transform(extra.begin(),extra.end(),extra.begin(),::tolower);
	    if (extra.find("vector") == string::npos) {
	      G4double internalValue =
		G4UIcommand::ConvertToDimensionedDouble(valueAndUnit);
	      AddValuesAndDefs
		(pValues,pDefinitions,
		 valueName,valueName,
		 G4UIcommand::ConvertToString(internalValue/valueOfUnit),
		 standardValueType);
	    } else {
	      G4ThreeVector internalValue =
		G4UIcommand::ConvertToDimensioned3Vector(valueAndUnit);
	      AddValuesAndDefs
		(pValues,pDefinitions,
		 valueName,valueName+"-X",
		 G4UIcommand::ConvertToString(internalValue.x()/valueOfUnit),
		 standardValueType,
		 fpDefinitions->find(valueName)->second.GetDesc()+"-X");
	      AddValuesAndDefs
		(pValues,pDefinitions,
		 valueName,valueName+"-Y",
		 G4UIcommand::ConvertToString(internalValue.y()/valueOfUnit),
		 standardValueType,
		 fpDefinitions->find(valueName)->second.GetDesc()+"-Y");
	      AddValuesAndDefs
		(pValues,pDefinitions,
		 valueName,valueName+"-Z",
		 G4UIcommand::ConvertToString(internalValue.z()/valueOfUnit),
		 standardValueType,
		 fpDefinitions->find(valueName)->second.GetDesc()+"-Z");
	    }
	  }
	}
      }
    }
  }
  return G4AttCheck(pValues,pDefinitions);
}

std::ostream& operator<< (std::ostream& os, const G4AttCheck& ac) {
  using namespace std;
  vector<G4AttValue>::const_iterator iValue;
  for (iValue = ac.fpValues->begin(); iValue != ac.fpValues->end(); ++iValue) {
    map<G4String,G4AttDef>::const_iterator iDef =
      ac.fpDefinitions->find(iValue->GetName());
    if (iDef == ac.fpDefinitions->end()) {
      os << "ERROR: No G4AttDef for G4AttValue \""
	 <<  iValue->GetName() << "\": " << iValue->GetValue() << endl;
    } else {
      if (ac.fValueTypes.find(iDef->second.GetValueType()) ==
	  ac.fValueTypes.end()) {
	os << "ERROR: Illegal value type \""
	   << iDef->second.GetValueType()
	   << "\" for G4AttValue \""
	   <<  iValue->GetName()
	   << "\": "
	   << iValue->GetValue() <<
	  "\n  Possible value types:";
	std::set<G4String>::iterator i;
	for (i = ac.fValueTypes.begin(); i != ac.fValueTypes.end(); ++i) {
	  os << ' ' << *i;
	}
	os << endl;
     } else {
	os << iDef->second.GetDesc() << ": "
	   << iValue->GetValue() << " ("
	   << iDef->second.GetValueType() << ")"
	   << endl;
      }
    }
  }
  return os;
}
