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

#include "G4AttCheck.hh"

#include "globals.hh"

#include "G4AttDef.hh"
#include "G4AttDefStore.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"
#include "G4UIcommand.hh"

G4ThreadLocal G4bool G4AttCheck::fFirst = true;

G4ThreadLocal std::set<G4String> *G4AttCheck::fUnitCategories = nullptr;

G4ThreadLocal std::map<G4String,G4String> *G4AttCheck::fStandardUnits = nullptr;

G4ThreadLocal std::set<G4String> *G4AttCheck::fCategories = nullptr;

G4ThreadLocal std::set<G4String> *G4AttCheck::fUnits = nullptr;

G4ThreadLocal std::set<G4String> *G4AttCheck::fValueTypes = nullptr;

G4AttCheck::G4AttCheck
(const std::vector<G4AttValue>* values,
 const std::map<G4String,G4AttDef>* definitions):
  fpValues(values),
  fpDefinitions(definitions)
{
  Init();

  if (fFirst) {  // Initialise static containers.
    fFirst = false;

    // Legal Unit Category Types...
    fUnitCategories->insert("Length");
    fUnitCategories->insert("Energy");
    fUnitCategories->insert("Time");
    fUnitCategories->insert("Electric charge");
    fUnitCategories->insert("Volumic Mass");  // (Density)

    // Corresponding Standard Units...
    (*fStandardUnits)["Length"] = "m";
    (*fStandardUnits)["Energy"] = "MeV";
    (*fStandardUnits)["Time"] = "ns";
    (*fStandardUnits)["Electric charge"] = "e+";
    (*fStandardUnits)["Volumic Mass"] = "kg/m3";

    // Legal Categories...
    fCategories->insert("Bookkeeping");
    fCategories->insert("Draw");
    fCategories->insert("Physics");
    fCategories->insert("PickAction");
    fCategories->insert("Association");

    // Legal units...
    fUnits->insert("");
    fUnits->insert("G4BestUnit");
    // ...plus any legal unit symbol ("MeV", "km", etc.)...
    G4UnitsTable& units = G4UnitDefinition::GetUnitsTable();
    for (size_t i = 0; i < units.size(); ++i) {
      if (fUnitCategories->find(units[i]->GetName()) !=
          fUnitCategories->end()) {
        //G4cout << units[i]->GetName() << G4endl;
        G4UnitsContainer& container = units[i]->GetUnitsList();
        for (auto & j : container) {
          //G4cout << container[j]->GetName() << ' '
          //       << container[j]->GetSymbol() << G4endl;
          fUnits->insert(j->GetSymbol());
        }
      }
    }

    // Legal Value Types...
    fValueTypes->insert("G4String");
    fValueTypes->insert("G4int");
    fValueTypes->insert("G4double");
    fValueTypes->insert("G4ThreeVector");
    fValueTypes->insert("G4bool");
  }
}

void G4AttCheck::Init()
{
  if (fValueTypes == nullptr) fValueTypes = new std::set<G4String>;
  if (fUnits == nullptr) fUnits = new std::set<G4String>;
  if (fCategories == nullptr) fCategories = new std::set<G4String>;
  if (fStandardUnits == nullptr) fStandardUnits = new std::map<G4String,G4String>;
  if (fUnitCategories == nullptr) fUnitCategories = new std::set<G4String>;
}

G4bool G4AttCheck::Check(const G4String& leader) const
{
  // Check only.  Silent unless error - then G4cerr.  Returns error.
  G4bool error = false;
  static G4ThreadLocal G4int iError = 0;
  G4bool print = false;
  if (iError < 10 || iError%100 == 0) {
    print = true;
  }
  using namespace std;
  if (fpValues == nullptr) return error;  // A null values vector is a valid situation.
  if (fpDefinitions == nullptr) {
    ++iError;
    error = true;
    if (print) {
      G4cerr <<
        "\n*******************************************************";
      if (!leader.empty()) {
        G4cerr << '\n' << leader;
      }
      G4cerr <<
        "\nG4AttCheck: ERROR " << iError << ": Null definitions pointer"
        "\n*******************************************************"
             << G4endl;
    }
    return error;
  }
  vector<G4AttValue>::const_iterator iValue;
  for (iValue = fpValues->begin(); iValue != fpValues->end(); ++iValue) {
    const G4String& valueName = iValue->GetName();
    const G4String& value = iValue->GetValue();
    // NOLINTNEXTLINE(modernize-use-auto): Explicitly want a const_iterator
    map<G4String,G4AttDef>::const_iterator iDef =
      fpDefinitions->find(valueName);
    if (iDef == fpDefinitions->end()) {
      ++iError;
      error = true;
      if (print) {
        G4cerr <<
          "\n*******************************************************";
        if (!leader.empty()) {
          G4cerr << '\n' << leader;
        }
        G4cerr <<
          "\nG4AttCheck: ERROR " << iError << ": No G4AttDef for G4AttValue \""
               <<  valueName << "\": " << value <<
          "\n*******************************************************"
               << G4endl;
      }
    } else {
      const G4String& category = iDef->second.GetCategory();
      const G4String& extra = iDef->second.GetExtra();
      const G4String& valueType = iDef->second.GetValueType();
      if (fCategories->find(category) == fCategories->end()) {
        ++iError;
        error = true;
        if (print) {
          G4cerr <<
            "\n*******************************************************";
          if (!leader.empty()) {
            G4cerr << '\n' << leader;
          }
          G4cerr <<
            "\nG4AttCheck: ERROR " << iError << ": Illegal Category Field \""
                 << category << "\" for G4AttValue \""
                 << valueName << "\": " << value <<
            "\n  Possible Categories:";
          set<G4String>::iterator i;
          for (i = fCategories->begin(); i != fCategories->end(); ++i) {
            G4cerr << ' ' << *i;
          }
          G4cerr <<
            "\n*******************************************************"
                 << G4endl;
        }
      }
      if(category == "Physics" && fUnits->find(extra) == fUnits->end()) {
        ++iError;
        error = true;
        if (print) {
          G4cerr <<
            "\n*******************************************************";
          if (!leader.empty()) {
            G4cerr << '\n' << leader;
          }
          G4cerr <<
            "\nG4AttCheck: ERROR " << iError << ": Illegal Extra field \""
                 << extra << "\" for G4AttValue \""
                 << valueName << "\": " << value <<
            "\n  Possible Extra fields if Category==\"Physics\":\n    ";
          set<G4String>::iterator i;
          for (i = fUnits->begin(); i != fUnits->end(); ++i) {
            G4cerr << ' ' << *i;
          }
          G4cerr <<
            "\n*******************************************************"
                 << G4endl;
        }
      }
      if (fValueTypes->find(valueType) == fValueTypes->end()) {
        ++iError;
        error = true;
        if (print) {
          G4cerr <<
            "\n*******************************************************";
          if (!leader.empty()) {
            G4cerr << '\n' << leader;
          }
          G4cerr <<
            "\nG4AttCheck: ERROR " << iError << ": Illegal Value Type field \""
                 << valueType << "\" for G4AttValue \""
                 << valueName << "\": " << value <<
            "\n  Possible Value Types:";
          set<G4String>::iterator i;
          for (i = fValueTypes->begin(); i != fValueTypes->end(); ++i) {
            G4cerr << ' ' << *i;
          }
          G4cerr <<
            "\n*******************************************************"
               << G4endl;
        }
      }
    }
  }
  return error;
}

std::ostream& operator<< (std::ostream& os, const G4AttCheck& ac)
{
  using namespace std;
  if (ac.fpDefinitions == nullptr) {
    os << "G4AttCheck: ERROR: zero definitions pointer." << endl;
    return os;
  }
  G4String storeKey;
  if (G4AttDefStore::GetStoreKey(ac.fpDefinitions, storeKey)) {
    os << storeKey << ':' << endl;
  }
  if (ac.fpValues == nullptr) {
    // A null values vector is a valid situation.
    os << "G4AttCheck: zero values pointer." << endl;
    return os;
  }
  vector<G4AttValue>::const_iterator iValue;
  for (iValue = ac.fpValues->begin(); iValue != ac.fpValues->end(); ++iValue) {
    const G4String& valueName = iValue->GetName();
    const G4String& value = iValue->GetValue();
    // NOLINTNEXTLINE(modernize-use-auto): Explicitly want a const_iterator
    map<G4String,G4AttDef>::const_iterator iDef = 
      ac.fpDefinitions->find(valueName);
    G4bool error = false;
    if (iDef == ac.fpDefinitions->end()) {
      error = true;
      os << "G4AttCheck: ERROR: No G4AttDef for G4AttValue \""
         << valueName << "\": " << value << endl;
    } else {
      const G4String& category = iDef->second.GetCategory();
      const G4String& extra = iDef->second.GetExtra();
      const G4String& valueType = iDef->second.GetValueType();
      if (ac.fCategories->find(category) == ac.fCategories->end()) {
        error = true;
        os <<
          "G4AttCheck: ERROR: Illegal Category Field \"" << category
           << "\" for G4AttValue \"" << valueName << "\": " << value <<
          "\n  Possible Categories:";
        set<G4String>::iterator i;
        for (i = ac.fCategories->begin(); i != ac.fCategories->end(); ++i) {
          os << ' ' << *i;
        }
        os << endl;
      }
      if(category == "Physics" && ac.fUnits->find(extra) == ac.fUnits->end()) {
        error = true;
        os <<
          "G4AttCheck: ERROR: Illegal Extra field \""<< extra
           << "\" for G4AttValue \"" << valueName << "\": " << value <<
          "\n  Possible Extra fields if Category==\"Physics\":\n    ";
        set<G4String>::iterator i;
        for (i = ac.fUnits->begin(); i != ac.fUnits->end(); ++i) {
          os << ' ' << *i;
        }
        os << endl;
      }
      if (ac.fValueTypes->find(valueType) == ac.fValueTypes->end()) {
        error = true;
        os <<
          "G4AttCheck: ERROR: Illegal Value Type field \"" << valueType
           << "\" for G4AttValue \"" << valueName << "\": " << value <<
          "\n  Possible Value Types:";
        set<G4String>::iterator i;
        for (i = ac.fValueTypes->begin(); i != ac.fValueTypes->end(); ++i) {
          os << ' ' << *i;
        }
        os << endl;
      }
    }
    if (!error) {
      os << iDef->second.GetDesc()
         << " (" << valueName
         << "): " << value;
      if (iDef->second.GetCategory() == "Physics" &&
          !iDef->second.GetExtra().empty()) {
        os << " (" << iDef->second.GetExtra() << ")";
      }
      os << endl;
    }
  }
  return os;
}

void G4AttCheck::AddValuesAndDefs
(std::vector<G4AttValue>* standardValues,
 std::map<G4String,G4AttDef>* standardDefinitions,
 const G4String& oldName,
 const G4String& name,
 const G4String& value,
 const G4String& extra,
 const G4String& description) const
{
  // Add new G4AttDeff...
  standardValues->push_back(G4AttValue(name,value,""));
  // Copy original G4AttDef...
  (*standardDefinitions)[name] = fpDefinitions->find(oldName)->second;
  // ...and make appropriate changes...
  (*standardDefinitions)[name].SetName(name);
  (*standardDefinitions)[name].SetExtra(extra);
  if (!description.empty()) (*standardDefinitions)[name].SetDesc(description);
}

G4bool G4AttCheck::Standard
(std::vector<G4AttValue>* standardValues,
 std::map<G4String,G4AttDef>* standardDefinitions) const
{
  // Places standard versions in provided vector and map and returns error.
  // Assumes valid input.  Use Check to check.
  using namespace std;
  G4bool error = false;
  vector<G4AttValue>::const_iterator iValue;
  for (iValue = fpValues->begin(); iValue != fpValues->end(); ++iValue) {
    const G4String& valueName = iValue->GetName();
    const G4String& value = iValue->GetValue();
    // NOLINTNEXTLINE(modernize-use-auto): Explicitly want a const_iterator
    map<G4String,G4AttDef>::const_iterator iDef =
      fpDefinitions->find(valueName);
    if (iDef == fpDefinitions->end()) {
      error = true;
    } else {
      const G4String& category = iDef->second.GetCategory();
      const G4String& extra = iDef->second.GetExtra();
      const G4String& valueType = iDef->second.GetValueType();
      if (fCategories->find(category) == fCategories->end() ||
          (category == "Physics" && fUnits->find(extra) == fUnits->end()) ||
          fValueTypes->find(valueType) == fValueTypes->end()) {
        error = true;
      } else {
        if (category != "Physics") {  // Simply copy...
          standardValues->push_back(*iValue);
          (*standardDefinitions)[valueName] =
            fpDefinitions->find(valueName)->second;
        } else {  // "Physics"...
          if (extra.empty()) {  // Dimensionless...
            if (valueType == "G4ThreeVector") {  // Split vector into 3...
              G4ThreeVector internalValue =
                G4UIcommand::ConvertTo3Vector(value);
              AddValuesAndDefs
                (standardValues,standardDefinitions,
                 valueName,valueName+"-X",
                 G4UIcommand::ConvertToString(internalValue.x()),"",
                 fpDefinitions->find(valueName)->second.GetDesc()+"-X");
              AddValuesAndDefs
                (standardValues,standardDefinitions,
                 valueName,valueName+"-Y",
                 G4UIcommand::ConvertToString(internalValue.y()),"",
                 fpDefinitions->find(valueName)->second.GetDesc()+"-Y");
              AddValuesAndDefs
                (standardValues,standardDefinitions,
                 valueName,valueName+"-Z",
                 G4UIcommand::ConvertToString(internalValue.z()),"",
                 fpDefinitions->find(valueName)->second.GetDesc()+"-Z");
            } else {  // Simply copy...
              standardValues->push_back(*iValue);
              (*standardDefinitions)[valueName] =
                fpDefinitions->find(valueName)->second;
            }
          } else {  // Dimensioned...
            G4String valueAndUnit;
            G4String unit;
            if (extra == "G4BestUnit") {
              valueAndUnit = G4StrUtil::strip_copy(value);
              unit = valueAndUnit.substr(valueAndUnit.rfind(' ')+1);
            } else {
              valueAndUnit = G4StrUtil::strip_copy(value + ' ' + extra);
              unit = extra;
            }
            G4String unitCategory = G4UnitDefinition::GetCategory(unit);
            if (fUnitCategories->find(unitCategory) != fUnitCategories->end()) {
              G4String standardUnit = (*fStandardUnits)[unitCategory];
              G4double valueOfStandardUnit =
                G4UnitDefinition::GetValueOf(standardUnit);
//              G4String exstr = iDef->second.GetExtra();
              if (valueType == "G4ThreeVector") {  // Split vector into 3...
                G4ThreeVector internalValue =
                  G4UIcommand::ConvertToDimensioned3Vector(valueAndUnit);
                AddValuesAndDefs
                  (standardValues,standardDefinitions,
                   valueName,valueName+"-X",
                   G4UIcommand::ConvertToString
                   (internalValue.x()/valueOfStandardUnit),
                   standardUnit,
                   fpDefinitions->find(valueName)->second.GetDesc()+"-X");
                AddValuesAndDefs
                  (standardValues,standardDefinitions,
                   valueName,valueName+"-Y",
                   G4UIcommand::ConvertToString
                   (internalValue.y()/valueOfStandardUnit),
                   standardUnit,
                   fpDefinitions->find(valueName)->second.GetDesc()+"-Y");
                AddValuesAndDefs
                  (standardValues,standardDefinitions,
                   valueName,valueName+"-Z",
                   G4UIcommand::ConvertToString
                   (internalValue.z()/valueOfStandardUnit),
                   standardUnit,
                   fpDefinitions->find(valueName)->second.GetDesc()+"-Z");
              } else {
                G4double internalValue =
                  G4UIcommand::ConvertToDimensionedDouble(valueAndUnit);
                AddValuesAndDefs
                  (standardValues,standardDefinitions,
                   valueName,valueName,
                   G4UIcommand::ConvertToString
                   (internalValue/valueOfStandardUnit),
                   standardUnit);
              }
            }
          }
        }
      }
    }
  }
  if (error) {
    G4cerr << "G4AttCheck::Standard: Conversion error." << G4endl;
  }
  return error;
}
