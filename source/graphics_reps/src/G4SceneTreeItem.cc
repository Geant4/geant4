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

#include "G4SceneTreeItem.hh"

#include "G4AttCheck.hh"

#include <iostream>

// A ghost is a touchable we know to be there but know only its path
std::map<G4SceneTreeItem::Type, G4String> G4SceneTreeItem::fTypeMap = {
  {G4SceneTreeItem::unidentified, "unidentified"},
  {G4SceneTreeItem::root, "root"},
  {G4SceneTreeItem::ghost, "ghost"},
  {G4SceneTreeItem::model, "model"},
  {G4SceneTreeItem::pvmodel, "pvmodel"},
  {G4SceneTreeItem::touchable, "touchable"}};

// Reset visibility of all objects to false - visible objects will then set to true
void G4SceneTreeItem::ResetVisibility()
{
  // Reset all but the root item, which is always visible (i.e., active)
  // Other items will be set visible if and when presented to the scene
  if (fType != root) fVisAttributes.SetVisibility(false);
  for (auto& child : fChildren)
    child.ResetVisibility();
}

// If found, returns "true" and places iterator in foundIter
G4bool G4SceneTreeItem::FindTouchableFromRoot(const G4String& fullPathString,
                                              std::list<G4SceneTreeItem>::iterator& foundIter)
{
  if (fType != root) {
    G4ExceptionDescription ed;
    ed << "Not a root item:\n";
    DumpSingleItem(ed);
    G4Exception("G4SceneTreeItem::FindTouchableFromRoot", "greps0011", JustWarning, ed);
    return false;
  }

  for (auto& aModel : fChildren) {
    if (aModel.fModelType == "G4PhysicalVolumeModel") {  // Top item, i.e., root of touchables
      // Work down the path - "name id", then "name id name id", etc.
      G4String partialPathString;
      auto iter = aModel.fChildren.begin();
      auto iterEnd = aModel.fChildren.end();
      std::istringstream iss(fullPathString);
      G4String name, copyNo;
      while (iss >> name >> copyNo) {
        partialPathString += ' ' + name + ' ' + copyNo;
        for (; iter != iterEnd; ++iter) {
          if (iter->fPVPath == partialPathString) {
            if (partialPathString != fullPathString) {
              // Go to next level
              iter = iter->fChildren.begin();
              iterEnd = iter->fChildren.end();
            }
            break;
          }
        }
        if (iter != iterEnd) {  // Found
          foundIter = iter;
          return true;
        }
      }
    }
  }
  return false;
}

// Dump single item, i.e., ignore any children
void G4SceneTreeItem::DumpSingleItem(std::ostream& os, G4int verbosity) const
{
  static G4bool first = true;
  if (first) {
    first = false;
    os << "  Verbosity actions:"
       << "\n  >=0 one line"
       << "\n  >=1 a few lines"
       << "\n  >=2 check G4Atts"
       << "\n  >=3 print G4Atts"
       << "\n  >=4 print some attValues"
       << '\n';
  }

  os << GetTypeString() << " (";
  G4String status;
  switch (fType) {
    case unidentified:
      status = "error";
      break;
    case root:
      status = "active";
      break;
    case model:
      [[fallthrough]];
    case pvmodel:
      status = (fVisAttributes.IsVisible() ? "active" : "inactive");
      break;
    case ghost:
      [[fallthrough]];
    case touchable:
      status = (fVisAttributes.IsVisible() ? "visible" : "invisible");
      break;
  }
  if (fExpanded) {
    status += ",expanded";
  } else {
    status += ",collapsed";
  }
  os << status << ')';

  G4String description;
  switch (fType) {
    case unidentified:
      break;
    case root:
      break;
    case model:
      [[fallthrough]];
    case pvmodel:
      description = " \"" + fModelDescription + '"';
      break;
    case ghost:
      [[fallthrough]];
    case touchable:
      description = " (" + fPVPath.substr(1, fPVPath.length() - 1) + ')';
      break;
  }
  os << description;

  if (fType == touchable) {
    os << ' ' << fVisAttributes.GetColour();
  }

  // clang-format off
  if (verbosity >= 1) {
    os << "\n  Description: " << GetDescription()
    << "\n  Model type: " << GetModelType()
    << "\n  Model description: " << GetModelDescription()
    << "\n  " << fChildren.size() << (fChildren.size()==1?" child":" children");
  }
  // clang-format on

  if (verbosity >= 2) {
    const auto& attDefs = GetAttDefs();
    const auto& attValues = GetAttValues();
    if (attDefs == nullptr || attValues == nullptr) {
      os << "\n  No G4Atts";  // Legitimate
    }
    else {
      G4AttCheck attCheck(attValues, attDefs);
      if (attCheck.Check("G4SceneTreeItem::Dump")) {
        G4ExceptionDescription ed;
        ed << "Item: " << attCheck;
        G4Exception("G4SceneTreeItem::Dump", "greps0010", JustWarning, ed,
                    "G4Atts don't check out");
        return;
      }

      if (verbosity >= 3) {
        os << "\n  G4Atts:\n" << attCheck;
        static G4bool first1 = true;
        if (first1) {
          first1 = false;
          os << "\n  Available G4Atts for touchable:";
          // clang-format off
          for (const auto& att: *GetAttDefs()) {
            os << "\n  " << att.first
            << ',' << att.second.GetName()
            << ',' << att.second.GetDesc()
            << ',' << att.second.GetCategory()
            << ',' << att.second.GetExtra()
            << ',' << att.second.GetValueType()
            << ',' << att.second.GetTypeKey();
          }
          // clang-format on
          os << '\n';
        }
      }

      if (verbosity >= 4) {
        for (G4String name : {"PVPath", "GlobalExtent"}) {
          G4String result;
          const auto& iterAttDef = attDefs->find(name);
          if (iterAttDef != attDefs->end()) {
            result = result + "\n  " + iterAttDef->second.GetName();
            result = result + ", " + iterAttDef->second.GetDesc();
          }
          const G4AttValue* pAttValue = nullptr;
          for (const auto& attValue : *attValues) {
            // Why are the attValues not in a map like the attDefs???
            if (attValue.GetName() == name) {
              pAttValue = &attValue;  // Avoid copy
              break;
            }
          }
          if (pAttValue) {
            result = result + ", " + pAttValue->GetValue();
          }

          os << result;
        }
      }
    }
  }

  os << std::endl;
}

// Dump whole tree
void G4SceneTreeItem::DumpTree(std::ostream& os, G4int verbosity) const
{
  static G4int depth = 0;
  for (G4int i = 0; i < depth; i++) os << "  ";
  DumpSingleItem(os, verbosity);
  for (auto& child : GetChildren()) {
    depth++;
    child.DumpTree(os, verbosity);
    depth--;
  }
}
