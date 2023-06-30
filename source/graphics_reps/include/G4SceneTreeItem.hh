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

#ifndef G4SceneTreeItem_hh
#define G4SceneTreeItem_hh

#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

#include <list>
#include <map>
#include <vector>

class G4SceneTreeItem
{
  public:
    // A ghost is a touchable we know to be there but know only its path
    enum Type
    {
      unidentified,
      root,
      model,
      pvmodel,  // Physical Volume model - special case 
      ghost,
      touchable
    };

    explicit G4SceneTreeItem(Type type = unidentified) { fType = type; }
    ~G4SceneTreeItem() = default;

    // Copy contructor copies the whole tree, i.e., children and all descendants
    G4SceneTreeItem(const G4SceneTreeItem&) = default;

    // Assigns the whole tree, i.e., children and all descendants
    G4SceneTreeItem& operator= (const G4SceneTreeItem&) = default;

    // Access functions

    Type GetType() const { return fType; }
    const G4String& GetTypeString() const { return fTypeMap[fType]; }
    void SetType(Type type) { fType = type; }

    const G4String& GetPVPath() const { return fPVPath; }
    void SetPVPath(const G4String& PVPath) { fPVPath = PVPath; }

    const G4String& GetDescription() const { return fDescription; }
    void SetDescription(const G4String& description) { fDescription = description; }

    const G4String& GetModelType() const { return fModelType; }
    void SetModelType(const G4String& modelType) { fModelType = modelType; }

    const G4String& GetModelDescription() const { return fModelDescription; }
    void SetModelDescription(const G4String& modelDescription)
    {
      fModelDescription = modelDescription;
    }

    const std::map<G4String, G4AttDef>* GetAttDefs() const { return fpAttDefs; }
    void SetAttDefs(const std::map<G4String, G4AttDef>* pAttDefs) { fpAttDefs = pAttDefs; };

    std::vector<G4AttValue>* GetAttValues() const { return fpAttValues; }
    void SetAttValues(std::vector<G4AttValue>* pAttValues) { fpAttValues = pAttValues; }

    const G4VisAttributes& GetVisAttributes() const { return fVisAttributes; }
    G4VisAttributes& AccessVisAttributes() { return fVisAttributes; }
    void SetVisAttributes(const G4VisAttributes& visAtts) { fVisAttributes = visAtts; }

    const std::list<G4SceneTreeItem>& GetChildren() const { return fChildren; }
    // Insert item at - or rather, just before - pos
    std::list<G4SceneTreeItem>::iterator InsertChild(std::list<G4SceneTreeItem>::iterator pos,
                                                     const G4SceneTreeItem& item)
    {
      return fChildren.insert(pos, item);
    }
    std::list<G4SceneTreeItem>& AccessChildren() { return fChildren; }

    G4bool IsExpanded() const         { return fExpanded; }
    void SetExpanded(G4bool expanded) { fExpanded = expanded; }

    // Utility functions

    // Reset visibility of all objects to false - visible objects will then set to true
    void ResetVisibility();

    // If found, returns "true" and places iterator in foundIter
    G4bool FindTouchableFromRoot(const G4String& fullPathString,
                                 std::list<G4SceneTreeItem>::iterator& foundIter);

    // Dump single item, i.e., ignore any children
    void DumpSingleItem(std::ostream&, G4int verbosity = 0) const;

    // Dump whole tree
    void DumpTree(std::ostream&, G4int verbosity = 0) const;

  private:
    Type fType = unidentified;
    static std::map<Type, G4String> fTypeMap;
    G4String fDescription;
    G4String fModelType = "none";
    G4String fModelDescription;
    G4String fPVPath;
    G4VisAttributes fVisAttributes;
    const std::map<G4String, G4AttDef>* fpAttDefs = nullptr;
    std::vector<G4AttValue>* fpAttValues = nullptr;
    std::list<G4SceneTreeItem> fChildren;
    G4bool fExpanded = true;
};

#endif  // G4SceneTreeItem_hh
