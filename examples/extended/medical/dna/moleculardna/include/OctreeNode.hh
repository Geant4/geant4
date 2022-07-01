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
#ifndef MOLECULAR_OCTREE_NODE_HH
#define MOLECULAR_OCTREE_NODE_HH

#include "globals.hh"
#include "G4ThreeVector.hh"

#include <vector>
#include <array>

class G4VPhysicalVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class OctreeNode
{
 public:
  // Uniform divisions constructor along each axis
  OctreeNode(G4ThreeVector, G4ThreeVector, G4int, OctreeNode* parent = nullptr);

  virtual ~OctreeNode();

  inline G4bool HasChildren() const
  {
    return (fChildren[0] != nullptr);
  }

  inline OctreeNode* GetParent() const { return fParent; };

  inline const auto& GetHalfLengths() const
  {
    return fHalfLengths;
  };

  inline G4double GetHalfLengthsMag() const
  {
    return fHalfLengthsMag;
  };

  inline const G4ThreeVector& GetPosition() const
  {
    return fPosition;
  };

  inline const auto& GetChildren() const { return fChildren; };

  const std::vector<G4VPhysicalVolume*> SearchOctree(
    const G4ThreeVector&, G4double _rad = 0) const;

  void SearchOctree(const G4ThreeVector& pos,
                    std::vector<G4VPhysicalVolume*>& out,
                    G4double _rad = 0) const;

  const std::vector<G4VPhysicalVolume*> SearchOctree(
    const G4ThreeVector&) const;

  G4int GetNumberOfTerminalNodes();

  void AddPhysicalVolume(G4VPhysicalVolume*);

  std::vector<G4VPhysicalVolume*> GetContents() const;

  inline G4int GetMaxContents() const { return fMaxContents; };

 protected:
  void Split();

  const OctreeNode* GetChildFromPosition(
    G4ThreeVector const&) const;

  OctreeNode* GetChildFromPosition(G4ThreeVector const& pos);

 private:
  G4ThreeVector fPosition, fHalfLengths;
  G4int fMaxContents;

  std::vector<G4VPhysicalVolume*> fContents;

  OctreeNode* fParent;
  // fChildren is arranged logically to save on queries
  // The scheme is defined by quadrant as follows:
  // X Y Z Index | X Y Z Index
  // + + +   0   | - + +   4
  // + + -   1   | - + -   5
  // + - +   2   | - - +   6
  // + - -   3   | - - -   7
  std::array<OctreeNode*, 8> fChildren;
  G4double fHalfLengthsMag;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif  // MOLECULAR_OCTREE_NODE_HH
