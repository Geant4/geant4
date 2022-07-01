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
#include "OctreeNode.hh"
#include "G4VPhysicalVolume.hh"
#include <cassert>
#include <utility>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OctreeNode::OctreeNode(G4ThreeVector position, G4ThreeVector halfLengths,
                       G4int maxContents, OctreeNode* parent)
  : fPosition(std::move(position))
  , fHalfLengths(std::move(halfLengths))
  , fMaxContents(maxContents)
  , fParent(parent)
  , fChildren(std::array<OctreeNode*, 8>())
{
  fHalfLengthsMag = fHalfLengths.mag();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OctreeNode::~OctreeNode()
{
  // Delete children
  for(auto& it : fChildren)
  {
    delete it;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const OctreeNode* OctreeNode::GetChildFromPosition(
  const G4ThreeVector& pos) const
{
  G4ThreeVector localPosition = pos - fPosition;

  // Get the right child
  // The scheme is defined by quadrant as follows:
  // Zero is treated as positive
  // X Y Z Index | X Y Z Index
  // + + +   0   | - + +   4
  // + + -   1   | - + -   5
  // + - +   2   | - - +   6
  // + - -   3   | - - -   7
  int bit2 = (localPosition.getX() < 0) ? 4 : 0;
  int bit1 = (localPosition.getY() < 0) ? 2 : 0;
  int bit0 = (localPosition.getZ() < 0) ? 1 : 0;
  return fChildren[bit0 + bit1 + bit2];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Non-const version of above function. Ugly but useful in this case
OctreeNode* OctreeNode::GetChildFromPosition(const G4ThreeVector& pos)
{
  G4ThreeVector localPosition = pos - fPosition;

  // Get the right child
  // The scheme is defined by quadrant as follows:
  // Zero is treated as positive
  // X Y Z Index | X Y Z Index
  // + + +   0   | - + +   4
  // + + -   1   | - + -   5
  // + - +   2   | - - +   6
  // + - -   3   | - - -   7
  int bit2 = (localPosition.getX() < 0) ? 4 : 0;
  int bit1 = (localPosition.getY() < 0) ? 2 : 0;
  int bit0 = (localPosition.getZ() < 0) ? 1 : 0;
  int idx  = bit0 + bit1 + bit2;

  assert(idx < (int) fChildren.size());
  return fChildren[idx];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OctreeNode::AddPhysicalVolume(G4VPhysicalVolume* pv)
{
  // Consider adding test for the validity of the PV position
  if(this->HasChildren())
  {
    OctreeNode* child = this->GetChildFromPosition(pv->GetTranslation());
    child->AddPhysicalVolume(pv);
  }
  else
  {
    // if there are fMaxContents elements in the bin, we need to split
    // the bin before adding a new element.
    if((G4int) fContents.size() == fMaxContents)
    {
      this->Split();
      this->AddPhysicalVolume(pv);
    }
    else
    {
      fContents.push_back(pv);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OctreeNode::Split()
{
  G4ThreeVector hl = 0.5 * fHalfLengths;
  G4double xp      = fPosition.getX();
  G4double yp      = fPosition.getY();
  G4double zp      = fPosition.getZ();
  G4double xl      = hl.getX();
  G4double yl      = hl.getY();
  G4double zl      = hl.getZ();
  G4ThreeVector newpos;

  // Construct children
  newpos       = G4ThreeVector(xp + xl, yp + yl, zp + zl);
  fChildren[0] = new OctreeNode(newpos, hl, fMaxContents, this);
  newpos       = G4ThreeVector(xp + xl, yp + yl, zp - zl);
  fChildren[1] = new OctreeNode(newpos, hl, fMaxContents, this);
  newpos       = G4ThreeVector(xp + xl, yp - yl, zp + zl);
  fChildren[2] = new OctreeNode(newpos, hl, fMaxContents, this);
  newpos       = G4ThreeVector(xp + xl, yp - yl, zp - zl);
  fChildren[3] = new OctreeNode(newpos, hl, fMaxContents, this);
  newpos       = G4ThreeVector(xp - xl, yp + yl, zp + zl);
  fChildren[4] = new OctreeNode(newpos, hl, fMaxContents, this);
  newpos       = G4ThreeVector(xp - xl, yp + yl, zp - zl);
  fChildren[5] = new OctreeNode(newpos, hl, fMaxContents, this);
  newpos       = G4ThreeVector(xp - xl, yp - yl, zp + zl);
  fChildren[6] = new OctreeNode(newpos, hl, fMaxContents, this);
  newpos       = G4ThreeVector(xp - xl, yp - yl, zp - zl);
  fChildren[7] = new OctreeNode(newpos, hl, fMaxContents, this);

  // Distribute contents to children

  for(int i = fContents.size() - 1; i >= 0; --i)
  {
    G4VPhysicalVolume* pv = fContents[i];
    OctreeNode* child     = this->GetChildFromPosition(pv->GetTranslation());
    assert(child != this);
    child->AddPhysicalVolume(pv);
  }
  fContents.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4VPhysicalVolume*> OctreeNode::GetContents() const
{
  if(this->HasChildren())  // if has sub-nodes
  {
    std::vector<G4VPhysicalVolume*> vec;
    std::vector<G4VPhysicalVolume*> childCont;
    for(auto it = fChildren.begin(); it != fChildren.end(); ++it)
    {
      childCont = (*it)->GetContents();
      for(auto jt = childCont.begin(); jt != childCont.end(); ++jt)
      {
        vec.push_back(*jt);
      }
    }
    return vec;
  }
  else  // if no sub-nodes
  {
    return fContents;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Search for Octree nodes in a volume
const std::vector<G4VPhysicalVolume*> OctreeNode::SearchOctree(
  const G4ThreeVector& pos, G4double _rad) const
{
  // Need to search based on absolute position of each volume rather than
  // relative  position, and maintain a list of candidate nodes and
  // final nodes within the radius

  std::vector<const OctreeNode*> nodes;
  std::vector<const OctreeNode*> candidates;
  nodes.reserve(512);
  candidates.reserve(512);

  if(this->HasChildren())
  {
    for(auto it = this->GetChildren().begin(); it != this->GetChildren().end();
        it++)
    {
      candidates.push_back(*it);
    }
  }
  else
  {
    candidates.push_back(this);
  }

  const OctreeNode* aNode;
  G4double d;
  while(!candidates.empty())
  {
    aNode = candidates.back();
    candidates.pop_back();
    d = (aNode->GetPosition() - pos).mag() - aNode->GetHalfLengthsMag();
    // if node within circle
    if(d < _rad)
    {
      if(aNode->HasChildren())
      {
        for(auto it = aNode->GetChildren().begin();
            it != aNode->GetChildren().end(); ++it)
        {
          candidates.push_back(*it);
        }
      }
      else
      {
        nodes.push_back(aNode);
      }
    }
  }

  // G4cout << "Found " << nodes.size() << " nodes" << G4endl;

  // Get the physical volumes
  G4double r2 = _rad * _rad;
  std::vector<G4VPhysicalVolume*> pvols;
  for(auto it = nodes.begin(); it != nodes.end(); ++it)
  {
    std::vector<G4VPhysicalVolume*> conts = (*it)->GetContents();
    for(auto jt = conts.begin(); jt != conts.end(); ++jt)
    {
      if((pos - (*jt)->GetTranslation()).mag2() < r2) {
        pvols.push_back((*jt));}
    }
  }
  // G4cout << "Found " << pvols.size() << " candidates" << G4endl;
  return pvols;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Search for Octree nodes in a volume
void OctreeNode::SearchOctree(const G4ThreeVector& pos,
                              std::vector<G4VPhysicalVolume*>& out,
                              G4double _rad) const
{
  // Need to search based on absolute position of each volume rather than
  // relative  position, and maintain a list of candidate nodes and
  // final nodes within the radius

  std::vector<const OctreeNode*> nodes;
  std::vector<const OctreeNode*> candidates;
  nodes.reserve(512);
  candidates.reserve(512);

  if(this->HasChildren())
  {
    for(auto it = this->GetChildren().begin(); it != this->GetChildren().end();
        ++it)
    {
      candidates.push_back(*it);
    }
  }
  else
  {
    candidates.push_back(this);
  }

  const OctreeNode* aNode;
  G4double d;
  while(!candidates.empty())
  {
    aNode = candidates.back();
    candidates.pop_back();
    d = (aNode->GetPosition() - pos).mag() - aNode->GetHalfLengthsMag();
    // if node within circle
    if(d < _rad)
    {
      if(aNode->HasChildren())
      {
        for(auto it = aNode->GetChildren().begin();
            it != aNode->GetChildren().end(); it++)
        {
          candidates.push_back(*it);
        }
      }
      else
      {
        nodes.push_back(aNode);
      }
    }
  }

  // Get the physical volumes
  G4double r2 = _rad * _rad;
  for(auto it = nodes.begin(); it != nodes.end(); ++it)
  {
    std::vector<G4VPhysicalVolume*> conts = (*it)->GetContents();
    for(auto jt = conts.begin(); jt != conts.end(); ++jt)
    {
      if((pos - (*jt)->GetTranslation()).mag2() < r2) {
        out.push_back((*jt));}
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::vector<G4VPhysicalVolume*> OctreeNode::SearchOctree(
  const G4ThreeVector& pos) const
{
  const OctreeNode* child = GetChildFromPosition(pos);
  return child->GetContents();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int OctreeNode::GetNumberOfTerminalNodes()
{
  if(!this->HasChildren()) {
    return 1;
}
  // sum the number of nodes of each child
  G4int n = 0;
  for(auto it = this->GetChildren().begin(); it != this->GetChildren().end();
      ++it)
  {
    n = n + (*it)->GetNumberOfTerminalNodes();
  }

  return n;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
