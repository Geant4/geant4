#ifndef OCTREE_NODE_H
#define OCTREE_NODE_H
#include "globals.hh"

class Octree;

enum OctreeNodeType { MIDDLE_NODE, TERMINAL_NODE };

//  OctreeNode represents a single node in the octree.

class OctreeNode
{
public:

  OctreeNode();
  ~OctreeNode();
  OctreeNode( OctreeNode* pParent );
  
  G4double Density(){return mDensity;}

  virtual OctreeNode*& operator []( G4int index ) = 0;

  virtual OctreeNodeType Type() = 0;

  const OctreeNode* Parent(){ return mParent; }

  virtual G4int FindChild( const OctreeNode* pNode ) = 0;

  static G4int InstanceCounter() { return mInstanceCounter; }

  virtual G4int MemSize() = 0;

private:
  static G4int mInstanceCounter;
  OctreeNode* mParent;
  G4double  mDensity;
};
#endif 
