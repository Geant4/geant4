#ifndef DICOMOCTREENODE_H
#define DICOMOCTREENODE_H
#include "globals.hh"

enum OctreeNodeType { MIDDLE_NODE, TERMINAL_NODE };

//  OctreeNode represents a single node in the Octree.
class DicomOctreeNode
{
public:

  DicomOctreeNode();
  ~DicomOctreeNode();

  DicomOctreeNode( DicomOctreeNode* pParent );
  
  G4double Density(){return mDensity;}

  virtual DicomOctreeNode*& operator []( G4int index ) = 0;

  virtual OctreeNodeType Type() = 0;

  const DicomOctreeNode* Parent(){ return mParent; }

  virtual G4int FindChild( const DicomOctreeNode* pNode ) = 0;

  static G4int InstanceCounter() { return mInstanceCounter; }

  virtual G4int MemSize() = 0;

private:
  static G4int mInstanceCounter;
  DicomOctreeNode* mParent;
  G4double  mDensity;
};
#endif 
