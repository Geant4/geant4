#ifndef DicomOctreeMiddleNode_h
#define DicomOctreeMiddleNode_h
#include "globals.hh"

class Octree;

class MiddleNode : public OctreeNode
{
  //    friend class Octree;

public:
  OctreeNode*  mChildren[8];

  void ResetFamily();
  MiddleNode();
  MiddleNode( OctreeNode* pParent );
  ~MiddleNode();

  OctreeNode*& operator []( G4int index );
  OctreeNodeType Type();
  G4int FindChild( const OctreeNode* pNode );
  G4int MemSize();
};
#endif
