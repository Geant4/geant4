#ifndef DicomOctreeMiddleNode_h
#define DicomOctreeMiddleNode_h
#include "globals.hh"

class Octree;

class MiddleNode : public OctreeNode
{
public:
  MiddleNode();
  ~MiddleNode();

public:
  void ResetFamily();
  MiddleNode( OctreeNode* pParent );
  G4int FindChild( const OctreeNode* pNode );
  G4int MemSize();
  OctreeNode*& operator []( G4int index )   {return mChildren[index];}
  OctreeNodeType Type()                     {return MIDDLE_NODE;}

private:
 OctreeNode*  mChildren[8];
};
#endif
