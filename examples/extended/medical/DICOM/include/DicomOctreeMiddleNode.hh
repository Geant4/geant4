#ifndef DicomOctreeMiddleNode_h
#define DicomOctreeMiddleNode_h
#include "globals.hh"
class DicomOctree;

class DicomOctreeMiddleNode : public DicomOctreeNode
{
public:
  DicomOctreeMiddleNode();
  ~DicomOctreeMiddleNode();

public:
  void ResetFamily();
  DicomOctreeMiddleNode( DicomOctreeNode* pParent );
  G4int FindChild( const DicomOctreeNode* pNode );
  G4int MemSize();
  DicomOctreeNode*& operator []( G4int index )   {return mChildren[index];}
  OctreeNodeType Type()                     {return MIDDLE_NODE;}

private:
  DicomOctreeNode*  mChildren[8];
  DicomOctree* octree;
};
#endif
