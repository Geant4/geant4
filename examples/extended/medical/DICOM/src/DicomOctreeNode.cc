
//#include <stddef.h>
#include "DicomOctreeNode.hh"

G4int OctreeNode::mInstanceCounter = 0;

OctreeNode::OctreeNode()
{
  mParent = 0;
  mInstanceCounter++;
}

OctreeNode::~OctreeNode()
{
  mInstanceCounter--;
}

OctreeNode::OctreeNode( const OctreeNode* pParent )
{
  mParent = pParent;
  mInstanceCounter++;
}

