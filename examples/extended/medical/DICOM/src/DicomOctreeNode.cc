
#include <stddef.h>
#include "DicomOctreeNode.hh"

G4int DicomOctreeNode::mInstanceCounter = 0;

DicomOctreeNode::DicomOctreeNode()
{
  mParent = 0;
  mInstanceCounter++;
}

DicomOctreeNode::~DicomOctreeNode()
{
  mInstanceCounter--;
}

DicomOctreeNode::DicomOctreeNode( DicomOctreeNode* pParent )
{
  mParent = pParent;
  mInstanceCounter++;
}

