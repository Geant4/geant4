#include "DicomOctreeNode.hh"
#include "DicomOctreeTerminalNode.hh"


OctreeNode* DicomOctreeTerminalNode::mNull = 0;
DicomOctreeTerminalNode::DicomOctreeTerminalNode( OctreeNode* pParent) : OctreeNode( pParent )
{

}
DicomOctreeTerminalNode::~DicomOctreeTerminalNode()
{

}

OctreeNode*& DicomOctreeTerminalNode::operator []( G4int )
{
  return mNull;
}

G4int DicomOctreeTerminalNode::MemSize()
{
  return sizeof(DicomOctreeTerminalNode);
}


