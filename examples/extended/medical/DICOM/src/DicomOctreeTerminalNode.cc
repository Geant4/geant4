#include "DicomOctreeNode.hh"
#include "DicomOctreeTerminalNode.hh"

DicomOctreeNode* DicomOctreeTerminalNode::mNull = 0;

DicomOctreeTerminalNode::DicomOctreeTerminalNode( DicomOctreeNode* pParent) : DicomOctreeNode( pParent )
{

}
DicomOctreeTerminalNode::~DicomOctreeTerminalNode()
{

}

DicomOctreeNode*& DicomOctreeTerminalNode::operator []( G4int )
{
  return mNull;
}

G4int DicomOctreeTerminalNode::MemSize()
{
  return sizeof(DicomOctreeTerminalNode);
}


