#include "DicomOctreeNode.hh"
#include "DicomOctreeTerminalNode.hh"


OctreeNode* TerminalNode::mNull = 0;
TerminalNode::TerminalNode( OctreeNode* pParent) : OctreeNode( pParent )
{

}
TerminalNode::~TerminalNode()
{

}

OctreeNode*& TerminalNode::operator []( G4int )
{
  return mNull;
}

G4int TerminalNode::MemSize()
{
  return sizeof(TerminalNode);
}


