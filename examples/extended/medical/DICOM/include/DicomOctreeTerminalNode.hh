#ifndef DicomOctreeTerminalNode_h
#define DicomOctreeTerminalNode_h

class Octree;
class TerminalNode : public OctreeNode
{
public:

  TerminalNode(const OctreeNode* pParent );
  ~TerminalNode();

public:

  // ---- MGP ---- Replace *& with proper design of the operator
  OctreeNode*& operator []( G4int index );

  G4int MemSize();

  OctreeNodeType Type()                      {return TERMINAL_NODE;}
  G4int FindChild( const OctreeNode* pNode ) {return -1;}

private:
  static OctreeNode* mNull;
};
#endif
