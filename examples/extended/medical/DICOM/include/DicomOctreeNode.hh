#ifndef __OCTREE_NODE_H__
#define __OCTREE_NODE_H__
//---------------------------------------------------------------------------
class Octree;
//---------------------------------------------------------------------------
enum OctreeNodeType { MIDDLE_NODE, TERMINAL_NODE };
//---------------------------------------------------------------------------
/// \brief OctreeNode represents a single node in the octree.
///
/// Longer description.
class OctreeNode
{
    friend class Octree;
private:
    static int       mInstanceCounter;
    OctreeNode*      mParent;
    float            mDensity;
public:
    OctreeNode();
    OctreeNode( OctreeNode* pParent );
    virtual ~OctreeNode();

    float&        Density();

    virtual OctreeNode*& operator []( int index ) = 0;
    virtual OctreeNodeType Type() = 0;
    OctreeNode* Parent();
    virtual int FindChild( OctreeNode* pNode ) = 0;

    static int InstanceCounter();

    virtual int MemSize() = 0;

};

//---------------------------------------------------------------------------
class MiddleNode : public OctreeNode
{
    friend class Octree;
private:

public:
    OctreeNode*  mChildren[8];

    void ResetFamily();
    MiddleNode();
    MiddleNode( OctreeNode* pParent );
    ~MiddleNode();

    OctreeNode*& operator []( int index );
    OctreeNodeType Type();
    int FindChild( OctreeNode* pNode );
    int MemSize();
};
//---------------------------------------------------------------------------
class TerminalNode : public OctreeNode
{
    friend class Octree;
private:
    static OctreeNode* mNull;
public:
    OctreeNode*& operator []( int index );
    OctreeNodeType Type();
    TerminalNode( OctreeNode* pParent );
    ~TerminalNode();
    int FindChild( OctreeNode* pNode );
    int MemSize();
};
//---------------------------------------------------------------------------

#endif // __OCTREE_NODE_H__
