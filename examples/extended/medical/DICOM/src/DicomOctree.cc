//#include "common.h"
//#pragma hdrstop
//---------------------------------------------------------------------------
#include <stddef.h>
//---------------------------------------------------------------------------
#include "DicomOctree.hh"
//---------------------------------------------------------------------------
Octree::Octree( int noLevels, float size )
{
    mNoLevels = noLevels;
    mSize = size;
    mRoot = new MiddleNode( NULL );
}
//---------------------------------------------------------------------------
Octree::~Octree()
{
    if ( mRoot != NULL ) delete mRoot; 
}
//---------------------------------------------------------------------------
OctreeNode* Octree::CreateNode( float nodeX, float nodeY, float nodeZ, int level )
{
    OctreeNode* current = mRoot;
    float curr_x = 0;
    float curr_y = 0;
    float curr_z = 0;
    // Make children
    for ( int i = 0; i < 8; i++ )
    {
        float child_level_resolution = (1 << (i+1) );
        float child_size = mSize / child_level_resolution;

        int dir_x = int( nodeX >= curr_x + child_size );
        int dir_y = int( nodeY >= curr_y + child_size );
        int dir_z = int( nodeZ >= curr_z + child_size );
        int direction = dir_x + ( dir_y << 1 ) + ( dir_z << 2 );

        if ( (*current)[direction] == 0 ) // NULL
        {
            if ( i < level - 1 )
            {
                (*current)[direction] = new MiddleNode( current );
            } else
            {
                (*current)[direction] = new TerminalNode( current );
            }
        }

        current = (*current)[direction];
        curr_x += dir_x*child_size;
        curr_y += dir_y*child_size;
        curr_z += dir_z*child_size;
    }
    return current;
}
//---------------------------------------------------------------------------
OctreeNode* Octree::operator()( float nodeX, float nodeY, float nodeZ, int level )
{
    OctreeNode* current = mRoot;
    float curr_x = 0;
    float curr_y = 0;
    float curr_z = 0;
    // Make children
    for ( int i = 0; i < level; i++ )
    {
        float child_level_resolution = ( 1 << (i+1) );
        float child_size = mSize / child_level_resolution;

        int dir_x = int( nodeX >= curr_x + child_size );
        int dir_y = int( nodeY >= curr_y + child_size );
        int dir_z = int( nodeZ >= curr_z + child_size );
        int direction = dir_x + ( dir_y << 1 ) + ( dir_z << 2 );

        if ( (*current)[direction] == 0/*NULL*/ ) return 0/*NULL*/;

        current = (*current)[direction];
        curr_x += dir_x*child_size;
        curr_y += dir_y*child_size;
        curr_z += dir_z*child_size;
    }
    return current;
}
//---------------------------------------------------------------------------
float Octree::Size()
{
    return mSize;
}
//---------------------------------------------------------------------------
OctreeNode* Octree::Root()
{
    return mRoot;
}
//---------------------------------------------------------------------------
int Octree::NoLevels()
{
    return mNoLevels;
}
//---------------------------------------------------------------------------
int Octree::Resolution()
{
    return ( 1 << mNoLevels );
}
//---------------------------------------------------------------------------
void Octree::DeleteTree()
{
    delete mRoot;
    mRoot = NULL;
}
//---------------------------------------------------------------------------
void Octree::CountRecursive( OctreeNode* pNode, int &rMiddle, int &rTerminal )
{
    if ( pNode->Type() == MIDDLE_NODE )
    {
        rMiddle++;
    } else
    {
        rTerminal++;
    }
    for ( unsigned int i = 0; i < 8; i++ )
    {
        if ( (*pNode)[i] != NULL )
        {
            CountRecursive( (*pNode)[i], rMiddle, rTerminal );
        }
    }
}
//---------------------------------------------------------------------------
int Octree::CountMemory( int &rMiddle, int &rTerminal )
{
    CountRecursive( mRoot, rMiddle, rTerminal );
    int total = rMiddle*sizeof(MiddleNode) + rTerminal*sizeof(TerminalNode);
    return total;
}
//---------------------------------------------------------------------------

