#ifndef __OCTREE_H__
#define __OCTREE_H__

#include "DicomOctreeNode.hh"

//---------------------------------------------------------------------------
/// \brief Octree encapsulates a volumetric compressed representation of the
/// vector field used for the reconstruction of a 3D model from multiple
/// range images or curves.
/// 
/// It is assumed that the space occupied by octree is a cube whose
/// bottom-left-front corner is at origin (0,0,0) while the oposite corner
/// (top,right,back) is at ( mSize,mSize,mSize). Resolution of an octree
/// denotes the number of voxels along each axis such
/// that the resolution = 2^mNoLevels.


class Octree
{
private:
  /// Root node of the tree
  MiddleNode* mRoot;
  
  /// In an octree, size denotes physical size of the cube,
  /// ie length of its sides (which are assumed equal).
  float      mSize;
  
  /// Number of levels denotes the maximal number of
  /// nodes in a single branch, starting from the root node.
  int        mNoLevels;

  void CountRecursive( OctreeNode* pNode, int &rMiddle, int &rTerminal );

public:
  
  Octree( int noLevels, float size );
  ~Octree();
  void DeleteTree();

  OctreeNode* CreateNode( float i, float j, float k, int level );
  OctreeNode* operator()( float nodeX, float nodeY, float nodeZ, int level );
    
  OctreeNode* Root();
  float Size();
  int NoLevels();
  int Resolution();
  int CountMemory( int &rMiddle, int &rTerminal );
};
#endif 

