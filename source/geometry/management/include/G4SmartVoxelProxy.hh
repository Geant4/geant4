// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SmartVoxelProxy.hh,v 1.3 2000-04-20 16:49:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4SmartVoxelProxy
//
// Class description:
//
// Class for proxying smart voxels. The class represents either a header
// (in turn refering to more VoxelProxies) or a node. If created as a node,
// calls to GetHeader cause an exception, and likewise GetNode when a header.
//
// Note that the proxy does NOT gain deletion responsibility for proxied
// objects.

// History:
// 12.07.95 P.Kent Initial version
// 03.08.95 P.Kent Updated to become non abstract class, removing
//                 HeaderProxy and NodeProxy derived classes

#ifndef G4SMARTVOXELPROXY_HH
#define G4SMARTVOXELPROXY_HH

#include "globals.hh"
#include <assert.h>

class G4SmartVoxelNode;
class G4SmartVoxelHeader;

class G4SmartVoxelProxy 
{

  public:  // with description

    G4SmartVoxelProxy(G4SmartVoxelHeader *pHeader)
      : fHeader(pHeader), fNode(0) {}
      // Proxy for the specified header.

    G4SmartVoxelProxy(G4SmartVoxelNode *pNode)
      : fHeader(0), fNode(pNode) {}
      // Proxy for the specified node.

    ~G4SmartVoxelProxy() {}
      // Destructor - do nothing. Not responsible for proxied objects.

    G4bool IsHeader() const;
      // Return true if proxying for a header, else false.

    G4bool IsNode() const;
      // Return true if proxying for a node, else false.

    G4SmartVoxelNode* GetNode() const;
      // Return ptr to proxied node, else call G4Exception.

    G4SmartVoxelHeader* GetHeader() const;
      // Return ptr to proxied header, else call G4Exception

    G4bool operator == (const G4SmartVoxelProxy& v) const;
      // Equality operator.
      // True when objects share same address.

  private:

    G4SmartVoxelNode* fNode;
    G4SmartVoxelHeader* fHeader;
};

#include "G4SmartVoxelProxy.icc"

#endif
