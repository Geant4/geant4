// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SmartVoxelProxy.hh,v 1.1 1999-01-07 16:07:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// class G4SmartVoxelProxy
//
// Class for proxying smart voxels. The class
// represents either a header (in turn refering to more VoxelProxies)
// or a node. If created as a node, calls to GetHeader cause an exception,
// and likewise GetNode when a header.
//
// Note that the proxy does NOT gain deletion responsibility for proxied
// objects.
//
// Member functions:
//
// G4SmartVoxelProxy(G4SmartVoxelHeader *pHeader);
//   Proxy for the specified header
// G4SmartVoxelProxy(G4SmartVoxelNode *pNode)
//   Proxy for the specified node
// G4bool IsHeader() const
//   Return true if proxying for a header, else false
// G4bool IsNode() const
//   Return true if proxying for a node, else false
// G4SmartVoxelNode* GetNode() const
//   Return ptr to proxied node, else call G4Exception
// G4SmartVoxelHeader* GetHeader() const
//   Return ptr to proxied header, else call G4Exception
//
//
// operator == (const G4SmartVoxelProxy& v)
//   True when objects share same address.
//
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

public:
    G4SmartVoxelProxy(G4SmartVoxelHeader *pHeader)
    {
	fHeader=pHeader;
	fNode=0;
    }

    G4SmartVoxelProxy(G4SmartVoxelNode *pNode)
    {
	fHeader=0;
	fNode=pNode;
    }

// Destructor - do nothing. Not responsible for proxied objects
    ~G4SmartVoxelProxy() {;}

    G4bool IsHeader() const
    {
	return (fHeader) ? true:false;
    }   

    G4bool IsNode() const
    {
	return (fNode) ? true:false;
    }
    
    G4SmartVoxelNode* GetNode() const
    {
	assert(fNode != 0);
	return fNode;
    }

    G4SmartVoxelHeader* GetHeader() const
    {
	assert(fHeader != 0);
	return fHeader;
    }

    G4bool operator == (const G4SmartVoxelProxy& v) const
    {
	return (this==&v) ? true : false;
    }

private:
    G4SmartVoxelNode* fNode;
    G4SmartVoxelHeader* fHeader;
};

#endif



