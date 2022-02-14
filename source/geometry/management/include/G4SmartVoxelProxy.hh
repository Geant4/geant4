//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4SmartVoxelProxy
//
// Class description:
//
// Class for proxying smart voxels. The class represents either a header
// (in turn refering to more VoxelProxies) or a node. If created as a node,
// calls to GetHeader cause an exception, and likewise GetNode when a header.
//
// Note that the proxy does NOT gain deletion responsibility for proxied
// objects.

// 12.07.95, P.Kent - Initial version
// --------------------------------------------------------------------
#ifndef G4SMARTVOXELPROXY_HH
#define G4SMARTVOXELPROXY_HH 1

#include <assert.h>

#include "geomwdefs.hh"
#include "G4Types.hh"
#include "G4Allocator.hh"

class G4SmartVoxelNode;
class G4SmartVoxelHeader;

class G4SmartVoxelProxy 
{

  public:

    G4SmartVoxelProxy(G4SmartVoxelHeader *pHeader)
      : fHeader(pHeader) {}
      // Proxy for the specified header.

    G4SmartVoxelProxy(G4SmartVoxelNode *pNode)
      : fNode(pNode) {}
      // Proxy for the specified node.

    ~G4SmartVoxelProxy() = default;
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

    inline void* operator new(std::size_t);
      // Override "new" for "G4Allocator".
    inline void operator delete(void* aProxy);
      // Override "delete" for "G4Allocator".

  private:

    G4SmartVoxelHeader* fHeader = nullptr;
    G4SmartVoxelNode* fNode = nullptr;
};

#include "G4SmartVoxelProxy.icc"

#endif
