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

// Author: Paul Kent (CERN), 12.07.1995 - Initial version
// --------------------------------------------------------------------
#ifndef G4SMARTVOXELPROXY_HH
#define G4SMARTVOXELPROXY_HH

#include <assert.h>

#include "geomwdefs.hh"
#include "G4Types.hh"

class G4SmartVoxelNode;
class G4SmartVoxelHeader;

/**
 * @brief G4SmartVoxelProxy is a class for proxying smart voxels. The class
 * represents either a header (in turn referring to more VoxelProxies) or a
 * node. If created as a node, calls to GetHeader() cause an exception, and
 * likewise GetNode() when a header.
 */

class G4SmartVoxelProxy 
{
  public:

    /**
     * Constructs a Proxy for the specified header.
     *  @param[in] pHeader Pointer to the voxel header.
     */
    inline G4SmartVoxelProxy(G4SmartVoxelHeader* pHeader);

    /**
     * Constructs a Proxy for the specified node.
     *  @param[in] pNode Pointer to the voxel node.
     */
    inline G4SmartVoxelProxy(G4SmartVoxelNode* pNode);

    /**
     * Default destructor. Not responsible for proxied objects.
     */
    ~G4SmartVoxelProxy() = default;

    /**
     * Equality operator. True when objects share the same address.
     */
    inline G4bool operator == (const G4SmartVoxelProxy& v) const;

    /**
     * Returns true if proxying for a header, else false.
     */
    inline G4bool IsHeader() const;

    /**
     * Returns true if proxying for a node, else false.
     */
    inline G4bool IsNode() const;

    /**
     * Returns the pointer to the proxied node, else throws a G4Exception.
     */
    inline G4SmartVoxelNode* GetNode() const;
 
    /**
     * Returns the pointer to the proxied header, else throws a G4Exception.
     */
    inline G4SmartVoxelHeader* GetHeader() const;

  private:

    G4SmartVoxelHeader* fHeader = nullptr;
    G4SmartVoxelNode* fNode = nullptr;
};

#include "G4SmartVoxelProxy.icc"

#endif
