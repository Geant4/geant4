// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GeometryManager.hh,v 1.2 1999-12-15 14:49:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4GeometryManager
//
// A class resposible for high level geometrical functions, and for
// high level objects in the geometry subdomain.
// The class is `singleton', with access via G4GeometryManager::GetInstance
//
// Member functions:
//
// G4bool CloseGeometry(G4bool pOptimise=true);
//   Close (`lock') the geometry: perform sanity and `completion' checks
//   and optionally [default=yes] Build optimisation information.
//  
// void OpenGeometry();
//   Open (`unlock') the geometry and remove optimisation information if
//   present.
//
// static G4GeometryManager* GetInstance()
//   Return ptr to singleton instance of the class.
//
// Member data:
//
// static G4GeometryManager* fgInstance
//   Ptr to the unique instance of class
//
// History:
// 26.07.95 P.Kent Initial version, incuding optimisation Build

#ifndef G4GEOMETRYMANAGER_HH
#define G4GEOMETRYMANAGER_HH

#include "globals.hh"

// Needed for building optimisations
#include "geomdefs.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4SmartVoxelHeader.hh"
#ifdef  G4GEOMETRY_VOXELDEBUG
#include "G4ios.hh"
#endif

class G4GeometryManager
{
public:
    G4bool CloseGeometry(G4bool pOptimise=true);
    void OpenGeometry();
    static G4GeometryManager* GetInstance();

protected:
    G4GeometryManager();
private:
    void BuildOptimisations(const G4bool allOpt);
    void DeleteOptimisations();

    static G4GeometryManager* fgInstance;
    G4bool fIsClosed;
};

#endif


