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
//
// $Id: G4GeometryManager.hh 103235 2017-03-22 15:53:48Z gcosmo $
//
// class G4GeometryManager
//
// Class description:
//
// A class responsible for high level geometrical functions, and for
// high level objects in the geometry subdomain.
// The class is a `singleton', with access via the static method
// G4GeometryManager::GetInstance().
//
// Member data:
//
//   - fgInstance
//     Ptr to the unique instance of class

// Author:
// 26.07.95 P.Kent Initial version, including optimisation Build
// --------------------------------------------------------------------
#ifndef G4GEOMETRYMANAGER_HH
#define G4GEOMETRYMANAGER_HH

#include <vector>
#include "G4Types.hh"
#include "G4SmartVoxelStat.hh"

class G4VPhysicalVolume;

class G4GeometryManager
{
  public: // with description
  
    G4bool CloseGeometry(G4bool pOptimise=true, G4bool verbose=false,
                         G4VPhysicalVolume* vol=0);
      // Close (`lock') the geometry: perform sanity and `completion' checks
      // and optionally [default=yes] build optimisation information.
      // Applies to just a specific subtree if a physical volume is specified.

    void OpenGeometry(G4VPhysicalVolume* vol=0);
      // Open (`unlock') the geometry and remove optimisation information if
      // present. Applies to just a specific subtree if a physical volume is
      // specified.

    static G4bool IsGeometryClosed();
      // Return true/false according to state of optimised geoemtry.

    void SetWorldMaximumExtent(G4double worldExtent);
      // Set the maximum extent of the world volume. The operation is
      // allowed only if NO solids have been created already.

    static G4GeometryManager* GetInstance();
      // Return ptr to singleton instance of the class, creating it if
      // not existing.

    static G4GeometryManager* GetInstanceIfExist();
      // Return ptr to singleton instance.

  public:  // without description

   ~G4GeometryManager();
      // Destructor.

  protected:

    G4GeometryManager();
      // Protected constructor

  private:

    void BuildOptimisations(G4bool allOpt, G4bool verbose=false);
    void BuildOptimisations(G4bool allOpt, G4VPhysicalVolume* vol);
    void DeleteOptimisations();
    void DeleteOptimisations(G4VPhysicalVolume* vol);
    static void ReportVoxelStats( std::vector<G4SmartVoxelStat> & stats,
                                  G4double totalCpuTime );
    static G4ThreadLocal G4GeometryManager* fgInstance;
    static G4ThreadLocal G4bool fIsClosed;
};

#endif
