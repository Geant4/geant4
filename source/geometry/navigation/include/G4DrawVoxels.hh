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
// $Id: G4DrawVoxels.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
// 
// class G4DrawVoxels
//
// Class description:
//
// Utility class for the visualization of voxels in the detector geometry.
// Define G4DrawVoxelsDebug in the environment at compilation for debugging
// information printed to G4cout.

// History:
// 03/08/1999 The G4VisAttributes have been made member data for lifetime
//            reasons / visualisation - L.G (ask John Allison for further
//            explanation).
// 29/07/1999 First comitted version - L.G.
// --------------------------------------------------------------------
#ifndef G4DrawVoxels_HH
#define G4DrawVoxels_HH

#include "G4VisAttributes.hh"
#include "G4VoxelLimits.hh"
#include "G4PlacedPolyhedron.hh"

class G4SmartVoxelHeader;
class G4LogicalVolume;

// ***********************************************************************

class G4DrawVoxels
{
  public: // with description

    G4DrawVoxels();
      // Constructor. It initialises the members data to default colors
      // Copy constructor and assignment operator not supported (array
      // fvoxelcolours ...).

    ~G4DrawVoxels();
      // Destructor NOT virtual. Not a base class.
    
    void DrawVoxels(const G4LogicalVolume* lv) const;
    G4PlacedPolyhedronList* CreatePlacedPolyhedra(const G4LogicalVolume*) const;

    void SetVoxelsVisAttributes(G4VisAttributes&,
                                G4VisAttributes&,
                                G4VisAttributes&);
    void SetBoundingBoxVisAttributes(G4VisAttributes&);

  private:
  
    // Member data
    G4VisAttributes fVoxelsVisAttributes[3];
    G4VisAttributes fBoundingBoxVisAttributes;
    
    void ComputeVoxelPolyhedra(const G4LogicalVolume*,
                               const G4SmartVoxelHeader*,
                                     G4VoxelLimits&,
                                     G4PlacedPolyhedronList*) const;
    
    // Copy constructor Assignment operator not supported
    // (array fvoxelcolours ...)
    G4DrawVoxels(const G4DrawVoxels&);	
    G4DrawVoxels operator=(const G4DrawVoxels&);	
};

#endif
