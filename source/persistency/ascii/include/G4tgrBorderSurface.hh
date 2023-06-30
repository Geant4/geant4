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
// G4tgrBorderSurface
//
// Class description:
// G4tgrBorderSurface is a class to register optical surface
// properties.
//
#ifndef G4tgrBorderSurface_hh
#define G4tgrBorderSurface_hh 1

#include "globals.hh"
#include "G4tgrMaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"

class G4tgrBorderSurface
{
  public:

    G4tgrBorderSurface(const std::vector<G4String>&);
    virtual ~G4tgrBorderSurface();

    G4String GetName() const { return theName; }
    G4String GetV1Name() const { return V1Name; }
    G4String GetV2Name() const { return V2Name; }
    G4int    GetCopyNo1() const { return copyNo1; }
    G4int    GetCopyNo2() const { return copyNo2; }
    G4OpticalSurface* GetOptic() const { return theOptic; }
    G4tgrMaterialPropertiesTable* 
      GetTgrMaterialPropertiesTable() const { return theMPT; }

    // Converts string to enum.
    G4SurfaceType          GetType(G4String);
    G4OpticalSurfaceFinish GetFinish(G4String);
    G4OpticalSurfaceModel  GetModel(G4String);

  protected:

    G4String theName; /// Name of surface
    G4String V1Name;  /// Name of volume 1
    G4String V2Name;  /// Name of volume 2
    G4int copyNo1;    /// Copy number of volume 1
    G4int copyNo2;    /// Copy number of volume 2

    /// Points to G4OpticalSurface object
    G4OpticalSurface* theOptic; 
    /// Border surface properties
    G4tgrMaterialPropertiesTable* theMPT = nullptr; 

};

#endif
