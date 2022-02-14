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
// G4STRead
//
// Class description:
//
// GDML class for import of triangularised geometry descriptions
// (.geom and .tree structures) generated out of STEP files from
// CAD systems (STWriter STEP Tool and similar...).

// Author: Zoltan Torzsok, November 2007
// --------------------------------------------------------------------
#ifndef G4STREAD_HH
#define G4STREAD_HH 1

#include <vector>
#include <map>

#include "G4ThreeVector.hh"

class G4Material;
class G4Box;
class G4TessellatedSolid;
class G4LogicalVolume;

class G4STRead
{
  public:

    G4LogicalVolume* Read(const G4String&, G4Material* mediumMaterial,
                          G4Material* solidMaterial);

  private:

    void TessellatedRead(const std::string&);
    void FacetRead(const std::string&);
    void PhysvolRead(const std::string&);
    void ReadGeom(const G4String&);
    void ReadTree(const G4String&);

  private:

    G4Box* world_box = nullptr;
    G4ThreeVector world_extent;
    G4Material* solid_material = nullptr;
    G4LogicalVolume* world_volume = nullptr;

    std::vector<G4TessellatedSolid*> tessellatedList;
    std::map<G4TessellatedSolid*, G4LogicalVolume*> volumeMap;
};

#endif
