#ifndef _G4STEPREAD_INCLUDED_
#define _G4STEPREAD_INCLUDED_

#include "G4TessellatedSolid.hh"
#include "G4QuadrangularFacet.hh"
#include "G4TriangularFacet.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4Box.hh"

#include <fstream>
#include <vector>
#include <map>

class G4STEPRead {
   G4Box* world_box;
   G4ThreeVector world_extent;
   G4Material* solid_material;
   G4Material* medium_material;
   G4LogicalVolume* world_volume;
   std::vector<G4TessellatedSolid*> tessellatedList;
   std::map<G4TessellatedSolid*,G4LogicalVolume*> volumeMap;

   void tessellatedRead(const std::string&);
   void facetRead(const std::string&);
   void physvolRead(const std::string&);
   void ReadGeom(const G4String&);
   void ReadTree(const G4String&);
public:
   void Read(const G4String&,G4Material* mediumMaterial,G4Material* solidMaterial);
   G4VPhysicalVolume* GetWorldVolume();
};

#endif
