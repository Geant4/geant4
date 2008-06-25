#ifndef _G4STEPREAD_INCLUDED_
#define _G4STEPREAD_INCLUDED_

#include "G4TessellatedSolid.hh"
#include "G4QuadrangularFacet.hh"
#include "G4TriangularFacet.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"

#include <fstream>

class G4STEPRead {
   std::vector<G4LogicalVolume*> volumeList;
   std::vector<G4TessellatedSolid*> tessellatedList;

   G4Material* solid_material;
   G4Material* medium_material;

   void tessellatedRead(const std::string&);
   void facetRead(const std::string&);
public:
   void ReadGeom(const G4String&);
};

#endif
